"""Arithmetic expressions shared by result calculation and preprocessing."""

from __future__ import annotations

import ast
import keyword
from collections.abc import Mapping
from dataclasses import dataclass, field
from math import isfinite
from types import MappingProxyType

import numpy as np
from numpy.typing import ArrayLike, NDArray

from elisaflow.configuration.validation import finite_number, validate_config

NumericResult = float | NDArray[np.float64]
_OPERATORS = {ast.Add: np.add, ast.Sub: np.subtract, ast.Mult: np.multiply, ast.Div: np.divide}


@dataclass(frozen=True)
class Formula:
    """An arithmetic expression using ``sample`` and named control means.

    References always resolve to original control means supplied by the caller.
    Expressions cannot call functions, access attributes, or execute Python code.
    """

    expression: str
    references: Mapping[str, str] = field(default_factory=dict)
    _tree: ast.expr = field(init=False, repr=False, compare=False)
    _names: frozenset[str] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.expression, str) or not self.expression.strip() or len(self.expression) > 2048:
            raise ValueError("Een formule moet tussen 1 en 2048 tekens bevatten.")
        if not isinstance(self.references, Mapping):
            raise ValueError("Formulereferenties moeten referentienamen aan controlenamen koppelen.")
        for alias, control in self.references.items():
            if not isinstance(alias, str) or not alias.isidentifier() or keyword.iskeyword(alias) or alias in {"sample", "x"} or alias.startswith("_"):
                raise ValueError(f"Ongeldige formulereferentienaam: {alias!r}.")
            if not isinstance(control, str) or not control.strip():
                raise ValueError("Controles in formulereferenties moeten een niet-lege naam hebben.")
        try:
            tree = ast.parse(self.expression.strip(), mode="eval").body
        except (SyntaxError, ValueError, RecursionError) as error:
            raise ValueError("Ongeldige rekenkundige formule.") from error
        nodes = list(ast.walk(tree))
        if len(nodes) > 128:
            raise ValueError("De formule is te complex; gebruik minder bewerkingen (maximaal 128 syntaxiselementen).")
        allowed = (ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Load, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.UAdd, ast.USub)
        for node in nodes:
            if not isinstance(node, allowed):
                raise ValueError("Formules ondersteunen alleen getallen, referenties, haakjes en +, -, *, /.")
            if isinstance(node, ast.Constant):
                finite_number(node.value, "expression constant")
        names = frozenset(node.id for node in nodes if isinstance(node, ast.Name))
        unknown = names - {"sample", *self.references}
        if unknown:
            raise ValueError(f"Onbekende formulereferenties: {', '.join(sorted(unknown))}.")
        object.__setattr__(self, "references", MappingProxyType(dict(self.references)))
        object.__setattr__(self, "_tree", tree)
        object.__setattr__(self, "_names", names)

    @property
    def control_names(self) -> frozenset[str]:
        return frozenset(self.references[name] for name in self._names if name != "sample")

    @property
    def used_references(self) -> dict[str, str]:
        return {alias: control for alias, control in self.references.items() if alias in self._names}

    @classmethod
    def from_dict(cls, config: object) -> Formula:
        values = validate_config(config, {"expression", "references"}, "Formula")
        return cls(values.get("expression"), values.get("references", {}))

    def to_dict(self) -> dict:
        return {"expression": self.expression, **({"references": dict(self.references)} if self.references else {})}

    def affine_coefficients(self, controls: Mapping[str, float]) -> tuple[float, float] | None:
        """Return slope and intercept only when a formula is provably affine."""
        variables = {"sample": (1.0, 0.0)}
        variables.update({alias: (0.0, controls[name]) for alias, name in self.used_references.items()})
        return _affine_coefficients(self._tree, variables)

    def evaluate(self, sample: ArrayLike, controls: Mapping[str, float] | None = None) -> NumericResult:
        try:
            values = np.array(sample, dtype=float, copy=True)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError("Formule-invoer moet numeriek zijn.") from error
        if values.ndim > 1 or np.isinf(values).any():
            raise ValueError("Formule-invoer moet een enkel getal of een eendimensionale reeks zijn, zonder oneindige waarden.")
        controls = {} if controls is None else controls
        if not isinstance(controls, Mapping):
            raise ValueError("Formulecontroles moeten namen aan meetgemiddelden koppelen.")
        variables = {"sample": values}
        for alias in self._names - {"sample"}:
            name = self.references[alias]
            if name not in controls:
                raise ValueError(f"Geen meetgemiddelde beschikbaar voor controle {name!r}.")
            variables[alias] = finite_number(controls[name], f"control {name!r} mean")
        result = np.asarray(_evaluate(self._tree, variables), dtype=float)
        # A constant formula must not invent results for missing measurements.
        result = np.where(np.isnan(values), np.nan, np.broadcast_to(result, values.shape))
        return float(result) if result.ndim == 0 else result


def _affine_coefficients(node: ast.expr, variables: Mapping[str, tuple[float, float]]) -> tuple[float, float] | None:
    if isinstance(node, ast.Constant):
        return 0.0, float(node.value)
    if isinstance(node, ast.Name):
        return variables[node.id]
    if isinstance(node, ast.UnaryOp):
        value = _affine_coefficients(node.operand, variables)
        return tuple(-item for item in value) if value is not None and isinstance(node.op, ast.USub) else value
    left = _affine_coefficients(node.left, variables)
    right = _affine_coefficients(node.right, variables)
    if left is None or right is None:
        return None
    a, b = left
    c, d = right
    if isinstance(node.op, ast.Add):
        result = a + c, b + d
    elif isinstance(node.op, ast.Sub):
        result = a - c, b - d
    elif isinstance(node.op, ast.Mult) and (a == 0 or c == 0):
        result = a * d + b * c, b * d
    elif isinstance(node.op, ast.Div) and c == 0 and d != 0:
        result = a / d, b / d
    else:
        return None
    return result if all(isfinite(item) for item in result) else None


def _evaluate(node: ast.expr, variables: Mapping[str, ArrayLike]) -> ArrayLike:
    if isinstance(node, ast.Constant):
        return float(node.value)
    if isinstance(node, ast.Name):
        return variables[node.id]
    if isinstance(node, ast.UnaryOp):
        value = _evaluate(node.operand, variables)
        return np.negative(value) if isinstance(node.op, ast.USub) else value
    left, right = _evaluate(node.left, variables), _evaluate(node.right, variables)
    if isinstance(node.op, ast.Div) and np.any(np.asarray(right) == 0):
        raise ValueError("De formule mag niet door nul delen.")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        result = _OPERATORS[type(node.op)](left, right)
    if np.isinf(result).any():
        raise ValueError("De formule leverde oneindige waarden op.")
    return result
