"""One notation for formula names in reports, summaries and editable equations."""

import keyword
import re
from ast import literal_eval
from collections.abc import Mapping
from io import StringIO
from tokenize import NAME, STRING, TokenError, generate_tokens, untokenize
from unicodedata import normalize

from elisaflow.formulas.expressions import Formula


def _symbol(name: str) -> str:
    """Quote names that could otherwise be mistaken for arithmetic or keywords."""
    return name if name.isidentifier() and not keyword.iskeyword(name) and name not in {"sample", "x"} else repr(name)


def _translate(expression: str, names: Mapping[str, str], *, parse: bool = False) -> str:
    try:
        tokens = []
        for token in generate_tokens(StringIO(expression).readline):
            if token.type == NAME or parse and token.type == STRING:
                name = repr(literal_eval(token.string)) if token.type == STRING else token.string
                if name in names:
                    token = token._replace(string=names[name])
                elif parse:
                    raise ValueError(f"Onbekende formulereferenties: {token.string}.")
            tokens.append(token)
        return untokenize(tokens)
    except (TokenError, IndentationError, SyntaxError) as error:
        raise ValueError("Ongeldige rekenkundige formule.") from error


class FormulaNames:
    """Translate actual sample/control names to and from stored formula aliases."""

    def __init__(self, sample_name: str, references: Mapping[str, str]):
        self.sample_name = sample_name.strip()
        self.references = dict(references)
        self._control_symbols: dict[str, str] = {}
        for name in self.references.values():
            self._control_symbol(name)

    @property
    def sample_symbol(self) -> str:
        return _symbol(self.sample_name)

    def _control_symbol(self, name: str) -> str:
        if name not in self._control_symbols:
            symbol = _symbol(name)
            used = {self.sample_symbol, *self._control_symbols.values()}
            suffix = 1
            while symbol in used:
                label = "controlegemiddelde" if suffix == 1 else f"controlegemiddelde {suffix}"
                symbol = _symbol(f"{name} ({label})")
                suffix += 1
            self._control_symbols[name] = symbol
        return self._control_symbols[name]

    def control_symbol(self, name: str) -> str:
        """Bind a selected control internally and return its visible name."""
        if name not in self.references.values():
            stem = re.sub(r"\W+", "_", normalize("NFKC", name).casefold()).strip("_") or "controle"
            if not stem.isidentifier() or keyword.iskeyword(stem) or stem in {"sample", "x"}:
                stem = f"controle_{stem}"
            alias = stem
            suffix = 2
            while alias in self.references:
                alias = f"{stem}_{suffix}"
                suffix += 1
            self.references[alias] = name
        return self._control_symbol(name)

    def _symbols(self) -> dict[str, str]:
        return {"sample": self.sample_symbol, **{alias: self._control_symbol(name) for alias, name in self.references.items()}}

    def format(self, expression: str) -> str:
        """Render stored arithmetic with the exact control names from its bindings."""
        return _translate(expression, self._symbols())

    def parse(self, expression: str) -> Formula:
        """Accept the same visible names that formatting and insertion produce."""
        aliases = {}
        for alias, symbol in self._symbols().items():
            aliases.setdefault(symbol, alias)
        return Formula(_translate(expression.strip(), aliases, parse=True), self.references)
