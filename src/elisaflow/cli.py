"""Open ElisaFlow or start a measurement directly."""

import argparse
import logging
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import BadZipFile

from elisaflow.configuration.application import ApplicationSettings, default_config_path
from elisaflow.diagnostics import configure_logging, show_error
from elisaflow.ui.desktop import prompt_analyst_name


def main(argv: list[str] | None = None) -> int:
    try:
        log_path = configure_logging()
    except OSError as error:
        show_error(f"Kan het logbestand niet openen: {error}")
        return 1

    parser = argparse.ArgumentParser(prog="elisaflow", description=__doc__)
    parser.add_argument("command", nargs="?", choices=("run", "settings"), help="Start a run or open assay settings; default: desktop.")
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--protocol", type=Path, help="Protocol for a live measurement; otherwise show the selection dialog.")
    inputs.add_argument("--experiment", type=Path, help="Reanalyze a saved experiment folder using its completed SKAX and settings snapshot.")
    parser.add_argument("--output", type=Path, help="PDF path; defaults to a new run folder.")
    args = parser.parse_args(argv)

    run_options = args.protocol is not None or args.experiment is not None or args.output is not None
    if args.command != "run" and run_options:
        parser.error("Run options require the 'run' command")

    try:
        config = ApplicationSettings.from_file(default_config_path())
        if args.command == "run":
            from elisaflow.runs.workflow import run_workflow

            analyst = prompt_analyst_name()
            if analyst is None:
                return 0
            run_workflow(
                config,
                analyst=analyst,
                protocol_template=args.protocol,
                experiment_dir=args.experiment,
                output=args.output,
            )
            return 0

        from elisaflow.ui.application import launch

        return launch(config, show_settings=args.command == "settings")
    except (OSError, ValueError, RuntimeError, KeyError, ET.ParseError, BadZipFile) as error:
        logging.getLogger(__name__).exception("%s", error)
        if log_path is not None:
            show_error(str(error), log_path)
        return 1
