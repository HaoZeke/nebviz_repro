#!/usr/bin/env python3
"""Write the package versions used by the figure workflow."""

from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import json
from pathlib import Path


PACKAGES = ["eon", "readcon", "rgpycrumbs", "chemparseplot", "matplotlib", "numpy"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    versions = {}
    for package in PACKAGES:
        importlib.import_module(package)
        versions[package] = importlib.metadata.version(package)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(versions, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
