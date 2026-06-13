#!/usr/bin/env python3
"""Write the package versions used by the figure workflow."""

from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import json
import shutil
import subprocess
from pathlib import Path


PACKAGES = ["eon", "readcon", "rgpycrumbs", "chemparseplot", "matplotlib", "numpy"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    versions = {}
    for package in PACKAGES:
        package_info = {}
        try:
            importlib.import_module(package)
            package_info["import"] = "ok"
        except ImportError as exc:
            package_info["import"] = f"failed: {exc}"
        try:
            package_info["version"] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            package_info["version"] = "metadata unavailable"
        versions[package] = package_info

    eonclient = shutil.which("eonclient")
    versions["eonclient"] = {"path": eonclient or "not found"}
    if eonclient:
        proc = subprocess.run(
            [eonclient, "--version"],
            check=False,
            text=True,
            capture_output=True,
        )
        versions["eonclient"]["version_output"] = (
            proc.stdout.strip() or proc.stderr.strip() or f"exit {proc.returncode}"
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(versions, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
