#!/usr/bin/env python3
"""File-level provenance for an RNA Nextflow run."""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def checksum(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": h.hexdigest()
    }


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True, default=str) + "\n")


def input_files(params):
    for sample in params["samples"]:
        yield "sample:" + sample["id"], Path(sample["fastq"])

    for key in (
        "reference_fasta", "reference_gtf",
        "mane_transcripts", "clinvar_variants", "report_script"
    ):
        if params.get(key):
            yield key, Path(params[key])

# Include completed and cached Nextflow tasks when collecting output files.
def completed_workdirs(trace):
    if not trace.exists():
        return []

    with trace.open() as f:
        header = f.readline().lstrip("#").strip()
        if not header:
            return []

        rows = csv.DictReader([header, *f], delimiter="\t")
        return [
            Path(row["workdir"])
            for row in rows
            if row.get("status", "").upper() in ("COMPLETED", "CACHED")
            and row.get("workdir")
        ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=["start", "finish"])
    parser.add_argument("params_file", type=Path)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("exit_code", nargs="?", type=int, default=0)
    args = parser.parse_args()

    root = args.run_root.resolve()
    prov = root / "reproducibility"
    prov.mkdir(parents=True, exist_ok=True)

    params = json.loads(args.params_file.read_text())

    if args.phase == "start":
        inputs = []
        for role, path in input_files(params):
            if not path.is_file():
                raise FileNotFoundError(f"Missing {role}: {path}")
            inputs.append({"role": role, **checksum(path)})

        save_json(prov / "input_manifest.json", inputs)

        with (prov / "input_checksums.sha256").open("w") as f:
            for item in inputs:
                f.write(f"{item['sha256']}  {item['path']}\n")

        (prov / "started_utc.txt").write_text(
            datetime.now(timezone.utc).isoformat() + "\n"
        )
        return

    # Capture published outputs and successful Nextflow task outputs.
    files = set()

    for path in root.rglob("*"):
        if prov not in path.parents and path.is_file():
            files.add(path)

    trace = prov / "nextflow_trace.txt"

    for workdir in completed_workdirs(trace):
        if workdir.is_dir():
            for path in workdir.rglob("*"):
                if (
                    path.is_file()
                    and not path.is_symlink()
                    and not any(
                        part.startswith(".")
                        for part in path.relative_to(workdir).parts
                    )
                ):
                    files.add(path)

    outputs = []
    hash_cache = {}

    for path in sorted(files, key=str):
        real_path = path.resolve()
        if real_path not in hash_cache:
            hash_cache[real_path] = checksum(path)

        info = hash_cache[real_path]
        outputs.append({
            "path": str(path),
            "bytes": info["bytes"],
            "sha256": info["sha256"]
        })

    save_json(prov / "output_manifest.json", outputs)

    with (prov / "output_checksums.sha256").open("w") as f:
        for item in outputs:
            f.write(f"{item['sha256']}  {item['path']}\n")

    metadata_path = prov / "run_metadata.json"
    if not metadata_path.exists():
        for workdir in completed_workdirs(trace):
            candidate = workdir / "run_metadata.json"
            if candidate.is_file():
                metadata_path = candidate
                break

    inputs_path = prov / "input_manifest.json"

    manifest = {
        "dataset": "GSE287954",
        "samples": [s["id"] for s in params["samples"]],
        "preprocessing": params.get("rna_preprocessing"),
        "kit": params.get("rna_pychopper_kit"),
        "slurm_run": root.name,
        "status": "COMPLETED" if args.exit_code == 0 else "FAILED",
        "exit_code": args.exit_code,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_parameters": params,
        "run_metadata": (
            json.loads(metadata_path.read_text())
            if metadata_path.exists() else None
        ),
        "inputs": (
            json.loads(inputs_path.read_text())
            if inputs_path.exists() else []
        ),
        "outputs": outputs,
        "execution_records": [
            str(p) for p in sorted(prov.glob("nextflow_*"))
        ]
    }

    save_json(prov / "reproducibility_manifest.json", manifest)


if __name__ == "__main__":
    main()
