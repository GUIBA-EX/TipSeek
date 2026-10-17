#!/usr/bin/env python3
"""Resume the interrupted frozen HX50 family queue without discarding its cache."""
import concurrent.futures
import csv
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

from validate_hx50_read_boundaries import atomic_json


def merge_tree(source, target, seen):
    target.mkdir(parents=True, exist_ok=True)
    for p in sorted(source.iterdir()):
        destination = target / p.name
        if p.is_dir():
            merge_tree(p, destination, seen)
        elif p.name == "execution.tsv":
            continue
        elif p.suffix == ".tsv":
            with p.open() as stream:
                header = stream.readline()
                if destination.exists():
                    with destination.open() as existing:
                        assert existing.readline() == header, destination
                with destination.open("a") as output:
                    if destination.stat().st_size == 0:
                        output.write(header)
                    unique = p.name in {"reference_translation.tsv", "backend.tsv"}
                    values = seen.setdefault(str(destination), set()) if unique else None
                    for line in stream:
                        if not unique or line not in values:
                            output.write(line)
                            if unique:
                                values.add(line)
        else:
            assert not destination.exists(), destination
            try:
                os.link(p, destination)
            except OSError:
                shutil.copyfile(p, destination)


def complete(task):
    out = task / "output"
    return (out / "manifest/long_multiple_candidates.tsv").exists() and not (out / ".work").exists()


def main(root):
    status_path = root / "status.json"
    state = json.loads(status_path.read_text())
    if state["status"] == "completed":
        print("Frozen full-library workflow is already complete.")
        return
    old_pid = state.get("runner_pid", state.get("child_pid"))
    if state["status"] == "running" and old_pid and Path("/proc", str(old_pid)).exists():
        raise RuntimeError(f"workflow process {old_pid} is still present; refusing concurrent resume")
    archive = root / "status-before-resume.json"
    if not archive.exists():
        atomic_json(archive, state)
    tasks = sorted((root / "exonerate-rs/.family_workers").glob("task_*"))
    assert len(tasks) == 6935, len(tasks)
    finished = [t for t in tasks if complete(t)]
    pending = [t for t in tasks if not complete(t)]
    state.update(status="running", phase="exonerate-resume", runner_pid=os.getpid(),
                 child_pid=os.getpid(), resume_started_unix=time.time(),
                 preserved_family_count=len(finished), remaining_family_count=len(pending),
                 completed_family_count=len(finished), error=None)
    atomic_json(status_path, state)
    bins = root / "components"
    environment = dict(os.environ, PATH=str(bins) + os.pathsep + os.environ["PATH"],
                       GM2_COMPONENT_DIR=str(bins))
    log_path = root / "exonerate-all7006.log"
    with log_path.open("a") as stream:
        stream.write(f"Resume: preserved {len(finished)} completed families; {len(pending)} remaining.\n")
    started = time.monotonic()

    def annotate(task):
        out = task / "output"
        if out.exists():
            out.rename(task / ("interrupted-output-" + str(time.time_ns())))
        command = [str(bins / "gene_workflow"), "annotate", "--backend", "exonerate-rs",
                   "--input", str(task / "input"), "--nucleotide-reference", str(task / "references"),
                   "--out", str(out), "--reads-root", str(root / "workflow"), "--threads", "1"]
        with (task / "resume.log").open("w") as log:
            process = subprocess.Popen(["/usr/bin/time", "-v", *command], env=environment,
                                       stdout=log, stderr=subprocess.STDOUT)
            atomic_json(task / "resume-process.json", {"pid": process.pid, "command": command,
                                                       "started_unix": time.time()})
            code = process.wait()
        if code or not complete(task):
            raise RuntimeError(f"{task.name}: annotate failed ({code}); see resume.log")
        atomic_json(task / "resume-completed.json", {"returncode": code, "finished_unix": time.time()})
        return task

    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
            futures = [pool.submit(annotate, task) for task in pending]
            count = len(finished)
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception:
                    state.update(status="failed", error=traceback.format_exc())
                    atomic_json(status_path, state)
                    for queued in futures:
                        queued.cancel()
                    raise
                count += 1
                state["completed_family_count"] = count
                atomic_json(status_path, state)
                if count % 100 == 0 or count == len(tasks):
                    with log_path.open("a") as stream:
                        stream.write(f"SCO annotation completed {count}/{len(tasks)} (resumed)\n")
        state["stages"].append({"name": "exonerate-all7006", "seconds": time.monotonic()-started,
            "returncode": 0, "seconds_scope": "remaining families after interruption, not full annotation",
            "preserved_families": len(finished), "realigned_families": len(pending),
            "full_annotation_wall_seconds": None})
        state["phase"] = "exonerate-merge-resumed"
        atomic_json(status_path, state)
        staging = root / "exonerate-resume-merge"
        if staging.exists():
            staging.rename(root / ("exonerate-resume-merge-incomplete-" + str(time.time_ns())))
        seen = {}
        for task in tasks:
            assert complete(task), task
            merge_tree(task / "output", staging, seen)
        (staging / "manifest/execution.tsv").write_text(
            "candidate_workers\tthreads_per_candidate\tfamily_workers\tscheduling\n"
            "16\t1\t16\tresumed_frozen_family_queue\n")
        atomic_json(staging / "manifest/resume_summary.json", {
            "preserved_family_count": len(finished), "remaining_family_count": len(pending),
            "original_interruption_time_unknown": True, "full_annotation_wall_seconds": None,
            "task_checkpoints_retained": True})
        destination = root / "exonerate-rs"
        for p in sorted(staging.iterdir()):
            assert not (destination / p.name).exists(), p
            p.rename(destination / p.name)
        staging.rmdir()
        old = Path("/home/yasuda2/tipseek_HX50_sco_projection_20261011/exonerate-rs-optimized/manifest/sco_recovery.tsv")
        with old.open() as stream:
            backbones = {r["family_id"]: r["protein_reference"] for r in csv.DictReader(stream, delimiter="\t")}
        with (destination / "manifest/sco_recovery.tsv").open() as stream:
            for r in csv.DictReader(stream, delimiter="\t"):
                backbones.setdefault(r["family_id"], r["protein_reference"])
        manifest = root / "exonerate-rs-backbones.tsv"
        with manifest.open("w") as stream:
            stream.write("family_id\tprotein_reference\n")
            for family, reference in sorted(backbones.items()):
                stream.write(f"{family}\t{reference}\n")
        state["phase"] = "exonerate-rs-fixed-recovery"
        atomic_json(status_path, state)
        command = [str(bins / "gene_workflow"), "recover-sco", "--annotation", str(destination),
                   "--input", str(root / "workflow/gene"), "--out", str(root / "exonerate-rs-fixed"),
                   "--threads", "16", "--sco-backbone-manifest", str(manifest)]
        started = time.monotonic()
        with (root / "exonerate-rs-fixed-recovery.log").open("w") as log:
            subprocess.run(["/usr/bin/time", "-v", *command], env=environment,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        state["stages"].append({"name": "exonerate-rs-fixed-recovery", "command": command,
                                "seconds": time.monotonic()-started, "returncode": 0})
        state.update(status="completed", phase="completed", finished_unix=time.time())
        atomic_json(status_path, state)
    except Exception:
        state.update(status="failed", error=traceback.format_exc(), finished_unix=time.time())
        atomic_json(status_path, state)
        raise


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
