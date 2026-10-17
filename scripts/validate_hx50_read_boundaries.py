#!/usr/bin/env python3
"""Repair ambiguous-read families and replay the remaining frozen HX50 alignments.

The original workflow is never modified. Each backend uses one fixed backbone
manifest for the old, expanded-recruitment, and read-boundary-corrected controls.
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import time
import traceback
from pathlib import Path

from report_full_sco_exons import score
from report_sco_projection import fasta, rows, summarize, verify


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def link(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        assert destination.resolve() == source.resolve(), destination
    elif not destination.exists():
        destination.symlink_to(source.resolve(), target_is_directory=source.is_dir())
    else:
        raise RuntimeError(f"unexpected existing destination: {destination}")


def atomic_json(path, data):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(path)


class Validation:
    def __init__(self, root, threads):
        self.root = root.resolve()
        self.work = self.root / "n-boundary-fix"
        self.gene = self.work / "components/gene_workflow"
        self.threads = str(threads)
        self.state_path = self.work / "validation-status.json"
        self.state = json.loads(self.state_path.read_text()) if self.state_path.exists() else {
            "status": "running", "started_unix": time.time(), "stages": [],
        }
        self.state.update(status="running", error=None)
        self.affected = {f for f, _ in rows_json(self.root / "ambiguous-read-audit.json")["families"]}
        self.environment = dict(os.environ)
        self.environment["PATH"] = str(self.root / "components") + os.pathsep + os.environ["PATH"]

    def save(self, phase):
        self.state["phase"] = phase
        atomic_json(self.state_path, self.state)

    def run(self, name, command):
        previous = next((s for s in self.state["stages"] if s["name"] == name), None)
        if previous:
            assert previous["returncode"] == 0 and previous["command"] == command, name
            return
        self.save(name)
        started = time.monotonic()
        with (self.work / (name + ".log")).open("w") as stream:
            process = subprocess.Popen(["/usr/bin/time", "-v", *command], stdout=stream,
                                       stderr=subprocess.STDOUT, env=self.environment)
            self.state["child_pid"] = process.pid
            self.save(name)
            code = process.wait()
        self.state["stages"].append({"name": name, "command": command,
                                     "seconds": time.monotonic() - started, "returncode": code})
        self.save(name)
        if code:
            raise RuntimeError(f"{name} failed ({code}); inspect its log")

    def wait_original(self, backend):
        name = "miniprot-fixed-recovery" if backend == "miniprot" else "exonerate-rs-fixed-recovery"
        while True:
            state = rows_json(self.root / "status.json")
            if any(s["name"] == name and s["returncode"] == 0 for s in state["stages"]):
                return
            if state["status"] == "failed":
                raise RuntimeError(f"original workflow failed: {state.get('error')}")
            pid = state.get("runner_pid", state.get("child_pid"))
            if state["status"] == "running" and pid:
                process = Path("/proc") / str(pid) / "status"
                if not process.exists() or any(
                    line.startswith("State:") and "Z" in line
                    for line in process.read_text().splitlines()
                ):
                    raise RuntimeError(f"workflow status is stale: process {pid} has exited")
            self.save("waiting-for-" + name)
            time.sleep(15)

    def prepare(self):
        while True:
            text = (self.work / "assembly.log").read_text()
            if "Exit status: 0" in text:
                break
            if "Exit status:" in text:
                raise RuntimeError("corrected assembly failed")
            self.save("waiting-for-corrected-assembly")
            time.sleep(5)
        reference = self.work / "reference"
        subset = self.work / "gene-653"
        self.run("classify-653", [str(self.gene), "classify", "--reference", str(reference),
            "--contigs", str(self.work / "assembly/contigs_all"), "--sample", "1_Hx50",
            "--out", str(subset), "--retain-exon-fragments"])
        self.run("cohort-653", [str(self.gene), "cohort", "--reference", str(reference),
            "--out", str(subset), "--sample", "1_Hx50"])
        old = self.root / "workflow/gene/samples/1_Hx50/candidates"
        new = subset / "samples/1_Hx50/candidates"
        flags = []
        for directory in [old, new]:
            unanchored = {}
            with (directory.parent / "candidate_assembly.tsv").open() as stream:
                for r in csv.DictReader(stream, delimiter="\t"):
                    if r["family_id"] in self.affected and r["assembler_seed_position"] == "-2":
                        unanchored.setdefault(r["family_id"], set()).add(r["candidate"])
            flags.append(unanchored)
        self.changed = set()
        comparison = []
        for family in sorted(self.affected):
            a, b = old / (family + ".fasta"), new / (family + ".fasta")
            aa, bb = fasta(a) if a.exists() else {}, fasta(b) if b.exists() else {}
            screening_changed = flags[0].get(family, set()) != flags[1].get(family, set())
            changed = aa != bb or screening_changed
            if changed:
                self.changed.add(family)
            comparison.append({"family": family, "old_candidates": len(aa), "new_candidates": len(bb),
                               "changed": changed, "screening_changed": screening_changed,
                               "old_sha256": digest(a) if a.exists() else None,
                               "new_sha256": digest(b) if b.exists() else None})
        atomic_json(self.work / "candidate-comparison.json", {
            "affected_families": len(self.affected), "changed_families": len(self.changed),
            "unchanged_families": len(self.affected - self.changed), "families": comparison})
        changed_gene = self.work / "gene-changed"
        for p in new.glob("*.fasta"):
            if p.stem in self.changed:
                link(p, changed_gene / "samples/1_Hx50/candidates" / p.name)
        # An empty candidate directory is still a valid input when all corrected
        # families disappear; no lost family is silently restored from its cache.
        (changed_gene / "samples/1_Hx50/candidates").mkdir(parents=True, exist_ok=True)
        # The annotation screen distinguishes unanchored components using this
        # table. Omitting it would change model acceptance in the subset run.
        link(subset / "samples/1_Hx50/candidate_assembly.tsv",
             changed_gene / "samples/1_Hx50/candidate_assembly.tsv")
        merged_gene = self.work / "gene-full"
        (merged_gene / "samples/1_Hx50/candidates").mkdir(parents=True, exist_ok=True)
        for source in [old, new]:
            for p in source.glob("*.fasta"):
                if (source == old and p.stem not in self.changed) or (source == new and p.stem in self.changed):
                    link(p, merged_gene / "samples/1_Hx50/candidates" / p.name)
        self.state.update(affected_families=len(self.affected), changed_families=len(self.changed))
        self.save("candidates-compared")

    def merge(self, backend, original, corrected):
        destination = self.work / (backend + "-merged-cache")
        family_sources = {}
        for source in [original, corrected]:
            for p in (source / "sco").glob("*.fasta"):
                if (source == original and p.stem not in self.changed) or (source == corrected and p.stem in self.changed):
                    family_sources[p.stem] = source
                    link(p, destination / "sco" / p.name)
        for family, source in family_sources.items():
            for raw in ["raw_miniprot", "raw_exonerate"]:
                p = source / raw / "1_Hx50" / family
                if p.is_dir():
                    link(p, destination / raw / "1_Hx50" / family)
            p = source / "manifest/derived_proteins" / (family + ".faa")
            if p.exists():
                link(p, destination / "manifest/derived_proteins" / p.name)
        manifest = destination / "manifest"
        manifest.mkdir(parents=True, exist_ok=True)
        # Merge all family-addressed evidence, not just the replay prerequisites,
        # so the read-link and conflict audit remains traceable to its source.
        for name in sorted({p.name for s in [original, corrected] for p in (s / "manifest").glob("*.tsv")}):
            tables = []
            for source in [original, corrected]:
                p = source / "manifest" / name
                if p.exists():
                    with p.open() as stream:
                        reader = csv.DictReader(stream, delimiter="\t")
                        header = reader.fieldnames
                        if "family_id" in header:
                            tables.append((source, header, list(reader)))
            if not tables:
                continue
            header = list(dict.fromkeys(h for _, hs, _ in tables for h in hs))
            with (manifest / name).open("w") as stream:
                writer = csv.DictWriter(stream, header, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for source, _, records in tables:
                    for r in records:
                        if family_sources.get(r["family_id"]) == source:
                            writer.writerow(r)
        atomic_json(destination / "cache-provenance.json", {
            "families": {f: str(s) for f, s in sorted(family_sources.items())},
            "reused_family_count": sum(s == original for s in family_sources.values()),
            "realigned_family_count": sum(s == corrected for s in family_sources.values()),
        })
        return destination

    def backend(self, backend):
        corrected = self.work / (backend + "-changed")
        if self.changed:
            self.run(backend + "-annotate-changed", [str(self.gene), "annotate", "--backend", backend,
                "--input", str(self.work / "gene-changed"), "--nucleotide-reference", str(self.work / "reference"),
                "--out", str(corrected), "--reads-root", str(self.root / "workflow"), "--threads", self.threads])
        else:
            (corrected / "sco").mkdir(parents=True, exist_ok=True)
            (corrected / "manifest").mkdir(parents=True, exist_ok=True)
        self.wait_original(backend)
        original = self.root / ("workflow/exon" if backend == "miniprot" else "exonerate-rs")
        merged = self.merge(backend, original, corrected)
        backbone = self.root / (backend + "-backbones.tsv")
        base_root = Path("/home/yasuda2/tipseek_HX50_all7006_familyqueue_20261011")
        baseline = base_root / ("miniprot" if backend == "miniprot" else "exonerate-rs")
        comparisons = []
        for name, annotation, gene in [
            ("old", baseline, Path("/home/yasuda2/tipseek_HX50_all7006_20261011/workflow/gene")),
            ("expanded", original, self.root / "workflow/gene"),
            ("corrected", merged, self.work / "gene-full"),
        ]:
            recovery = self.work / (backend + "-" + name + "-recovery")
            self.run(backend + "-replay-" + name, [str(self.gene), "recover-sco", "--annotation", str(annotation),
                "--input", str(gene), "--out", str(recovery), "--threads", self.threads,
                "--sco-backbone-manifest", str(backbone)])
            audits = rows(recovery / "manifest/sco_recovery.tsv")
            verify(recovery, audits)
            report = score(self.root / "exon_coordinates.json", recovery, annotation, backbone)
            atomic_json(self.work / (backend + "-" + name + "-exons.json"), report)
            comparisons.append({"name": name, "sco": summarize(audits, 7006),
                                "exon_totals": report["totals"], "short_exons": report["short_exons_below_150nt"],
                                "complete_exon_interval_fraction": report["complete_exon_interval_fraction"],
                                "coordinate_conflict_codons": sum(int(r["coordinate_conflict_codons"]) for r in audits),
                                "dominated_models": sum(int(r["dominated_models"]) for r in audits)})
        atomic_json(self.work / (backend + "-comparison.json"), {
            "backend": backend, "backbone_manifest": str(backbone), "backbone_sha256": digest(backbone),
            "scope": "Reference exon coordinate completeness, not HX50 biological truth or false-merge rate.",
            "comparisons": comparisons,
        })

    def execute(self):
        try:
            self.prepare()
            for backend in ["miniprot", "exonerate-rs"]:
                self.backend(backend)
            self.state.update(status="completed", finished_unix=time.time())
            self.save("completed")
        except Exception:
            self.state.update(status="failed", error=traceback.format_exc(), finished_unix=time.time())
            self.save("failed")
            raise


def rows_json(path):
    return json.loads(path.read_text())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    Validation(args.root, args.threads).execute()
