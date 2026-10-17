#!/usr/bin/env python3
"""Check frozen-cache invariants and summarize HX50 read-boundary controls."""
import argparse
import json
import re
import subprocess
import time
from collections import Counter
from pathlib import Path

from report_sco_projection import rows
from report_sco_projection import summarize, verify
from report_full_sco_exons import score


def compare_common_backbone(root):
    work = root / "n-boundary-fix"
    annotation = work / "exonerate-rs-merged-cache"
    backbone = root / "miniprot-backbones.tsv"
    recovered = work / "exonerate-rs-common-miniprot-backbone"
    manifest = recovered / "manifest/sco_recovery.tsv"
    if not manifest.exists():
        if recovered.exists():
            recovered.rename(work / (recovered.name + "-interrupted-" + str(time.time_ns())))
        command = [str(work / "components/gene_workflow"), "recover-sco", "--annotation", str(annotation),
                   "--input", str(work / "gene-full"), "--out", str(recovered), "--threads", "8",
                   "--sco-backbone-manifest", str(backbone)]
        with (work / "common-backbone-recovery.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    controls = {}
    mapped = {}
    for backend, recovery, source in [
        ("miniprot", work / "miniprot-corrected-recovery", root / "workflow/exon"),
        ("exonerate-rs", recovered, annotation),
    ]:
        audits = rows(recovery / "manifest/sco_recovery.tsv")
        verify(recovery, audits)
        mapped[backend] = {r["family_id"]: r for r in audits}
        exons = score(root / "exon_coordinates.json", recovery, source, backbone)
        assert exons["totals"]["scorable_families"] == 7006
        controls[backend] = {"sco": summarize(audits, 7006), "exon_totals": exons["totals"],
                             "complete_exon_interval_fraction": exons["complete_exon_interval_fraction"],
                             "short_exons": exons["short_exons_below_150nt"]}
    assert controls["miniprot"]["exon_totals"]["reference_exons"] == controls["exonerate-rs"]["exon_totals"]["reference_exons"]
    for family in set(mapped["miniprot"]) & set(mapped["exonerate-rs"]):
        a, b = mapped["miniprot"][family], mapped["exonerate-rs"][family]
        assert a["protein_reference"] == b["protein_reference"] and a["length_nt"] == b["length_nt"], family
    return {"backbone_manifest": str(backbone), "alignment_reruns": 0, "backends": controls,
            "scope": "Both corrected backends replayed on the same miniprot backbone; missing or conflicting sample evidence remains N."}


def load(path):
    return json.loads(path.read_text())


def timing(path):
    text = path.read_text()
    clock = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", text)
    parts = [float(p) for p in clock[1].split(":")] if clock else []
    seconds = 0
    for p in parts:
        seconds = seconds * 60 + p
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    return {"wall_seconds": seconds if parts else None, "peak_rss_kib": int(rss[1]) if rss else None}


def graph(annotation):
    totals = Counter()
    p = annotation / "manifest/chain_graph_audit.tsv"
    if p.exists():
        records = rows(p)
        totals["graph_records"] = len(records)
        for r in records:
            for field in ["eligible_sources", "ordered_edges", "compatible_edges", "coding_overlap_conflicts",
                          "coding_coordinate_map_missing", "other_rejected_edges", "read_vetoed_edges"]:
                totals[field] += int(r[field])
    p = annotation / "manifest/fragment_read_support.tsv"
    if p.exists():
        for r in rows(p):
            totals["read_link_status_" + r["status"]] += 1
            for field in ["concordant_pairs", "discordant_pairs", "orientation_conflicts"]:
                totals[field] += int(r[field])
    return dict(totals)


def report(root):
    work = root / "n-boundary-fix"
    candidates = load(work / "candidate-comparison.json")
    changed = {r["family"] for r in candidates["families"] if r["changed"]}
    result = {"reference_families": 7006, "candidate_comparison": {
        k: candidates[k] for k in ["affected_families", "changed_families", "unchanged_families"]},
        "scope": "Exact source CDS/GFF exon coordinates, not HX50 biological truth; conflict counts are not a measured false-merge rate.",
        "binaries": load(work / "binary-hashes.json"), "backends": {},
        "recruitment_output_equality": load(root / "recruit-optimized-equality.json"),
        "assembly_653_same_input": {"before": timing(work / "assembly-before.log"),
                                   "corrected": timing(work / "assembly.log")},
        "validation_stages": load(work / "validation-status.json")["stages"],
        "original_full_workflow_stages": load(root / "status.json")["stages"]}
    performance = work / "assembly-performance.json"
    if performance.exists():
        result["assembly_653_same_input"] = load(performance)
    repeat = work / "assembly-repeat-performance.json"
    if repeat.exists():
        result["assembly_653_repeated_controls"] = load(repeat)
    for backend in ["miniprot", "exonerate-rs"]:
        comparison_path = work / (backend + "-comparison.json")
        if not comparison_path.exists():
            result["backends"][backend] = {"status": "pending"}
            continue
        comparison = load(comparison_path)
        expanded = work / (backend + "-expanded-recovery")
        corrected = work / (backend + "-corrected-recovery")
        audit = [{r["family_id"]: r for r in rows(p / "manifest/sco_recovery.tsv")}
                 for p in [expanded, corrected]]
        unchanged = (set(audit[0]) | set(audit[1])) - changed
        for family in unchanged:
            assert family in audit[0] and family in audit[1], family
            assert (expanded / "sco" / (family + ".fasta")).read_bytes() == (corrected / "sco" / (family + ".fasta")).read_bytes(), family
            assert (expanded / "sco_masks" / (family + ".tsv")).read_bytes() == (corrected / "sco_masks" / (family + ".tsv")).read_bytes(), family
        for family in set(audit[0]) & set(audit[1]):
            assert audit[0][family]["protein_reference"] == audit[1][family]["protein_reference"]
            assert audit[0][family]["length_nt"] == audit[1][family]["length_nt"]
        exon_tables = [load(work / (backend + "-" + name + "-exons.json"))
                       for name in ["old", "expanded", "corrected"]]
        assert all(r["totals"]["scorable_families"] == 7006 for r in exon_tables)
        assert len({r["totals"]["reference_exons"] for r in exon_tables}) == 1
        family_tables = [{r["family"]: r for r in t["families"]} for t in exon_tables]
        improved, regressed = [], []
        for family in sorted(changed):
            a, b = [t.get(family, {}).get("complete_reference_exon_intervals", 0) for t in family_tables[1:]]
            if b != a:
                (improved if b > a else regressed).append({"family": family, "before": a, "after": b})
        original = root / ("workflow/exon" if backend == "miniprot" else "exonerate-rs")
        baseline = Path("/home/yasuda2/tipseek_HX50_all7006_familyqueue_20261011") / (
            "miniprot" if backend == "miniprot" else "exonerate-rs")
        merged = work / (backend + "-merged-cache")
        result["backends"][backend] = {"status": "completed", **comparison,
            "unchanged_family_outputs_and_masks_verified": len(unchanged),
            "cache_provenance": {k: v for k, v in load(merged / "cache-provenance.json").items() if k != "families"},
            "boundary_fix_complete_exons_improved": improved,
            "boundary_fix_complete_exons_regressed": regressed,
            "old_graph_and_read_links": graph(baseline),
            "expanded_graph_and_read_links": graph(original), "corrected_graph_and_read_links": graph(merged)}
    result["status"] = "completed" if all(r["status"] == "completed" for r in result["backends"].values()) else "partial"
    timings = {s["name"]: s["seconds"] for s in result["validation_stages"]}
    if all(b + "-annotate-changed" in timings for b in ["miniprot", "exonerate-rs"]):
        result["matched_29_family_annotation_performance"] = {
            "scope": "Same 29 corrected families, reads, references, binary, and 8-worker budget; not full-workflow timing.",
            "miniprot_seconds": timings["miniprot-annotate-changed"],
            "exonerate_rs_seconds": timings["exonerate-rs-annotate-changed"],
            "exonerate_to_miniprot_time_ratio": timings["exonerate-rs-annotate-changed"] / timings["miniprot-annotate-changed"],
        }
    if result["status"] == "completed":
        result["common_backbone_backend_comparison"] = compare_common_backbone(root)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = report(args.root)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "backends": {
        b: {"status": r["status"], "comparisons": r.get("comparisons")} for b, r in result["backends"].items()}}, indent=2))
