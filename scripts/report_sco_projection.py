#!/usr/bin/env python3
"""Audit saved-alignment SCO recovery with selection-only and fixed-backbone controls.
Usage: report_sco_projection.py BASELINE RECOVERED FIXED_OLD REFERENCES OUT.json
"""
import csv
import hashlib
import json
import statistics
import sys
from pathlib import Path


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def fasta(path):
    records = {}
    header = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:]
            assert header not in records, (path, header)
            records[header] = ""
        elif header is not None:
            records[header] += line.strip().upper()
    return records


def summarize(audits, total, observed="observed_nt"):
    fractions = [int(r[observed]) / int(r["length_nt"]) for r in audits]
    return {
        "reference_families": total,
        "outputs": len(audits),
        "no_output": total - len(audits),
        "any_observed": sum(x > 0 for x in fractions),
        "all_N": sum(x == 0 for x in fractions),
        "at_least_50_percent": sum(x >= .5 for x in fractions),
        "at_least_80_percent": sum(x >= .8 for x in fractions),
        "at_least_95_percent": sum(x >= .95 for x in fractions),
        "entire_sequence_no_N": sum(x == 1 for x in fractions),
        "observed_nt": sum(int(r[observed]) for r in audits),
        "mean_observed_fraction_including_no_output": sum(fractions) / total,
        "total_coordinate_nt": sum(int(r["length_nt"]) for r in audits),
    }


def verify(root, audits):
    assert len({(r["sample"], r["family_id"]) for r in audits}) == len(audits)
    for r in audits:
        records = fasta(root / "sco" / (r["family_id"] + ".fasta"))
        matches = [(h,s) for h,s in records.items() if h.startswith(r["sample"] + "|")]
        assert len(matches) == 1, r
        h,s = matches[0]
        assert "reference=" + r["protein_reference"] + " " in h
        assert len(s) == int(r["length_nt"])
        assert set(s) <= set("ACGTN")
        assert len(s) - s.count("N") == int(r["observed_nt"])
        assert s.count("N") == int(r["missing_nt"]) + int(r["conflict_nt"])


def main():
    baseline, recovered, fixed, references, destination = map(Path, sys.argv[1:])
    old = rows(baseline / "manifest/sco_recovery.tsv")
    new = rows(recovered / "manifest/sco_recovery.tsv")
    frozen = rows(fixed / "manifest/sco_recovery.tsv")
    direct = rows(recovered / "manifest/sco_selection_control.tsv")
    for root,audits in [(baseline,old),(recovered,new),(fixed,frozen)]:
        verify(root,audits)
    keys = lambda rs: {(r["sample"],r["family_id"]):r for r in rs}
    o,n,f = map(keys,[old,new,frozen])
    assert o.keys() == n.keys() == f.keys()
    assert all(o[k]["protein_reference"] == f[k]["protein_reference"] and o[k]["length_nt"] == f[k]["length_nt"] for k in o)
    ref_files = sorted(references.glob("*.fasta"))
    total = len(ref_files)
    maps = rows(recovered / "manifest/reference_projection.tsv")
    statuses = {status:sum(r["status"] == status for r in maps) for status in sorted({r["status"] for r in maps})}
    report = {
        "scope": "HX50 real-read assemblies; reference-coordinate sample evidence composites, not validated same-copy CDS, exon boundaries or biological false-merge rates.",
        "reference_root": str(references),
        "reused_alignment_root": str(baseline),
        "recovered_root": str(recovered),
        "fixed_backbone_root": str(fixed),
        "baseline": summarize(old,total),
        "selection_only": summarize(direct,total,"direct_only_observed_nt"),
        "selection_and_projection": summarize(new,total),
        "projection_on_unchanged_old_backbone": summarize(frozen,total),
        "changed_backbones": sum(o[k]["protein_reference"] != n[k]["protein_reference"] for k in o),
        "median_coordinate_length_ratio_new_to_old": statistics.median(int(n[k]["length_nt"]) / int(o[k]["length_nt"]) for k in o),
        "fixed_backbone_families_with_more_observed_nt": sum(int(f[k]["observed_nt"]) > int(o[k]["observed_nt"]) for k in o),
        "fixed_backbone_families_with_fewer_observed_nt": sum(int(f[k]["observed_nt"]) < int(o[k]["observed_nt"]) for k in o),
        "coordinate_conflict_codons": sum(int(r["coordinate_conflict_codons"]) for r in new),
        "masked_conflict_nt": sum(int(r["conflict_nt"]) for r in new),
        "projection_statuses": statuses,
        "biological_wrong_merge_rate": None,
        "validation": "Verified one record per sample/SCO, reference IDs, lengths, ACGTN alphabet, observed/missing/conflict accounting, and unchanged fixed-backbone lengths/IDs in every output.",
    }
    sources = ["rust/gm2_tools/src/bin/gene_workflow.rs","rust/gm2_tools/src/gene_annotation/exon_recovery.rs","rust/gm2_tools/src/gene_annotation/protein_projection.rs","rust/gm2_tools/src/gene_annotation/sco_rebuild.rs","rust/gm2_tools/src/gene_annotation/family_scheduler.rs","Captus/captus/bioformats.py","HybPiper/hybpiper/distribute_targets.py","HybPiper/hybpiper/exonerate_hits.py"]
    report["source_sha256"] = {p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources}
    destination.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:report[k] for k in ["baseline","selection_only","selection_and_projection","projection_on_unchanged_old_backbone"]},indent=2))


if __name__ == "__main__":
    main()
