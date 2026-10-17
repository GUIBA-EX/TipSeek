#!/usr/bin/env python3
"""Validate one sample sequence per SCO against the user-supplied reference.
Usage: report_sco_exon_recovery.py PANEL RUN BASELINE REFERENCE OUTPUT.json
All full-exon counts are cross-species coordinate proxies, not HX50 truth.
"""
import csv, json, statistics, sys
from pathlib import Path
from audit_exonerate_reference import fasta
from benchmark_exon_shadow import digest
from report_starfish_sco_recovery import summarize


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    panel_dir, run, baseline, references, destination = map(Path, sys.argv[1:])
    panel = json.loads((panel_dir / "panel.json").read_text())
    report = {
        "goal": "Recover exons; one reference-coordinate sample sequence per SCO; missing/conflicting bases are N.",
        "reference_root": str(references),
        "reference_files_in_root": len(list(references.glob("*.fasta"))),
        "scope": "12 preselected HX50 real-read SCO families, 185 ppect reference exons. Unchanged recruitment reused. Not a 7006-family validation. Full observation is not proof of same-copy linkage or exact biological exon boundaries in HX50.",
        "runs": json.loads((run / "runs.json").read_text()),
        "baseline_runs": json.loads((baseline / "runs.json").read_text()),
        "results": {},
    }
    for backend in ["miniprot", "exonerate-rs"]:
        root = run / backend
        audits = rows(root / "manifest/sco_recovery.tsv")
        details = []
        exon_rows = []
        for item in panel["families"]:
            family = item["family"]
            reference = references / (family + ".fasta")
            assert digest(reference) == digest(panel_dir / "references" / reference.name)
            ref = fasta(reference)["ppect"]
            recovered = fasta(root / "sco" / reference.name)
            assert len(recovered) == 1, family
            header, sequence = next(iter(recovered.items()))
            assert header.startswith("1_Hx50|" + family), header
            assert len(sequence) == len(ref), (family, len(sequence), len(ref))
            assert set(sequence) <= set("ACGTN"), family
            audit = next(r for r in audits if r["family_id"] == family)
            assert int(audit["length_nt"]) == len(sequence)
            assert sequence.count("N") == int(audit["missing_nt"]) + int(audit["conflict_nt"])
            exons = sorted(item["exons"], key=lambda x: x[1], reverse=item["exons"][0][3] == "-")
            offset = 0
            complete = 0
            for index, exon in enumerate(exons, 1):
                end = offset + exon[2] - exon[1]
                segment = sequence[offset:end]
                assert len(segment) == end - offset
                observed = len(segment) - segment.count("N")
                complete += observed == len(segment)
                exon_rows.append({"backend": backend, "family": family, "reference_exon": index, "reference_cds_start_0": offset, "reference_cds_end_0": end, "observed_nt": observed, "masked_nt": segment.count("N"), "fully_observed_reference_interval": observed == len(segment)})
                offset = end
            assert offset == len(ref)
            details.append({**audit, "fully_observed_reference_exons": complete, "reference_exons": len(exons), "observed_fraction": int(audit["observed_nt"]) / len(sequence), "reference_sha256": digest(reference), "output_sha256": digest(root / "sco" / reference.name)})
        evidence = rows(root / "manifest/exon_candidates.tsv")
        graph = rows(root / "manifest/chain_graph_audit.tsv")
        graph = [r for r in graph if r["protein_reference"].startswith("ppect__")]
        result = {"families": details, "reference_exons": sum(r["reference_exons"] for r in details), "fully_observed_reference_exons": sum(r["fully_observed_reference_exons"] for r in details), "mean_observed_nucleotide_fraction": statistics.mean(r["observed_fraction"] for r in details), "masked_conflict_nt": sum(int(r["conflict_nt"]) for r in details), "masked_missing_nt": sum(int(r["missing_nt"]) for r in details), "unmappable_models": sum(int(r["unmappable_models"]) for r in details), "exon_coordinate_candidates": len(evidence), "candidates_with_two_unique_mapped_reads": sum(int(r["unique_mapped_reads"]) >= 2 for r in evidence), "candidates_fully_covered_by_mapped_reads": sum(int(r["mapped_read_covered_nt"]) == int(r["length_nt"]) for r in evidence), "ppect_graph_audit": graph, "previous_alignment_proxy": summarize(baseline / backend, panel, backend), "current_alignment_proxy": summarize(root, panel, backend), "exon_intervals": exon_rows, "biological_wrong_merge_rate": None}
        report["results"][backend] = result
        print(backend, "full_reference_exons", result["fully_observed_reference_exons"], "/", result["reference_exons"], "observed_fraction", result["mean_observed_nucleotide_fraction"], "conflict_nt", result["masked_conflict_nt"], "unmappable", result["unmappable_models"])
    source_paths = ["rust/gm2_tools/src/gene_annotation/exon_recovery.rs", "rust/gm2_tools/src/gene_annotation/evidence.rs", "rust/gm2_tools/src/gene_annotation/exonerate.rs", "rust/gm2_tools/src/gene_annotation/model.rs", "rust/gm2_tools/src/gene_annotation/reads.rs", "rust/gm2_tools/src/gene_annotation/chain.rs", "rust/gm2_tools/src/bin/gene_workflow.rs", "rust/main_assembler_original/src/main.rs", "rust/tipseek_cli/src/main.rs"]
    report["source_sha256"] = {name: digest(Path(name)) for name in source_paths}
    report["reference_logic_sha256"] = {name: digest(Path(name)) for name in ["Captus/captus/bioformats.py", "HybPiper/hybpiper/exonerate_hits.py"]}
    destination.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
