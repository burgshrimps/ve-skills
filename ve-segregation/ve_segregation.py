#!/usr/bin/env python3
"""ve-segregation: parent-side segregation for filtered pedigree VCFs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable


VERSION = "0.1.0"
DEFAULT_CASE = "ISDBM322015"
DEFAULT_FATHER = "ISDBM322016"
DEFAULT_MOTHER = "ISDBM322018"
DEFAULT_SIBLING = "ISDBM322017"
DEFAULT_ASSEMBLY = "GRCh37/b37"


VARIANT_COLUMNS = [
    "variant_key",
    "chrom",
    "pos",
    "id",
    "ref",
    "alt",
    "gene",
    "transcript",
    "effect",
    "effect_impact",
    "segregation_side",
    "case_sample",
    "father_sample",
    "mother_sample",
    "sibling_sample",
    "case_gt",
    "father_gt",
    "mother_gt",
    "sibling_gt",
    "case_dp",
    "father_dp",
    "mother_dp",
    "sibling_dp",
    "case_gq",
    "father_gq",
    "mother_gq",
    "sibling_gq",
    "passes_dp_gq",
    "passes_carrier_logic",
    "gene_both_parental_sides_high_effect",
    "segregation_flag",
]

GENE_COLUMNS = [
    "gene",
    "both_parental_sides_high_effect",
    "segregation_flag",
    "paternal_high_effect_variant_count",
    "maternal_high_effect_variant_count",
    "total_high_effect_variant_count",
    "paternal_variant_ids",
    "maternal_variant_ids",
    "paternal_loci",
    "maternal_loci",
    "paternal_effects",
    "maternal_effects",
    "interpretation_boundary",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def parse_info(info: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = "true"
    return parsed


def parse_eff(info: str) -> list[dict[str, str]]:
    eff = parse_info(info).get("EFF", "")
    annotations: list[dict[str, str]] = []
    if not eff:
        return annotations
    for raw in eff.split(","):
        if "(" not in raw or not raw.endswith(")"):
            continue
        effect, rest = raw.split("(", 1)
        fields = rest[:-1].split("|")
        impact = fields[0] if len(fields) > 0 else ""
        gene = fields[5] if len(fields) > 5 else ""
        transcript = fields[8] if len(fields) > 8 else ""
        if not gene:
            continue
        annotations.append(
            {
                "gene": gene,
                "effect": effect,
                "impact": impact,
                "transcript": transcript,
                "raw": raw,
            }
        )
    return annotations


def parse_sample(format_keys: list[str], sample_value: str) -> dict[str, str | int | None]:
    values = sample_value.split(":")
    data = {key: values[idx] if idx < len(values) else "" for idx, key in enumerate(format_keys)}
    return {
        "GT": data.get("GT", ""),
        "DP": parse_int(data.get("DP", "")),
        "GQ": parse_int(data.get("GQ", "")),
    }


def parse_int(value: str | None) -> int | None:
    if value in (None, "", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def carries_alt(gt: str | None) -> bool:
    if not gt or gt in {".", "./.", ".|."}:
        return False
    alleles = gt.replace("|", "/").split("/")
    return any(allele not in {"0", ".", ""} for allele in alleles)


def sample_passes_qc(sample: dict[str, str | int | None], min_dp: int, min_gq: int) -> bool:
    dp = sample.get("DP")
    gq = sample.get("GQ")
    return isinstance(dp, int) and isinstance(gq, int) and dp >= min_dp and gq >= min_gq


def sorted_join(values: Iterable[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def segregation_flag(paternal_count: int, maternal_count: int) -> str:
    if paternal_count > 0 and maternal_count > 0:
        return "BIPARENTAL_HIGH_EFFECT_LURE"
    if paternal_count > 0:
        return "PATERNAL_ONLY_HIGH_EFFECT"
    if maternal_count > 0:
        return "MATERNAL_ONLY_HIGH_EFFECT"
    return "NO_HIGH_EFFECT_GENE_EXTRACTED"


def parse_vcf(
    input_path: Path,
    assembly: str,
    case_sample: str,
    father_sample: str,
    mother_sample: str,
    sibling_sample: str | None,
    min_dp: int,
    min_gq: int,
) -> tuple[list[dict[str, str]], dict[str, int], list[str]]:
    sample_names: list[str] = []
    sample_indices: dict[str, int] = {}
    variant_rows: list[dict[str, str]] = []
    record_counts = {
        "input_records": 0,
        "paternal_records": 0,
        "maternal_records": 0,
        "ambiguous_records": 0,
        "records_with_high_effect_annotations": 0,
    }

    with open_text(input_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                sample_names = fields[9:]
                for sample in [case_sample, father_sample, mother_sample]:
                    if sample not in sample_names:
                        raise ValueError(f"Sample {sample!r} not present in VCF header")
                if sibling_sample and sibling_sample not in sample_names:
                    raise ValueError(f"Sibling sample {sibling_sample!r} not present in VCF header")
                sample_indices = {name: 9 + sample_names.index(name) for name in sample_names}
                continue
            if line.startswith("#"):
                continue
            if not sample_indices:
                raise ValueError("VCF header with sample names was not found before records")

            record_counts["input_records"] += 1
            fields = line.split("\t")
            if len(fields) < 10:
                raise ValueError(f"VCF record has too few columns: {line[:120]}")
            chrom, pos, variant_id, ref, alt, _qual, _filter, info, fmt = fields[:9]
            format_keys = fmt.split(":")
            case = parse_sample(format_keys, fields[sample_indices[case_sample]])
            father = parse_sample(format_keys, fields[sample_indices[father_sample]])
            mother = parse_sample(format_keys, fields[sample_indices[mother_sample]])
            sibling = parse_sample(format_keys, fields[sample_indices[sibling_sample]]) if sibling_sample else {"GT": "", "DP": None, "GQ": None}

            case_alt = carries_alt(str(case["GT"]))
            father_alt = carries_alt(str(father["GT"]))
            mother_alt = carries_alt(str(mother["GT"]))
            passes_carrier_logic = case_alt and ((1 if father_alt else 0) + (1 if mother_alt else 0) == 1)
            passes_dp_gq = all(
                sample_passes_qc(sample, min_dp, min_gq)
                for sample in [case, father, mother] + ([sibling] if sibling_sample else [])
            )

            if passes_carrier_logic and father_alt:
                side = "paternal"
                record_counts["paternal_records"] += 1
            elif passes_carrier_logic and mother_alt:
                side = "maternal"
                record_counts["maternal_records"] += 1
            else:
                side = "ambiguous"
                record_counts["ambiguous_records"] += 1

            variant_key = f"{assembly}:{chrom}:{pos}:{ref}:{alt}"
            locus = f"{chrom}:{pos} {ref}>{alt}"
            high_annotations = [ann for ann in parse_eff(info) if ann["impact"] == "HIGH"]
            if high_annotations:
                record_counts["records_with_high_effect_annotations"] += 1

            for ann in high_annotations:
                variant_rows.append(
                    {
                        "variant_key": variant_key,
                        "chrom": chrom,
                        "pos": pos,
                        "id": variant_id,
                        "ref": ref,
                        "alt": alt,
                        "gene": ann["gene"],
                        "transcript": ann["transcript"],
                        "effect": ann["effect"],
                        "effect_impact": ann["impact"],
                        "segregation_side": side,
                        "case_sample": case_sample,
                        "father_sample": father_sample,
                        "mother_sample": mother_sample,
                        "sibling_sample": sibling_sample or "",
                        "case_gt": str(case["GT"] or ""),
                        "father_gt": str(father["GT"] or ""),
                        "mother_gt": str(mother["GT"] or ""),
                        "sibling_gt": str(sibling["GT"] or ""),
                        "case_dp": str(case["DP"] or ""),
                        "father_dp": str(father["DP"] or ""),
                        "mother_dp": str(mother["DP"] or ""),
                        "sibling_dp": str(sibling["DP"] or ""),
                        "case_gq": str(case["GQ"] or ""),
                        "father_gq": str(father["GQ"] or ""),
                        "mother_gq": str(mother["GQ"] or ""),
                        "sibling_gq": str(sibling["GQ"] or ""),
                        "passes_dp_gq": str(passes_dp_gq).lower(),
                        "passes_carrier_logic": str(passes_carrier_logic).lower(),
                        "gene_both_parental_sides_high_effect": "false",
                        "segregation_flag": "",
                        "_locus": locus,
                    }
                )
    return variant_rows, record_counts, sample_names


def build_gene_rows(variant_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    genes: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in variant_rows:
        gene = row["gene"]
        side = row["segregation_side"]
        if side not in {"paternal", "maternal"}:
            continue
        genes[gene][f"{side}_variant_keys"].add(row["variant_key"])
        genes[gene][f"{side}_variant_ids"].add(row["id"])
        genes[gene][f"{side}_loci"].add(row["_locus"])
        genes[gene][f"{side}_effects"].add(row["effect"])

    gene_flags: dict[str, str] = {}
    gene_rows: list[dict[str, str]] = []
    for gene, data in genes.items():
        paternal_count = len(data["paternal_variant_keys"])
        maternal_count = len(data["maternal_variant_keys"])
        total = paternal_count + maternal_count
        both = paternal_count > 0 and maternal_count > 0
        flag = segregation_flag(paternal_count, maternal_count)
        gene_flags[gene] = flag
        gene_rows.append(
            {
                "gene": gene,
                "both_parental_sides_high_effect": str(both).lower(),
                "segregation_flag": flag,
                "paternal_high_effect_variant_count": str(paternal_count),
                "maternal_high_effect_variant_count": str(maternal_count),
                "total_high_effect_variant_count": str(total),
                "paternal_variant_ids": sorted_join(data["paternal_variant_ids"]),
                "maternal_variant_ids": sorted_join(data["maternal_variant_ids"]),
                "paternal_loci": sorted_join(data["paternal_loci"]),
                "maternal_loci": sorted_join(data["maternal_loci"]),
                "paternal_effects": sorted_join(data["paternal_effects"]),
                "maternal_effects": sorted_join(data["maternal_effects"]),
                "interpretation_boundary": "NO_DIAGNOSIS; segregation side is not molecular phase; no rarity or pathogenicity inferred",
            }
        )

    for row in variant_rows:
        flag = gene_flags.get(row["gene"], "NO_HIGH_EFFECT_GENE_EXTRACTED")
        row["segregation_flag"] = flag
        row["gene_both_parental_sides_high_effect"] = str(flag == "BIPARENTAL_HIGH_EFFECT_LURE").lower()

    gene_rows.sort(
        key=lambda row: (
            row["both_parental_sides_high_effect"] != "true",
            -int(row["total_high_effect_variant_count"]),
            row["gene"],
        )
    )
    return gene_rows


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    input_path = Path(args.input).resolve()
    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(input_path)
    if not (input_path.name.endswith(".vcf") or input_path.name.endswith(".vcf.gz")):
        raise ValueError("ve-segregation currently expects .vcf or .vcf.gz input")

    variant_rows, record_counts, sample_names = parse_vcf(
        input_path=input_path,
        assembly=args.assembly,
        case_sample=args.case_sample,
        father_sample=args.father_sample,
        mother_sample=args.mother_sample,
        sibling_sample=args.sibling_sample,
        min_dp=args.min_dp,
        min_gq=args.min_gq,
    )
    gene_rows = build_gene_rows(variant_rows)

    variant_path = output_dir / "variant_segregation.tsv"
    gene_path = output_dir / "gene_segregation.tsv"
    result_path = output_dir / "result.json"
    write_tsv(variant_path, variant_rows, VARIANT_COLUMNS)
    write_tsv(gene_path, gene_rows, GENE_COLUMNS)

    flagged_genes = [row["gene"] for row in gene_rows if row["both_parental_sides_high_effect"] == "true"]
    summary = {
        **record_counts,
        "high_effect_annotation_rows": len(variant_rows),
        "genes_total": len(gene_rows),
        "genes_paternal_only": sum(1 for row in gene_rows if row["segregation_flag"] == "PATERNAL_ONLY_HIGH_EFFECT"),
        "genes_maternal_only": sum(1 for row in gene_rows if row["segregation_flag"] == "MATERNAL_ONLY_HIGH_EFFECT"),
        "genes_both_parental_sides": len(flagged_genes),
    }
    result: dict[str, object] = {
        "skill": "ve-segregation",
        "version": VERSION,
        "input": {
            "path": str(input_path),
            "format": "vcf.gz" if input_path.name.endswith(".vcf.gz") else "vcf",
            "sha256": sha256(input_path),
            "assembly": args.assembly,
        },
        "parameters": {
            "min_dp": args.min_dp,
            "min_gq": args.min_gq,
        },
        "samples": {
            "case": args.case_sample,
            "father": args.father_sample,
            "mother": args.mother_sample,
            "sibling": args.sibling_sample,
            "vcf_samples": sample_names,
        },
        "summary": summary,
        "flagged_genes": flagged_genes,
        "outputs": {
            "variant_segregation": "variant_segregation.tsv",
            "gene_segregation": "gene_segregation.tsv",
        },
        "interpretation_boundary": {
            "diagnosis": "NO_DIAGNOSIS",
            "notes": [
                "Segregation side is not molecular phase.",
                "Biparental high-effect flag is not a compound-heterozygous call.",
                "Historical SnpEff EFF is not modern effect validation.",
                "No rarity, pathogenicity, or diagnosis is inferred.",
            ],
        },
    }
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Parent-side segregation for filtered pedigree VCFs")
    parser.add_argument("--input", required=True, help="Input .vcf or .vcf.gz")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--case-sample", default=DEFAULT_CASE)
    parser.add_argument("--father-sample", default=DEFAULT_FATHER)
    parser.add_argument("--mother-sample", default=DEFAULT_MOTHER)
    parser.add_argument("--sibling-sample", default=DEFAULT_SIBLING)
    parser.add_argument("--min-dp", type=int, default=10)
    parser.add_argument("--min-gq", type=int, default=20)
    parser.add_argument("--assembly", default=DEFAULT_ASSEMBLY)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run(args)


if __name__ == "__main__":
    main()
