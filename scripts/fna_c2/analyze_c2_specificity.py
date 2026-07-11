#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from analyze_crc_association import (
    COHORT_ORDER,
    fit_ols,
    paper_standardized_effect,
    paule_mandel_meta,
    write_tsv,
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _index_tool_rows(
    rows: list[dict[str, str]],
    *,
    tool: str,
    label: str,
) -> dict[str, dict[str, str]]:
    selected = [row for row in rows if row.get("tool") == tool]
    indexed: dict[str, dict[str, str]] = {}
    for row in selected:
        sample = row.get("sample", "")
        if not sample:
            raise ValueError(f"{label} contains an empty sample")
        if sample in indexed:
            raise ValueError(f"duplicate {label} row for {tool}/{sample}")
        indexed[sample] = row
    return indexed


def prepare_rows(
    manifest_rows: list[dict[str, str]],
    audit_rows: list[dict[str, str]],
    signal_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    manifest_by_sample: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        sample = row.get("sample", "")
        if not sample:
            raise ValueError("sample manifest contains an empty sample")
        if sample in manifest_by_sample:
            raise ValueError(f"duplicate sample manifest row: {sample}")
        manifest_by_sample[sample] = row

    audit_by_sample = _index_tool_rows(audit_rows, tool="Chimera", label="audit")
    signal_by_sample = _index_tool_rows(signal_rows, tool="Chimera", label="signal")
    expected_samples = set(manifest_by_sample)
    for label, observed in (("audit", set(audit_by_sample)), ("signal", set(signal_by_sample))):
        if observed != expected_samples:
            raise ValueError(
                f"Chimera {label} sample set mismatch: "
                f"missing={len(expected_samples - observed)} extra={len(observed - expected_samples)}"
            )

    prepared: list[dict[str, object]] = []
    for sample, manifest in manifest_by_sample.items():
        audit = audit_by_sample[sample]
        signal = signal_by_sample[sample]
        for field in ("cohort", "condition", "role"):
            values = {manifest.get(field, ""), audit.get(field, ""), signal.get(field, "")}
            if len(values) != 1 or "" in values:
                raise ValueError(f"{field} mismatch for {sample}: {sorted(values)}")
        if signal.get("signal_unit") != "reads_per_million":
            raise ValueError(f"Chimera signal must use reads_per_million for {sample}")

        depth = int(manifest["input_reads_used"])
        if depth <= 0 or int(audit["input_reads_used"]) != depth or int(signal["input_reads_used"]) != depth:
            raise ValueError(f"input depth mismatch for {sample}")
        strict_reads = int(audit["strict_c2_best_reads"])
        if not 0 <= strict_reads <= int(audit["candidate_reads"]):
            raise ValueError(f"invalid strict C2 count for {sample}")

        condition = manifest["condition"].lower()
        sex = manifest["sex"].lower()
        if condition not in {"disease", "control"}:
            raise ValueError(f"unsupported condition for {sample}: {condition}")
        if sex not in {"female", "male"}:
            raise ValueError(f"unsupported sex for {sample}: {sex}")

        paper_c1_fraction = float(manifest["paper_fna_c1_percent"]) / 100.0
        chimera_non_c2_fraction = (
            float(signal["fna_c1_signal"]) + float(signal["non_c1c2_signal"])
        ) / 1_000_000.0
        prepared.append(
            {
                **manifest,
                "case": int(condition == "disease"),
                "male": int(sex == "male"),
                "response": asin_sqrt_fraction(strict_reads / depth),
                "published_c1": asin_sqrt_fraction(paper_c1_fraction),
                "published_fna_c1_percent": float(manifest["paper_fna_c1_percent"]),
                "published_fna_c2_percent": float(manifest["paper_fna_c2_percent"]),
                "chimera_non_c2": asin_sqrt_fraction(chimera_non_c2_fraction),
            }
        )
    if len(prepared) != 760:
        raise ValueError(f"expected 760 frozen samples, found {len(prepared)}")
    if {str(row["cohort"]) for row in prepared} != set(COHORT_ORDER):
        raise ValueError("sample manifest does not contain the frozen three cohorts")
    return prepared


def asin_sqrt_fraction(value: float) -> float:
    return math.asin(math.sqrt(min(1.0, max(0.0, value))))


def zscore(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    if variance <= 0.0:
        raise ValueError("cannot standardize a constant variable")
    scale = math.sqrt(variance)
    return [(value - mean) / scale for value in values]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Test whether frozen Chimera strict evidence follows published Fna C2 rather than C1."
    )
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--audit-table", required=True, type=Path)
    parser.add_argument("--signal-table", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    prepared = prepare_rows(
        read_tsv(args.sample_manifest),
        read_tsv(args.audit_table),
        read_tsv(args.signal_table),
    )

    models = {
        "base": (),
        "plus_published_fna_c1": ("published_c1",),
        "plus_chimera_c1_and_non_c1c2": ("chimera_non_c2",),
        "plus_both_non_c2_covariates": ("published_c1", "chimera_non_c2"),
    }
    cohort_rows: list[dict[str, object]] = []
    meta_rows: list[dict[str, object]] = []
    for model_name, extra_fields in models.items():
        effects: list[float] = []
        variances: list[float] = []
        for cohort in COHORT_ORDER:
            group = [row for row in prepared if row["cohort"] == cohort]
            design = [
                [
                    1.0,
                    float(row["case"]),
                    float(row["male"]),
                    float(row["age"]),
                    float(row["bmi"]),
                    float(row["source_reads"]) / 10_000_000.0,
                    *(float(row[field]) for field in extra_fields),
                ]
                for row in group
            ]
            fitted = fit_ols(
                [float(row["response"]) for row in group],
                design,
                coefficient_index=1,
            )
            n_case = sum(int(row["case"]) for row in group)
            n_control = len(group) - n_case
            effect, effect_se = paper_standardized_effect(
                statistic=float(fitted["statistic"]),
                n_case=n_case,
                n_control=n_control,
            )
            effects.append(effect)
            variances.append(effect_se * effect_se)
            cohort_rows.append(
                {
                    "model": model_name,
                    "cohort": cohort,
                    "n_case": n_case,
                    "n_control": n_control,
                    "effect": f"{effect:.6f}",
                    "effect_se": f"{effect_se:.6f}",
                    "condition_beta": f"{float(fitted['coefficient']):.6f}",
                    "condition_p_value": f"{float(fitted['p_value']):.6g}",
                }
            )
        meta = paule_mandel_meta(effects, variances)
        meta_rows.append(
            {
                "model": model_name,
                "cohorts": 3,
                "samples": len(prepared),
                "effect": f"{float(meta['effect']):.6f}",
                "ci_low": f"{float(meta['ci_low']):.6f}",
                "ci_high": f"{float(meta['ci_high']):.6f}",
                "p_value": f"{float(meta['wald_p_value']):.6g}",
                "tau_squared": f"{float(meta['tau_squared']):.6f}",
                "i_squared": f"{float(meta['i_squared']):.6f}",
            }
        )

    abundance_cohort_rows: list[dict[str, object]] = []
    abundance_effects = {"published_fna_c2": [], "published_fna_c1": []}
    abundance_variances = {"published_fna_c2": [], "published_fna_c1": []}
    for cohort in COHORT_ORDER:
        group = [row for row in prepared if row["cohort"] == cohort]
        response = zscore([float(row["response"]) for row in group])
        c2 = zscore(
            [
                asin_sqrt_fraction(float(row["published_fna_c2_percent"]) / 100.0)
                for row in group
            ]
        )
        c1 = zscore([float(row["published_c1"]) for row in group])
        design = [
            [
                1.0,
                c2_value,
                c1_value,
                float(row["case"]),
                float(row["male"]),
                float(row["age"]),
                float(row["bmi"]),
                float(row["source_reads"]) / 10_000_000.0,
            ]
            for row, c2_value, c1_value in zip(group, c2, c1)
        ]
        for predictor, coefficient_index in (
            ("published_fna_c2", 1),
            ("published_fna_c1", 2),
        ):
            fitted = fit_ols(response, design, coefficient_index=coefficient_index)
            effect = float(fitted["coefficient"])
            standard_error = float(fitted["standard_error"])
            abundance_effects[predictor].append(effect)
            abundance_variances[predictor].append(standard_error * standard_error)
            abundance_cohort_rows.append(
                {
                    "cohort": cohort,
                    "predictor": predictor,
                    "samples": len(group),
                    "standardized_beta": f"{effect:.6f}",
                    "standard_error": f"{standard_error:.6f}",
                    "p_value": f"{float(fitted['p_value']):.6g}",
                }
            )

    abundance_meta_rows: list[dict[str, object]] = []
    for predictor in ("published_fna_c2", "published_fna_c1"):
        meta = paule_mandel_meta(
            abundance_effects[predictor], abundance_variances[predictor]
        )
        abundance_meta_rows.append(
            {
                "predictor": predictor,
                "cohorts": 3,
                "samples": len(prepared),
                "standardized_beta": f"{float(meta['effect']):.6f}",
                "ci_low": f"{float(meta['ci_low']):.6f}",
                "ci_high": f"{float(meta['ci_high']):.6f}",
                "p_value": f"{float(meta['wald_p_value']):.6g}",
                "tau_squared": f"{float(meta['tau_squared']):.6f}",
                "i_squared": f"{float(meta['i_squared']):.6f}",
            }
        )

    baseline = next(row for row in meta_rows if row["model"] == "base")
    if abs(float(baseline["effect"]) - 0.336435) > 5e-7:
        raise AssertionError(f"base model did not reproduce primary effect: {baseline}")

    write_tsv(
        args.out_dir / "c2_specificity_adjusted_cohort.tsv",
        cohort_rows,
        [
            "model",
            "cohort",
            "n_case",
            "n_control",
            "effect",
            "effect_se",
            "condition_beta",
            "condition_p_value",
        ],
    )
    write_tsv(
        args.out_dir / "c2_specificity_adjusted_meta.tsv",
        meta_rows,
        [
            "model",
            "cohorts",
            "samples",
            "effect",
            "ci_low",
            "ci_high",
            "p_value",
            "tau_squared",
            "i_squared",
        ],
    )
    write_tsv(
        args.out_dir / "strict_c2_vs_paper_clades_cohort.tsv",
        abundance_cohort_rows,
        [
            "cohort",
            "predictor",
            "samples",
            "standardized_beta",
            "standard_error",
            "p_value",
        ],
    )
    write_tsv(
        args.out_dir / "strict_c2_vs_paper_clades_meta.tsv",
        abundance_meta_rows,
        [
            "predictor",
            "cohorts",
            "samples",
            "standardized_beta",
            "ci_low",
            "ci_high",
            "p_value",
            "tau_squared",
            "i_squared",
        ],
    )
    combined_abundance_rows: list[dict[str, object]] = []
    for row in abundance_cohort_rows:
        effect = float(row["standardized_beta"])
        standard_error = float(row["standard_error"])
        combined_abundance_rows.append(
            {
                "scope": "cohort",
                "cohort": row["cohort"],
                "predictor": row["predictor"],
                "samples": row["samples"],
                "standardized_beta": row["standardized_beta"],
                "ci_low": f"{effect - 1.959963984540054 * standard_error:.6f}",
                "ci_high": f"{effect + 1.959963984540054 * standard_error:.6f}",
                "p_value": row["p_value"],
                "tau_squared": "NA",
                "i_squared": "NA",
            }
        )
    for row in abundance_meta_rows:
        combined_abundance_rows.append(
            {
                "scope": "three_cohort_meta",
                "cohort": "all",
                "predictor": row["predictor"],
                "samples": row["samples"],
                "standardized_beta": row["standardized_beta"],
                "ci_low": row["ci_low"],
                "ci_high": row["ci_high"],
                "p_value": row["p_value"],
                "tau_squared": row["tau_squared"],
                "i_squared": row["i_squared"],
            }
        )
    write_tsv(
        args.out_dir / "c2_specificity_partial_association.tsv",
        combined_abundance_rows,
        [
            "scope",
            "cohort",
            "predictor",
            "samples",
            "standardized_beta",
            "ci_low",
            "ci_high",
            "p_value",
            "tau_squared",
            "i_squared",
        ],
    )
    print(f"wrote C2 specificity sensitivity analysis to {args.out_dir}")


if __name__ == "__main__":
    main()
