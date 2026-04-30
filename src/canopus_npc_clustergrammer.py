#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats


FIXED_NPC_PATHWAY_FAMILIES: Dict[str, str] = {
    "Terpenoids": "micro_cvd_purple",
    "Fatty acids": "micro_cvd_blue",
    "Polyketides": "micro_cvd_orange",
    "Alkaloids": "micro_cvd_green",
    "Shikimates and Phenylpropanoids": "micro_cvd_turquoise",
    "Amino acids and Peptides": "micro_orange",
    "Carbohydrates": "micro_purple",
    "Other": "micro_cvd_gray",
    "Unclassified": "micro_cvd_gray",
}

FIXED_NPC_SHADES: Dict[str, List[str]] = {
    "micro_cvd_gray": ["#616161", "#8B8B8B", "#B7B7B7", "#D6D6D6", "#F5F5F5"],
    "micro_cvd_purple": ["#7D3560", "#A1527F", "#CC79A7", "#E794C1", "#EFB6D6"],
    "micro_cvd_blue": ["#098BD9", "#56B4E9", "#7DCCFF", "#BCE1FF", "#E7F4FF"],
    "micro_cvd_orange": ["#9D654C", "#C17754", "#F09163", "#FCB076", "#FFD5AF"],
    "micro_cvd_green": ["#4E7705", "#6D9F06", "#97CE2F", "#BDEC6F", "#DDFFA0"],
    "micro_cvd_turquoise": ["#148F77", "#009E73", "#43BA8F", "#48C9B0", "#A3E4D7"],
    "micro_orange": ["#ff7f00", "#fe9929", "#fdae6b", "#fec44f", "#feeda0"],
    "micro_purple": ["#6a51a3", "#807dba", "#9e9ac8", "#bcbddc", "#dadaeb"],
}

DEFAULT_SAMPLE_COLORS = [
    "#4E79A7",
    "#E15759",
    "#76B7B2",
    "#F28E2B",
    "#59A14F",
    "#EDC948",
    "#B07AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
]

DEFAULT_CLUSTERGRAMMER_JS_URL = (
    "https://cdn.jsdelivr.net/gh/MaayanLab/clustergrammer@master/clustergrammer.min.js"
)
DEFAULT_JQUERY_URL = "https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js"
DEFAULT_UNDERSCORE_URL = "https://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.13.7/underscore-min.js"
DEFAULT_D3_URL = "https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.17/d3.min.js"
DEFAULT_BOOTSTRAP_CSS_URL = (
    "https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css"
)
DEFAULT_BOOTSTRAP_JS_URL = (
    "https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a Clustergrammer interactive heatmap from CANOPUS NPC annotations "
            "and the paired MZmine quant table."
        )
    )
    parser.add_argument(
        "--canopus-file",
        "-c",
        required=True,
        help="Path to canopus_structure_summary.tsv",
    )
    parser.add_argument(
        "--quant-file",
        "-q",
        default=None,
        help="Optional override for the MZmine quant CSV",
    )
    parser.add_argument(
        "--metadata-file",
        "-m",
        default=None,
        help="Optional override for the treated metadata TSV",
    )
    parser.add_argument(
        "--output-html",
        "-o",
        default=None,
        help="Output HTML file [default inferred next to the CANOPUS file]",
    )
    parser.add_argument(
        "--output-table",
        default=None,
        help="Optional TSV with the plotted matrix and annotations [default inferred]",
    )
    parser.add_argument(
        "--sample-type",
        default="sample",
        help="Metadata sample_type value to keep; use 'all' to keep every column [default: sample]",
    )
    parser.add_argument(
        "--sample-annotation",
        default="ATTRIBUTE_condition",
        help="Metadata column used for the main sample category [default: ATTRIBUTE_condition]",
    )
    parser.add_argument(
        "--preset",
        choices=["default", "metabolomics_sota", "metabolomics_sota_feature"],
        default="default",
        help="Apply a curated metabolomics heatmap preset [default: default]",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=300,
        help="Top N annotated features ranked by total intensity; use 0 for all [default: 300]",
    )
    parser.add_argument(
        "--top-by",
        choices=["discriminant", "intensity", "variance"],
        default="discriminant",
        help="Feature ranking used before --top-n filtering [default: discriminant]",
    )
    parser.add_argument(
        "--aggregate-by",
        choices=["none", "class", "superclass", "pathway"],
        default="none",
        help="Aggregate normalized intensities by NPC family before ranking/clustering [default: none]",
    )
    parser.add_argument(
        "--aggregate-method",
        choices=["sum", "mean", "median"],
        default="sum",
        help="How to aggregate family intensities when --aggregate-by is enabled [default: sum]",
    )
    parser.add_argument(
        "--transform",
        choices=["log10", "none"],
        default="log10",
        help="Intensity transform [default: log10]",
    )
    parser.add_argument(
        "--sample-normalization",
        choices=["median_sum", "none"],
        default="median_sum",
        help=(
            "Sample-wise normalization applied before filtering/ranking/clustering. "
            "'median_sum' rescales each sample to the median total signal across samples "
            "using all quantified features [default: median_sum]"
        ),
    )
    parser.add_argument(
        "--scale",
        choices=["row_zscore", "pareto", "none"],
        default="row_zscore",
        help="Matrix scaling [default: row_zscore]",
    )
    parser.add_argument(
        "--cluster-rows",
        choices=["true", "false"],
        default="true",
        help="Run hierarchical clustering on rows [default: true]",
    )
    parser.add_argument(
        "--cluster-cols",
        choices=["true", "false"],
        default="true",
        help="Run hierarchical clustering on columns [default: true]",
    )
    parser.add_argument(
        "--dist-type",
        default="cosine",
        help="Distance metric passed to Clustergrammer clustering [default: cosine]",
    )
    parser.add_argument(
        "--linkage-type",
        default="average",
        help="Linkage method passed to Clustergrammer clustering [default: average]",
    )
    parser.add_argument(
        "--npc-prob-threshold",
        type=float,
        default=0.85,
        help=(
            "Minimum CANOPUS probability required for pathway, superclass, and class. "
            "Use 0 to disable this filter [default: 0.85]"
        ),
    )
    parser.add_argument(
        "--min-effect-size",
        type=float,
        default=0.0,
        help="Minimum eta-squared effect size required after ranking-stat computation [default: 0]",
    )
    parser.add_argument(
        "--min-abs-log2-fc",
        type=float,
        default=0.0,
        help="Minimum maximum absolute pairwise log2 fold-change between groups [default: 0]",
    )
    parser.add_argument(
        "--max-pvalue",
        type=float,
        default=1.0,
        help="Maximum allowed ANOVA p-value before clustering [default: 1]",
    )
    parser.add_argument(
        "--max-fdr",
        type=float,
        default=1.0,
        help="Maximum allowed Benjamini-Hochberg FDR before clustering [default: 1]",
    )
    parser.add_argument(
        "--clustergrammer-js-url",
        default=DEFAULT_CLUSTERGRAMMER_JS_URL,
        help="Clustergrammer JS asset URL [default: jsDelivr GitHub mirror]",
    )
    parser.add_argument(
        "--jquery-url",
        default=DEFAULT_JQUERY_URL,
        help="jQuery asset URL [default: cdnjs]",
    )
    parser.add_argument(
        "--underscore-url",
        default=DEFAULT_UNDERSCORE_URL,
        help="Underscore asset URL [default: cdnjs]",
    )
    parser.add_argument(
        "--d3-url",
        default=DEFAULT_D3_URL,
        help="D3 v3 asset URL [default: cdnjs]",
    )
    parser.add_argument(
        "--bootstrap-css-url",
        default=DEFAULT_BOOTSTRAP_CSS_URL,
        help="Bootstrap CSS asset URL [default: cdnjs]",
    )
    parser.add_argument(
        "--bootstrap-js-url",
        default=DEFAULT_BOOTSTRAP_JS_URL,
        help="Bootstrap JS asset URL [default: cdnjs]",
    )
    return parser.parse_args()


def as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def infer_batch_dir(canopus_file: Path) -> Path:
    return canopus_file.parent.parent.parent


def infer_quant_file(canopus_file: Path) -> Path:
    batch_dir = infer_batch_dir(canopus_file)
    return batch_dir / "results" / "mzmine" / f"{batch_dir.name}_quant.csv"


def infer_metadata_file(canopus_file: Path) -> Path:
    batch_dir = infer_batch_dir(canopus_file)
    return batch_dir / "metadata" / "treated" / f"{batch_dir.name}_metadata.tsv"


def infer_output_html(canopus_file: Path) -> Path:
    return canopus_file.parent / "canopus_npc_clustergrammer.html"


def default_output_table(output_html: Path) -> Path:
    return output_html.with_name(f"{output_html.stem}_data.tsv")


def ensure_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def clean_intensity_name(column_name: str) -> str:
    suffix = " Peak height"
    return column_name[: -len(suffix)] if column_name.endswith(suffix) else column_name


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    means = df.mean(axis=1)
    stds = df.std(axis=1, ddof=1).replace(0, np.nan)
    scaled = df.sub(means, axis=0).div(stds, axis=0)
    return scaled.replace([np.inf, -np.inf], np.nan).fillna(0.0)


def pareto_scale_rows(df: pd.DataFrame) -> pd.DataFrame:
    means = df.mean(axis=1)
    stds = df.std(axis=1, ddof=1).replace(0, np.nan)
    scaled = df.sub(means, axis=0).div(np.sqrt(stds), axis=0)
    return scaled.replace([np.inf, -np.inf], np.nan).fillna(0.0)


def scale_rows(df: pd.DataFrame, method: str) -> pd.DataFrame:
    if method == "none":
        return df.copy()
    if method == "row_zscore":
        return zscore_rows(df)
    if method == "pareto":
        return pareto_scale_rows(df)
    raise ValueError(f"Unsupported row scaling method: {method}")


def apply_preset(args: argparse.Namespace) -> argparse.Namespace:
    if args.preset == "default":
        return args

    if args.preset == "metabolomics_sota":
        args.aggregate_by = "superclass"
        args.aggregate_method = "sum"
        args.sample_normalization = "median_sum"
        args.transform = "log10"
        args.scale = "pareto"
        args.top_by = "discriminant"
        args.npc_prob_threshold = max(args.npc_prob_threshold, 0.85)
        args.min_effect_size = max(args.min_effect_size, 0.05)
        args.min_abs_log2_fc = max(args.min_abs_log2_fc, 0.4)
        args.max_fdr = min(args.max_fdr, 0.7)
        args.dist_type = "correlation"
        args.linkage_type = "complete"
        args.cluster_rows = "true"
        args.cluster_cols = "true"
        if args.top_n == 300:
            args.top_n = 40
        return args

    if args.preset == "metabolomics_sota_feature":
        args.aggregate_by = "none"
        args.aggregate_method = "sum"
        args.sample_normalization = "median_sum"
        args.transform = "log10"
        args.scale = "pareto"
        args.top_by = "discriminant"
        args.npc_prob_threshold = max(args.npc_prob_threshold, 0.85)
        args.min_effect_size = max(args.min_effect_size, 0.12)
        args.min_abs_log2_fc = max(args.min_abs_log2_fc, 0.8)
        args.max_fdr = min(args.max_fdr, 0.2)
        args.dist_type = "correlation"
        args.linkage_type = "complete"
        args.cluster_rows = "true"
        args.cluster_cols = "true"
        if args.top_n == 300:
            args.top_n = 150
        return args

    raise ValueError(f"Unsupported preset: {args.preset}")


def normalize_samples_by_median_sum(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    col_sums = df.sum(axis=0)
    positive_sums = col_sums[col_sums > 0]
    if positive_sums.empty:
        return df.copy(), pd.Series(1.0, index=df.columns, dtype=float)
    target_sum = float(np.median(positive_sums.to_numpy(dtype=float)))
    scale_factors = pd.Series(1.0, index=df.columns, dtype=float)
    scale_factors.loc[positive_sums.index] = target_sum / positive_sums
    normalized = df.mul(scale_factors, axis=1)
    return normalized, scale_factors


def log_matrix_for_ranking(df: pd.DataFrame, transform: str) -> pd.DataFrame:
    if transform == "log10":
        return np.log10(df + 1.0)
    return df.copy()


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full_like(p_values, np.nan, dtype=float)
    valid = np.isfinite(p_values)
    if not np.any(valid):
        return adjusted
    valid_p = p_values[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    n = len(ranked)
    bh = ranked * n / np.arange(1, n + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    out_valid = np.empty_like(valid_p)
    out_valid[order] = bh
    adjusted[valid] = out_valid
    return adjusted


def compute_discriminant_stats(
    ranking_df: pd.DataFrame,
    groups: Sequence[str],
    raw_df: pd.DataFrame,
) -> pd.DataFrame:
    matrix = ranking_df.to_numpy(dtype=float)
    raw_matrix = raw_df.to_numpy(dtype=float)
    labels = pd.Series(groups, dtype="object").fillna("Missing")
    unique_groups = list(pd.unique(labels))

    if len(unique_groups) < 2 or matrix.shape[1] <= len(unique_groups):
        variance_score = np.nanvar(matrix, axis=1)
        return pd.DataFrame(
            {
                "discriminant_score": variance_score,
                "variance_score": variance_score,
                "effect_size_eta2": np.zeros(matrix.shape[0], dtype=float),
                "max_abs_log2_fc": np.zeros(matrix.shape[0], dtype=float),
                "p_value": np.ones(matrix.shape[0], dtype=float),
                "fdr": np.ones(matrix.shape[0], dtype=float),
            },
            index=ranking_df.index,
        )

    overall_mean = matrix.mean(axis=1)
    ss_between = np.zeros(matrix.shape[0], dtype=float)
    ss_within = np.zeros(matrix.shape[0], dtype=float)
    raw_group_means: list[np.ndarray] = []

    for group in unique_groups:
        mask = labels.to_numpy() == group
        group_matrix = matrix[:, mask]
        if group_matrix.size == 0:
            continue
        group_mean = group_matrix.mean(axis=1)
        ss_between += group_matrix.shape[1] * np.square(group_mean - overall_mean)
        ss_within += np.square(group_matrix - group_mean[:, None]).sum(axis=1)
        raw_group_means.append(raw_matrix[:, mask].mean(axis=1))

    df_between = len(unique_groups) - 1
    df_within = matrix.shape[1] - len(unique_groups)

    if df_between <= 0 or df_within <= 0:
        variance_score = np.nanvar(matrix, axis=1)
        return pd.DataFrame(
            {
                "discriminant_score": variance_score,
                "variance_score": variance_score,
                "effect_size_eta2": np.zeros(matrix.shape[0], dtype=float),
                "max_abs_log2_fc": np.zeros(matrix.shape[0], dtype=float),
                "p_value": np.ones(matrix.shape[0], dtype=float),
                "fdr": np.ones(matrix.shape[0], dtype=float),
            },
            index=ranking_df.index,
        )

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    scores = np.zeros(matrix.shape[0], dtype=float)
    positive_within = ms_within > 0
    scores[positive_within] = ms_between[positive_within] / ms_within[positive_within]
    infinite_mask = (~positive_within) & (ms_between > 0)
    if np.any(infinite_mask):
        finite_scores = scores[np.isfinite(scores)]
        replacement = finite_scores.max() * 1.05 if finite_scores.size else float(ms_between[infinite_mask].max())
        scores[infinite_mask] = replacement if replacement > 0 else 1.0

    total_ss = ss_between + ss_within
    effect_size = np.divide(
        ss_between,
        total_ss,
        out=np.zeros(matrix.shape[0], dtype=float),
        where=total_ss > 0,
    )

    p_values = scipy_stats.f.sf(scores, df_between, df_within)
    p_values = np.clip(p_values, 0, 1)
    fdr = benjamini_hochberg(p_values)

    max_abs_log2_fc = np.zeros(matrix.shape[0], dtype=float)
    if len(raw_group_means) >= 2:
        for i in range(len(raw_group_means)):
            for j in range(i + 1, len(raw_group_means)):
                fc = np.abs(np.log2((raw_group_means[i] + 1.0) / (raw_group_means[j] + 1.0)))
                max_abs_log2_fc = np.maximum(max_abs_log2_fc, fc)

    return pd.DataFrame(
        {
            "discriminant_score": scores,
            "variance_score": np.nanvar(matrix, axis=1),
            "effect_size_eta2": effect_size,
            "max_abs_log2_fc": max_abs_log2_fc,
            "p_value": p_values,
            "fdr": fdr,
        },
        index=ranking_df.index,
    )


def aggregate_feature_table(
    feature_table: pd.DataFrame,
    intensity_cols: Sequence[str],
    aggregate_by: str,
    aggregate_method: str,
) -> pd.DataFrame:
    if aggregate_by == "none":
        out = feature_table.copy()
        out["row_kind"] = "feature"
        out["row_id"] = out["feature_id"].astype(str)
        out["row_label"] = "Feature: F" + out["feature_id"].astype(str) + " | " + out["npc_class"].astype(str)
        out["n_features"] = 1
        return out

    group_cols_map = {
        "class": ["npc_pathway", "npc_superclass", "npc_class"],
        "superclass": ["npc_pathway", "npc_superclass"],
        "pathway": ["npc_pathway"],
    }
    group_cols = group_cols_map[aggregate_by]

    agg_func = aggregate_method
    grouped = feature_table.groupby(group_cols, dropna=False)
    summary = grouped[list(intensity_cols)].agg(agg_func).reset_index()
    summary["n_features"] = grouped.size().to_numpy()
    summary["total_intensity"] = summary[list(intensity_cols)].sum(axis=1)
    summary["npc_pathway_prob"] = grouped["npc_pathway_prob"].mean().to_numpy()
    summary["npc_superclass_prob"] = grouped["npc_superclass_prob"].mean().to_numpy()
    summary["npc_class_prob"] = grouped["npc_class_prob"].mean().to_numpy()
    summary["feature_id"] = pd.NA
    summary["aligned_feature_id"] = pd.NA
    summary["ion_mass"] = pd.NA
    summary["retention_time_min"] = pd.NA
    summary["row_kind"] = aggregate_by
    if aggregate_by == "class":
        summary["row_id"] = summary["npc_class"].astype(str)
        summary["row_label"] = (
            "Class: " + summary["npc_class"].astype(str) + " (n=" + summary["n_features"].astype(str) + ")"
        )
    elif aggregate_by == "superclass":
        summary["npc_class"] = "Aggregated"
        summary["row_id"] = summary["npc_superclass"].astype(str)
        summary["row_label"] = (
            "Superclass: " + summary["npc_superclass"].astype(str) + " (n=" + summary["n_features"].astype(str) + ")"
        )
    else:
        summary["npc_superclass"] = "Aggregated"
        summary["npc_class"] = "Aggregated"
        summary["row_id"] = summary["npc_pathway"].astype(str)
        summary["row_label"] = (
            "Pathway: " + summary["npc_pathway"].astype(str) + " (n=" + summary["n_features"].astype(str) + ")"
        )
    return summary


def pick_sample_colors(values: Sequence[str]) -> Dict[str, str]:
    unique_values = list(dict.fromkeys(v for v in values if pd.notna(v)))
    return {
        value: DEFAULT_SAMPLE_COLORS[idx % len(DEFAULT_SAMPLE_COLORS)]
        for idx, value in enumerate(unique_values)
    }


def palette_family_for_pathway(pathway: str) -> str:
    return FIXED_NPC_PATHWAY_FAMILIES.get(pathway, "micro_cvd_gray")


def build_fixed_npc_color_maps(feature_meta: pd.DataFrame) -> dict[str, dict[str, str]]:
    feature_meta = feature_meta.copy()
    feature_meta["npc_pathway"] = feature_meta["npc_pathway"].fillna("Unclassified")
    feature_meta["npc_superclass"] = feature_meta["npc_superclass"].fillna("Other")
    feature_meta["npc_class"] = feature_meta["npc_class"].fillna("Other")

    pathway_colors: dict[str, str] = {}
    for pathway in feature_meta["npc_pathway"].drop_duplicates():
        family = palette_family_for_pathway(pathway)
        pathway_colors[pathway] = FIXED_NPC_SHADES[family][0]

    rank_df = (
        feature_meta.groupby(["npc_pathway", "npc_superclass"], dropna=False)
        .size()
        .reset_index(name="abundance")
        .sort_values(["npc_pathway", "abundance", "npc_superclass"], ascending=[True, False, True])
    )
    rank_df["within_pathway_rank"] = rank_df.groupby("npc_pathway").cumcount() + 1
    rank_df["shade_index"] = rank_df["within_pathway_rank"].clip(upper=5)
    rank_df.loc[rank_df["within_pathway_rank"] > 4, "shade_index"] = 5
    rank_df["shade_index"] = rank_df["shade_index"].astype(int)
    rank_df["family"] = rank_df["npc_pathway"].map(palette_family_for_pathway)
    rank_df["hex"] = rank_df.apply(
        lambda row: FIXED_NPC_SHADES[row["family"]][row["shade_index"] - 1], axis=1
    )

    superclass_colors = {
        row["npc_superclass"]: row["hex"] for _, row in rank_df.iterrows()
    }

    class_df = feature_meta.merge(
        rank_df[["npc_pathway", "npc_superclass", "hex"]],
        on=["npc_pathway", "npc_superclass"],
        how="left",
    )
    class_df["hex"] = class_df["hex"].fillna(class_df["npc_pathway"].map(pathway_colors))
    class_colors = {
        row["npc_class"]: row["hex"]
        for _, row in class_df[["npc_class", "hex"]].drop_duplicates().iterrows()
    }

    return {
        "pathway": pathway_colors,
        "superclass": superclass_colors,
        "class": class_colors,
    }


def build_row_tuples(feature_meta: pd.DataFrame) -> List[tuple[str, ...]]:
    tuples: List[tuple[str, ...]] = []
    for _, row in feature_meta.iterrows():
        tuples.append(
            (
                str(row.row_label),
                f"NPC Pathway: {row.npc_pathway}",
                f"NPC Superclass: {row.npc_superclass}",
                f"NPC Class: {row.npc_class}",
            )
        )
    return tuples


def build_col_tuples(sample_manifest: pd.DataFrame) -> List[tuple[str, ...]]:
    include_sample_type = sample_manifest["sample_type"].nunique(dropna=True) > 1
    tuples: List[tuple[str, ...]] = []
    for _, row in sample_manifest.iterrows():
        labels = [
            f"Sample: {row.sample_label}",
            f"Group: {row.sample_annotation}",
        ]
        if include_sample_type:
            labels.append(f"Sample Type: {row.sample_type}")
        tuples.append(tuple(labels))
    return tuples


def export_table(
    feature_meta: pd.DataFrame,
    matrix_df: pd.DataFrame,
    sample_manifest: pd.DataFrame,
    output_table: Path,
) -> None:
    preferred_cols = [
        "row_kind",
        "row_id",
        "row_label",
        "feature_id",
        "aligned_feature_id",
        "n_features",
        "ion_mass",
        "retention_time_min",
        "npc_pathway",
        "npc_superclass",
        "npc_class",
        "npc_pathway_prob",
        "npc_superclass_prob",
        "npc_class_prob",
        "discriminant_score",
        "variance_score",
        "effect_size_eta2",
        "max_abs_log2_fc",
        "p_value",
        "fdr",
        "total_intensity",
    ]
    out_df = feature_meta[[col for col in preferred_cols if col in feature_meta.columns]].copy()
    for intensity_col, filename in zip(sample_manifest["intensity_col"], sample_manifest["filename"]):
        out_df[filename] = matrix_df[intensity_col].to_numpy()
    out_df.to_csv(output_table, sep="\t", index=False)


def render_html(
    network_json: str,
    output_html: Path,
    *,
    title: str,
    about_html: str,
    input_domain: float,
    js_url: str,
    jquery_url: str,
    underscore_url: str,
    d3_url: str,
    bootstrap_css_url: str,
    bootstrap_js_url: str,
    row_label: str = "Features",
    col_label: str = "Samples",
    row_order: str = "clust",
    col_order: str = "clust",
) -> None:
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape_html(title)}</title>
  <link rel="stylesheet" href="{escape_html(bootstrap_css_url)}">
  <style>
    html, body {{
      height: 100%;
      margin: 0;
      background: #f7f7f5;
      color: #1f1f1f;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    .page {{
      min-height: 100%;
      display: grid;
      grid-template-rows: auto 1fr;
    }}
    .header {{
      padding: 14px 18px 8px 18px;
      border-bottom: 1px solid #ddd;
      background: #fbfbf9;
    }}
    .header h1 {{
      margin: 0 0 6px 0;
      font-size: 22px;
      font-weight: 700;
    }}
    .header p {{
      margin: 0;
      color: #555;
      font-size: 13px;
    }}
    #clustergrammer_container {{
      width: 100%;
      height: calc(100vh - 78px);
    }}
  </style>
</head>
<body>
  <div class="page">
    <div class="header">
      <h1>{escape_html(title)}</h1>
      <p>Cluster rows/columns, zoom, search, reorder, crop, and use the category bars to browse NPC structure.</p>
    </div>
    <div id="clustergrammer_container"></div>
  </div>

  <script src="{escape_html(jquery_url)}"></script>
  <script src="{escape_html(underscore_url)}"></script>
  <script src="{escape_html(d3_url)}"></script>
  <script src="{escape_html(bootstrap_js_url)}"></script>
  <script src="{escape_html(js_url)}"></script>
  <script>
    const networkData = {network_json};
    const args = {{
      root: '#clustergrammer_container',
      network_data: networkData,
      row_label: {json.dumps(row_label)},
      col_label: {json.dumps(col_label)},
      row_order: {json.dumps(row_order)},
      col_order: {json.dumps(col_order)},
      ini_expand: true,
      sidebar_width: 280,
      input_domain: {input_domain:.6f},
      tile_colors: ['#b2182b', '#2166ac'],
      about: {json.dumps(about_html)}
    }};

    const cgm = Clustergrammer(args);
    window.addEventListener('resize', function () {{
      if (cgm && typeof cgm.resize_viz === 'function') {{
        cgm.resize_viz();
      }}
    }});
  </script>
</body>
</html>
"""
    output_html.write_text(html, encoding="utf-8")


def escape_html(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def compute_input_domain(df: pd.DataFrame) -> float:
    values = np.abs(df.to_numpy(dtype=float)).ravel()
    values = values[np.isfinite(values)]
    if not len(values):
        return 1.0
    domain = float(np.quantile(values, 0.95))
    return domain if domain > 0 else float(values.max()) if values.max() > 0 else 1.0


def maybe_import_clustergrammer() -> object:
    if not hasattr(pd.DataFrame, "ix"):
        pd.DataFrame.ix = property(lambda self: self.loc)  # type: ignore[attr-defined]
    if not hasattr(pd.Series, "ix"):
        pd.Series.ix = property(lambda self: self.loc)  # type: ignore[attr-defined]
    try:
        from clustergrammer import Network
    except ImportError as exc:
        raise SystemExit(
            "clustergrammer is required for this script. Install it with "
            "`pip install --upgrade clustergrammer`."
        ) from exc
    return Network


def main() -> None:
    args = parse_args()
    args = apply_preset(args)

    canopus_file = Path(args.canopus_file).expanduser().resolve()
    quant_file = Path(args.quant_file).expanduser().resolve() if args.quant_file else infer_quant_file(canopus_file)
    metadata_file = (
        Path(args.metadata_file).expanduser().resolve()
        if args.metadata_file
        else infer_metadata_file(canopus_file)
    )
    output_html = (
        Path(args.output_html).expanduser().resolve()
        if args.output_html
        else infer_output_html(canopus_file)
    )
    output_table = (
        Path(args.output_table).expanduser().resolve()
        if args.output_table
        else default_output_table(output_html)
    )

    ensure_exists(canopus_file, "CANOPUS file")
    ensure_exists(quant_file, "Quant file")
    ensure_exists(metadata_file, "Metadata file")

    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_table.parent.mkdir(parents=True, exist_ok=True)

    if args.top_n < 0:
        raise SystemExit("--top-n must be a non-negative integer.")

    canopus_df = (
        pd.read_csv(canopus_file, sep="\t", usecols=[
            "mappingFeatureId",
            "alignedFeatureId",
            "ionMass",
            "retentionTimeInMinutes",
            "NPC#pathway",
            "NPC#pathway Probability",
            "NPC#superclass",
            "NPC#superclass Probability",
            "NPC#class",
            "NPC#class Probability",
        ])
        .rename(
            columns={
                "mappingFeatureId": "feature_id",
                "alignedFeatureId": "aligned_feature_id",
                "ionMass": "ion_mass",
                "retentionTimeInMinutes": "retention_time_min",
                "NPC#pathway": "npc_pathway",
                "NPC#pathway Probability": "npc_pathway_prob",
                "NPC#superclass": "npc_superclass",
                "NPC#superclass Probability": "npc_superclass_prob",
                "NPC#class": "npc_class",
                "NPC#class Probability": "npc_class_prob",
            }
        )
    )
    canopus_df["feature_id"] = pd.to_numeric(canopus_df["feature_id"], errors="coerce").astype("Int64")
    canopus_df = (
        canopus_df.dropna(subset=["feature_id"])
        .drop_duplicates(subset=["feature_id"])
        .copy()
    )
    canopus_df["feature_id"] = canopus_df["feature_id"].astype(int)
    canopus_df["npc_pathway"] = canopus_df["npc_pathway"].fillna("Unclassified")
    canopus_df["npc_superclass"] = canopus_df["npc_superclass"].fillna("Other")
    canopus_df["npc_class"] = canopus_df["npc_class"].fillna("Other")
    for prob_col in ["npc_pathway_prob", "npc_superclass_prob", "npc_class_prob"]:
        canopus_df[prob_col] = pd.to_numeric(canopus_df[prob_col], errors="coerce")
    if args.npc_prob_threshold > 0:
        canopus_df = canopus_df[
            (canopus_df["npc_pathway_prob"] >= args.npc_prob_threshold)
            & (canopus_df["npc_superclass_prob"] >= args.npc_prob_threshold)
            & (canopus_df["npc_class_prob"] >= args.npc_prob_threshold)
        ].copy()
        if canopus_df.empty:
            raise SystemExit(
                "No CANOPUS features passed the NPC probability threshold. "
                "Lower --npc-prob-threshold or disable it with 0."
            )

    quant_df = pd.read_csv(quant_file)
    intensity_cols = [col for col in quant_df.columns if col.endswith(" Peak height")]
    if not intensity_cols:
        raise SystemExit("No MZmine intensity columns ending with ' Peak height' were found.")
    for col in intensity_cols:
        quant_df[col] = pd.to_numeric(quant_df[col], errors="coerce").fillna(0.0)
    if args.sample_normalization == "median_sum":
        normalized_matrix, sample_scale_factors = normalize_samples_by_median_sum(quant_df[intensity_cols])
        quant_df.loc[:, intensity_cols] = normalized_matrix
    else:
        sample_scale_factors = pd.Series(1.0, index=intensity_cols, dtype=float)

    metadata_df = pd.read_csv(metadata_file, sep="\t")
    required_metadata_cols = {"filename", "sample_type", args.sample_annotation}
    missing_metadata_cols = required_metadata_cols.difference(metadata_df.columns)
    if missing_metadata_cols:
        raise SystemExit(
            "Metadata file is missing required columns: "
            + ", ".join(sorted(missing_metadata_cols))
        )

    sample_manifest = pd.DataFrame({"intensity_col": intensity_cols})
    sample_manifest["filename"] = sample_manifest["intensity_col"].map(clean_intensity_name)
    sample_manifest = sample_manifest.merge(metadata_df, on="filename", how="left")
    if sample_manifest["sample_type"].isna().any():
        missing_filename = sample_manifest.loc[sample_manifest["sample_type"].isna(), "filename"].iloc[0]
        raise SystemExit(f"Metadata is missing a row for quant column filename: {missing_filename}")

    if args.sample_type.lower() not in {"all", "*"}:
        sample_manifest = sample_manifest[
            sample_manifest["sample_type"].astype(str).str.lower() == args.sample_type.lower()
        ].copy()
        if sample_manifest.empty:
            raise SystemExit(
                f"No columns remain after filtering sample_type == {args.sample_type!r}."
            )

    sample_manifest["sample_type"] = sample_manifest["sample_type"].fillna("Missing").astype(str)
    sample_manifest["sample_annotation"] = (
        sample_manifest[args.sample_annotation].fillna("Missing").astype(str)
    )
    sample_manifest["sample_label"] = sample_manifest.get("sample_id", sample_manifest["filename"]).fillna(
        sample_manifest["filename"]
    ).astype(str)
    sample_manifest = sample_manifest.sort_values(
        by=["sample_annotation", "sample_type", "sample_label"]
    ).reset_index(drop=True)
    if args.top_by == "discriminant" and sample_manifest["sample_annotation"].nunique(dropna=True) < 2:
        print(
            "Only one sample annotation group remains after filtering; "
            "discriminant ranking falls back to variance ranking.",
            file=sys.stderr,
        )

    matrix_df = quant_df[["row ID", *sample_manifest["intensity_col"].tolist()]].copy()
    matrix_df = matrix_df.rename(columns={"row ID": "feature_id"})
    matrix_df["feature_id"] = pd.to_numeric(matrix_df["feature_id"], errors="coerce")
    matrix_df = matrix_df.dropna(subset=["feature_id"]).copy()
    matrix_df["feature_id"] = matrix_df["feature_id"].astype(int)

    feature_table = canopus_df.merge(matrix_df, on="feature_id", how="inner")
    feature_table["total_intensity"] = feature_table[sample_manifest["intensity_col"]].sum(axis=1)
    analysis_table = aggregate_feature_table(
        feature_table,
        sample_manifest["intensity_col"].tolist(),
        args.aggregate_by,
        args.aggregate_method,
    )
    ranking_matrix = log_matrix_for_ranking(analysis_table[sample_manifest["intensity_col"]], args.transform)
    stats_df = compute_discriminant_stats(
        ranking_matrix,
        sample_manifest["sample_annotation"].tolist(),
        analysis_table[sample_manifest["intensity_col"]],
    )
    analysis_table = pd.concat([analysis_table.reset_index(drop=True), stats_df.reset_index(drop=True)], axis=1)
    analysis_table = analysis_table[
        (analysis_table["effect_size_eta2"] >= args.min_effect_size)
        & (analysis_table["max_abs_log2_fc"] >= args.min_abs_log2_fc)
        & (analysis_table["p_value"] <= args.max_pvalue)
        & (analysis_table["fdr"] <= args.max_fdr)
    ].copy()
    if analysis_table.empty:
        raise SystemExit(
            "No rows remained after the discriminant filters. "
            "Relax --min-effect-size, --min-abs-log2-fc, --max-pvalue, or --max-fdr."
        )

    if args.top_by == "discriminant":
        rank_column = "discriminant_score"
    elif args.top_by == "variance":
        rank_column = "variance_score"
    else:
        rank_column = "total_intensity"

    analysis_table = analysis_table.sort_values(
        by=[rank_column, "effect_size_eta2", "max_abs_log2_fc", "npc_pathway", "npc_superclass", "npc_class", "row_label"],
        ascending=[False, False, False, True, True, True, True],
    ).reset_index(drop=True)

    if args.top_n > 0:
        analysis_table = analysis_table.head(args.top_n).copy()
    if analysis_table.empty:
        raise SystemExit("No annotated features remained after joining CANOPUS and the quant table.")

    export_table(analysis_table, analysis_table, sample_manifest, output_table)

    matrix_for_cluster = analysis_table[sample_manifest["intensity_col"]].copy()
    if args.transform == "log10":
        matrix_for_cluster = np.log10(matrix_for_cluster + 1.0)
    matrix_for_cluster = scale_rows(matrix_for_cluster, args.scale)

    feature_meta = analysis_table[
        [
            "row_kind",
            "row_id",
            "row_label",
            "feature_id",
            "aligned_feature_id",
            "n_features",
            "ion_mass",
            "retention_time_min",
            "npc_pathway",
            "npc_superclass",
            "npc_class",
            "npc_pathway_prob",
            "npc_superclass_prob",
            "npc_class_prob",
            "discriminant_score",
            "variance_score",
            "effect_size_eta2",
            "max_abs_log2_fc",
            "p_value",
            "fdr",
            "total_intensity",
        ]
    ].copy()

    cluster_df = pd.DataFrame(
        matrix_for_cluster.to_numpy(dtype=float),
        index=build_row_tuples(feature_meta),
        columns=build_col_tuples(sample_manifest),
    )

    Network = maybe_import_clustergrammer()
    net = Network()
    net.load_df(cluster_df)

    npc_colors = build_fixed_npc_color_maps(feature_meta)
    for pathway, color in npc_colors["pathway"].items():
        net.set_cat_color("row", 1, pathway, color)
    for superclass, color in npc_colors["superclass"].items():
        net.set_cat_color("row", 2, superclass, color)
    for npc_class, color in npc_colors["class"].items():
        net.set_cat_color("row", 3, npc_class, color)

    sample_group_colors = pick_sample_colors(sample_manifest["sample_annotation"].tolist())
    for group_name, color in sample_group_colors.items():
        net.set_cat_color("col", 1, group_name, color)
    if sample_manifest["sample_type"].nunique(dropna=True) > 1:
        sample_type_colors = pick_sample_colors(sample_manifest["sample_type"].tolist())
        for sample_type, color in sample_type_colors.items():
            net.set_cat_color("col", 2, sample_type, color)

    net.cluster(
        dist_type=args.dist_type,
        linkage_type=args.linkage_type,
        views=["N_row_sum", "N_row_var"],
        run_clustering=True,
    )

    network_json = net.export_net_json("viz")
    title = "CANOPUS NPC Clustergrammer heatmap"
    about_html = (
        f"<div><strong>Preset:</strong> {args.preset}<br>"
        f"<strong>Aggregation:</strong> {args.aggregate_by} ({args.aggregate_method})<br>"
        f"<strong>Top-by:</strong> {args.top_by}<br>"
        f"<strong>Sample normalization:</strong> {args.sample_normalization}<br>"
        f"<strong>Transform:</strong> {args.transform}<br>"
        f"<strong>Scale:</strong> {args.scale}<br>"
        f"<strong>Distance:</strong> {args.dist_type}<br>"
        f"<strong>Linkage:</strong> {args.linkage_type}<br>"
        f"<strong>NPC probability threshold:</strong> {args.npc_prob_threshold}<br>"
        f"<strong>Min eta²:</strong> {args.min_effect_size}<br>"
        f"<strong>Min max |log2FC|:</strong> {args.min_abs_log2_fc}<br>"
        f"<strong>Max p-value:</strong> {args.max_pvalue}<br>"
        f"<strong>Max FDR:</strong> {args.max_fdr}<br>"
        f"<strong>Rows:</strong> {len(analysis_table)}<br>"
        f"<strong>Samples:</strong> {len(sample_manifest)}<br>"
        f"<strong>CANOPUS:</strong> {canopus_file.name}</div>"
    )
    render_html(
        network_json,
        output_html,
        title=title,
        about_html=about_html,
        input_domain=compute_input_domain(matrix_for_cluster),
        js_url=args.clustergrammer_js_url,
        jquery_url=args.jquery_url,
        underscore_url=args.underscore_url,
        d3_url=args.d3_url,
        bootstrap_css_url=args.bootstrap_css_url,
        bootstrap_js_url=args.bootstrap_js_url,
        row_order="clust" if as_bool(args.cluster_rows) else "alpha",
        col_order="clust" if as_bool(args.cluster_cols) else "alpha",
    )

    print(f"Saved Clustergrammer HTML to: {output_html}")
    print(f"Saved plotted data to: {output_table}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
