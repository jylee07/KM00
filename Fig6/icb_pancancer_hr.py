"""
Figure 6D: pan-cancer hazard ratio (HR) analysis of ICB (immune checkpoint
inhibitor) outcomes.

For the ICB-treated cohort, tests association between TTNT (time to next
treatment) and a set of clinical/genomic features (cancer type, treatment
line, TMB, MSI, MATH, pathway-level mutation status, gene-level mutation
status, mutational signature presence) using univariate Cox proportional
hazards models.

Input (public data release, identified by KM_ID / KM_SAMPLE_ID):
    - clinical_info (KM00_Supple1_Cohortinfo.xlsx): KM_ID, Abbreviation
    - ngs_meta_info (KM00_Supple2_NGSinfo.xlsx): KM_ID, KM_SAMPLE_ID, PANEL,
      Selected -- used to find each subject's selected tissue sample
    - icb (KM00_Supple7-1_ICB.xlsx): KM_ID, Line_of_Treatment, TTNT_months,
      TTNT_event
    - icb_genomic (KM00_Supple7-2_ICB_GenomicFeatures.xlsx): KM_ID,
      TMB_status, MSI_status, MATH_score, SBS1..SBS110 (proportions,
      already restricted to samples with Total_Mutations >= 20)
    - panel gene list (KM00_Panel_GeneList.csv): per-panel gene coverage,
      used to find the gene set common to the tissue panels actually used
      by the ICB cohort (CS, KM1.0, KM1.1 -- no ICB patient used FIRST)
    - maf (KM00_AllSamples_Combined_Mutations.maf): mutation calls, used to derive pathway- and
      gene-level mutation status (tissue panels only) -- a pathway/gene
      counts as mutated if it (or any of its genes) carries a mutation in
      the subject's tissue sample

Output:
    - icb_pancancer_hr.csv: one row per feature with HR, CI, p-value, FDR
"""
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

DATA_DIR = "../Data"
OUTPUT_PATH = "output/icb_pancancer_hr.csv"

CLINICAL_PATH = f"{DATA_DIR}/KM00_Supple1_Cohortinfo.xlsx"
NGS_META_PATH = f"{DATA_DIR}/KM00_Supple2_NGSinfo.xlsx"
ICB_PATH = f"{DATA_DIR}/KM00_Supple7-1_ICB.xlsx"
GENOMIC_PATH = f"{DATA_DIR}/KM00_Supple7-2_ICB_GenomicFeatures.xlsx"
MAF_PATH = f"{DATA_DIR}/KM00_AllSamples_Combined_Mutations.maf"
PANEL_GENELIST_PATH = f"{DATA_DIR}/KM00_Panel_GeneList.csv"

TISSUE_PANELS = ["CS", "FIRST", "KM1.0", "KM1.1"]
TIME_COL = "TTNT_months"
EVENT_COL = "TTNT_event"

# genes that define each pathway (restricted to the set well covered across
# all tissue panels)
PATHWAY_GENE_MAP = {
    "mmr": ["MSH2", "MSH6", "MLH1", "PMS2"],
    "hrd": [
        "ATM", "ATR", "BRCA1", "BRCA2", "PALB2", "RAD50", "RAD51", "RAD51B",
        "RAD51C", "RAD51D", "MRE11A", "BLM", "CHEK1", "CHEK2", "FANCA", "FANCC", "BARD1",
    ],
    "cell_cycle_TP53": ["CDKN2A", "MDM4", "CCND1", "CCNE1", "MDM2", "TP53", "RB1", "FBXW7"],
    "rtkras": [
        "ALK", "ARAF", "BRAF", "CBL", "EGFR", "ERBB2", "ERBB3", "ERBB4",
        "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FLT3", "HRAS", "IGF1R", "KIT",
        "KRAS", "MAP2K1", "MAP2K2", "MET", "NF1", "NRAS", "NTRK1", "NTRK2",
        "PDGFRA", "PTPN11", "RAC1", "RAF1", "RET", "ROS1",
    ],
    "pi3k": [
        "AKT1", "AKT2", "AKT3", "INPP4B", "MTOR", "PIK3CA", "PIK3R1", "PIK3R2",
        "PTEN", "RICTOR", "RPTOR", "STK11", "TSC1",
    ],
    "wnt": ["APC", "CTNNB1", "AXIN1", "LRP6", "RNF43"],
    "tgfbeta": ["SMAD2", "SMAD4", "TGFBR2"],
    "notch": ["CREBBP", "EP300", "KDM5A", "NCOR1", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "SPEN"],
}
PATHWAY_GENE_MAP["ddr"] = sorted(set(PATHWAY_GENE_MAP["hrd"]) | set(PATHWAY_GENE_MAP["mmr"]))


def load_icb_cohort():
    """ICB cohort with cancer type, treatment line, TTNT outcome, and genomic features."""
    clinical = pd.read_excel(CLINICAL_PATH)[["KM_ID", "Abbreviation"]]
    icb = pd.read_excel(ICB_PATH)
    genomic = pd.read_excel(GENOMIC_PATH).drop(columns="KM_SAMPLE_ID")

    icb = pd.merge(icb, clinical, on="KM_ID", how="left")
    icb = pd.merge(icb, genomic, on="KM_ID", how="left")
    return icb


def selected_tissue_sample_by_km_id():
    """KM_ID -> KM_SAMPLE_ID for each subject's selected tissue-panel sample."""
    ngs_meta = pd.read_excel(NGS_META_PATH)
    ngs_meta = ngs_meta[ngs_meta["PANEL"].isin(TISSUE_PANELS) & (ngs_meta["Selected"] == "O")]
    return ngs_meta.drop_duplicates(subset="KM_ID").set_index("KM_ID")["KM_SAMPLE_ID"]


def icb_common_panel_genes(icb_km_ids):
    """Gene set common to all tissue panels actually used by the ICB cohort."""
    ngs_meta = pd.read_excel(NGS_META_PATH)
    used_panels = sorted(ngs_meta[
        ngs_meta["PANEL"].isin(TISSUE_PANELS)
        & (ngs_meta["Selected"] == "O")
        & ngs_meta["KM_ID"].isin(icb_km_ids)
    ]["PANEL"].unique())

    panel = pd.read_csv(PANEL_GENELIST_PATH)
    gene_sets = [set(panel[col].dropna().drop_duplicates()) for col in used_panels]
    return sorted(set.intersection(*gene_sets))


def add_mutation_features(icb):
    """Pathway-<name> / Genetic-<gene> = 1 if mutated in the subject's tissue sample."""
    icb = icb.copy()
    icb["KM_SAMPLE_ID"] = icb["KM_ID"].map(selected_tissue_sample_by_km_id())
    has_sample = icb["KM_SAMPLE_ID"].notna()

    genes = icb_common_panel_genes(icb["KM_ID"])
    pathway_genes = {g for gs in PATHWAY_GENE_MAP.values() for g in gs}
    genes_to_load = set(genes) | pathway_genes

    maf = pd.read_csv(MAF_PATH, sep="\t")
    maf = maf[maf["KM_SAMPLE_ID"].isin(icb["KM_SAMPLE_ID"].dropna()) & maf["Hugo_Symbol"].isin(genes_to_load)]
    mutated_by_gene = maf.groupby("Hugo_Symbol")["KM_SAMPLE_ID"].apply(set).to_dict()

    for pathway, pw_genes in PATHWAY_GENE_MAP.items():
        mutated_samples = set().union(*(mutated_by_gene.get(g, set()) for g in pw_genes))
        icb[f"Pathway-{pathway}"] = np.where(
            icb["KM_SAMPLE_ID"].isin(mutated_samples), 1,
            np.where(has_sample, 0, np.nan),
        )

    for gene in genes:
        mutated_samples = mutated_by_gene.get(gene, set())
        icb[f"Genetic-{gene}"] = np.where(
            icb["KM_SAMPLE_ID"].isin(mutated_samples), 1,
            np.where(has_sample, 0, np.nan),
        )

    return icb.drop(columns="KM_SAMPLE_ID")


def add_signature_features(icb):
    """SBS<x>_present = 1 if that signature's proportion is > 0 (already QC-restricted)."""
    icb = icb.copy()
    sig_cols = [c for c in icb.columns if c.startswith("SBS")]
    for col in sig_cols:
        icb[f"{col}_present"] = np.where(icb[col].notna(), (icb[col] > 0).astype(float), np.nan)
    return icb.drop(columns=sig_cols)


def run_univariate_cox(df, feature, feature_type, min_n=20, min_event=5, min_positive=3):
    """Fit a univariate Cox model for one feature against TIME_COL/EVENT_COL."""
    tmp = df[[TIME_COL, EVENT_COL, feature]].dropna()

    result = {
        "feature": feature, "feature_type": feature_type, "n": len(tmp),
        "n_event": tmp[EVENT_COL].sum() if len(tmp) > 0 else np.nan,
        "n_positive": np.nan, "HR": np.nan, "CI_lower": np.nan, "CI_upper": np.nan,
        "p_value": np.nan, "status": "failed",
    }

    if len(tmp) < min_n:
        result["status"] = "too_few_samples"
        return result
    if tmp[EVENT_COL].sum() < min_event:
        result["status"] = "too_few_events"
        return result

    if feature_type == "binary":
        n_pos = (tmp[feature] == 1).sum()
        n_neg = (tmp[feature] == 0).sum()
        result["n_positive"] = n_pos
        if tmp[feature].nunique() < 2:
            result["status"] = "single_class"
            return result
        if n_pos < min_positive or n_neg < min_positive:
            result["status"] = "too_imbalanced"
            return result

    if feature_type == "continuous":
        if tmp[feature].nunique() < 3:
            result["status"] = "low_variance"
            return result
        std = tmp[feature].std()
        if std == 0 or pd.isna(std):
            result["status"] = "zero_variance"
            return result
        tmp[feature] = (tmp[feature] - tmp[feature].mean()) / std

    try:
        cph = CoxPHFitter()
        cph.fit(tmp, duration_col=TIME_COL, event_col=EVENT_COL)
        s = cph.summary.loc[feature]
        result.update({
            "HR": s["exp(coef)"], "CI_lower": s["exp(coef) lower 95%"],
            "CI_upper": s["exp(coef) upper 95%"], "p_value": s["p"], "status": "ok",
        })
    except Exception as e:
        result["status"] = f"error: {e}"

    return result


def run_categorical_cox(df, cat_col, min_n=20, min_event=5, min_positive=3):
    """Dummy-encode a categorical feature and fit a univariate Cox model per level."""
    tmp = df[[TIME_COL, EVENT_COL, cat_col]].dropna()
    if len(tmp) < min_n:
        return [{
            "feature": cat_col, "feature_type": "categorical", "n": len(tmp),
            "n_event": np.nan, "n_positive": np.nan, "HR": np.nan,
            "CI_lower": np.nan, "CI_upper": np.nan, "p_value": np.nan,
            "status": "too_few_samples", "original_feature": cat_col,
        }]

    dummies = pd.get_dummies(tmp[cat_col], prefix=cat_col)
    tmp2 = pd.concat([tmp[[TIME_COL, EVENT_COL]], dummies], axis=1)

    results = []
    for dummy_col in dummies.columns:
        res = run_univariate_cox(
            tmp2, dummy_col, "binary",
            min_n=min_n, min_event=min_event, min_positive=min_positive,
        )
        res["original_feature"] = cat_col
        results.append(res)
    return results


def main():
    icb = load_icb_cohort()
    icb = add_mutation_features(icb)
    icb = add_signature_features(icb)

    icb["TMB_High"] = (icb["TMB_status"] == "TMB-H").astype(float)
    icb["MSI_High"] = (icb["MSI_status"] == "MSI").astype(float)

    continuous_vars = ["MATH_score"]
    binary_vars = (
        ["TMB_High", "MSI_High"]
        + [c for c in icb.columns if c.startswith("Pathway-")]
        + [c for c in icb.columns if c.startswith("Genetic-")]
        + [c for c in icb.columns if c.endswith("_present")]
    )
    categorical_vars = ["Abbreviation", "Line_of_Treatment"]

    results = [run_univariate_cox(icb, f, "continuous") for f in continuous_vars]
    results += [run_univariate_cox(icb, f, "binary") for f in binary_vars]
    for f in categorical_vars:
        results.extend(run_categorical_cox(icb, f))

    result_df = pd.DataFrame(results)
    ok = result_df["status"] == "ok"
    result_df.loc[ok, "FDR"] = multipletests(result_df.loc[ok, "p_value"], method="fdr_bh")[1]
    result_df = result_df.sort_values("p_value", na_position="last").reset_index(drop=True)
    result_df.to_csv(OUTPUT_PATH, index=False)


if __name__ == "__main__":
    main()
