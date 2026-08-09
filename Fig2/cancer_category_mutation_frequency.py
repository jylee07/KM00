"""
Figure 2B: mutation frequency comparison across cancer functional categories.

Restricted to tissue-panel samples only (AXEN excluded), one selected sample
per subject. For each functional category (e.g. Thoracic, Upper GI, Intestinal
GI, ...), tests whether each gene's mutation frequency in that category
differs from the rest of the cohort (Fisher's exact test, BH-FDR corrected
within each category).

Input (public data release, identified by KM_ID / KM_SAMPLE_ID):
    - clinical_info (KM00_Supple1_Cohortinfo.xlsx): one row per subject
      (KM_ID, Abbreviation, ...); cancer functional category is derived here
      from Abbreviation
    - ngs_meta_info (KM00_Supple2_NGSinfo.xlsx): one row per NGS sample
      (KM_ID, KM_SAMPLE_ID, PANEL, Selected); Selected == 'O' marks which
      sample was chosen for analysis when a subject has more than one
    - gene panel list (per-panel gene coverage, used to get the gene set
      common to all tissue panels)
    - maf (KM00_AllSamples_Combined_Mutations.maf): combined all-sample mutation calls

Output:
    - fisher_cancer_category.csv: one row per (Cancer_Category, Gene) with
      mutation frequency in-category vs rest, p-value, FDR, and frequency diff.
"""
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

DATA_DIR = "data"
OUTPUT_PATH = "output/fisher_cancer_category.csv"

CLINICAL_PATH = f"{DATA_DIR}/KM00_Supple1_Cohortinfo.xlsx"
NGS_META_PATH = f"{DATA_DIR}/KM00_Supple2_NGSinfo.xlsx"
PANEL_GENELIST_PATH = f"{DATA_DIR}/KM00_Panel_GeneList.csv"
MAF_PATH = f"{DATA_DIR}/KM00_AllSamples_Combined_Mutations.maf"

# NGS panels considered "tissue" panels (AXEN is a separate, non-tissue panel)
TISSUE_PANELS = ["CS", "FIRST", "KM1.0", "KM1.1"]

# maps each cancer type abbreviation to its functional category
CANCER_CATEGORY_MAP = {
    "Thoracic": ["THYM", "MESO", "LUCA"],
    "UpperGI": ["ESCA", "GC"],
    "IntGI": ["CRC", "SBC", "APCA", "ANCA"],
    "DevGI": ["HCC", "CHOL", "GBC", "AMPCA", "PACA"],
    "Genitourinary": ["KIRC", "ADCA", "URTC", "PRCA", "TECA", "PSCC"],
    "Gynecologic": ["OVCA", "UTER", "CC", "VC"],
    "Skin": ["MEL", "SKCA"],
}


def assign_cancer_category(clinical):
    clinical = clinical.copy()
    for category, abbreviations in CANCER_CATEGORY_MAP.items():
        clinical.loc[clinical["Abbreviation"].isin(abbreviations), "Cancer_Category"] = category
    return clinical


def load_tissue_sample_info():
    """Selected, tissue-panel NGS samples merged with each subject's cancer functional category."""
    clinical = assign_cancer_category(pd.read_excel(CLINICAL_PATH))

    ngs_meta = pd.read_excel(NGS_META_PATH)
    ngs_meta = ngs_meta[ngs_meta["PANEL"].isin(TISSUE_PANELS) & (ngs_meta["Selected"] == "O")]

    return pd.merge(
        ngs_meta,
        clinical[["KM_ID", "Abbreviation", "Cancer_Category"]],
        on="KM_ID",
        how="inner",
    )


def common_tissue_panel_genes():
    """Gene set covered by all tissue-panel versions (CS, FIRST, KM1.0, KM1.1)."""
    panel = pd.read_csv(PANEL_GENELIST_PATH)
    gene_sets = [panel[col].dropna().drop_duplicates().tolist() for col in TISSUE_PANELS]
    return sorted(set.intersection(*map(set, gene_sets)))


def fisher_test_by_category(sample_info, maf, genes):
    """Fisher's exact test of mutation frequency: each category vs the rest of the cohort."""
    sample_category = (
        sample_info.drop_duplicates(subset="KM_SAMPLE_ID").set_index("KM_SAMPLE_ID")["Cancer_Category"]
    )
    all_samples = set(sample_info["KM_SAMPLE_ID"])
    categories = sorted(sample_category.dropna().unique())

    mutated_samples_by_gene = (
        maf[maf["Hugo_Symbol"].isin(genes) & maf["KM_SAMPLE_ID"].isin(all_samples)]
        [["KM_SAMPLE_ID", "Hugo_Symbol"]]
        .drop_duplicates()
        .groupby("Hugo_Symbol")["KM_SAMPLE_ID"]
        .apply(set)
        .to_dict()
    )

    results = []
    for category in categories:
        cat_samples = set(sample_category[sample_category == category].index)
        rest_samples = all_samples - cat_samples  # includes samples with no assigned category
        n_cat, n_rest = len(cat_samples), len(rest_samples)

        for gene in genes:
            gene_mut = mutated_samples_by_gene.get(gene, set())
            mut_cat = len(gene_mut & cat_samples)
            mut_rest = len(gene_mut & rest_samples)
            wt_cat, wt_rest = n_cat - mut_cat, n_rest - mut_rest

            _, pvalue = fisher_exact([[mut_cat, wt_cat], [mut_rest, wt_rest]], alternative="two-sided")
            results.append(
                {
                    "Cancer_Category": category,
                    "Gene": gene,
                    "mut_cat": mut_cat,
                    "wt_cat": wt_cat,
                    "mut_rest": mut_rest,
                    "wt_rest": wt_rest,
                    "pvalue": pvalue,
                }
            )

    fisher_df = pd.DataFrame(results)
    fisher_df["fdr"] = fisher_df.groupby("Cancer_Category")["pvalue"].transform(
        lambda p: multipletests(p, method="fdr_bh")[1]
    )
    fisher_df["mut_cat_freq"] = fisher_df["mut_cat"] / (fisher_df["mut_cat"] + fisher_df["wt_cat"])
    fisher_df["mut_rest_freq"] = fisher_df["mut_rest"] / (fisher_df["mut_rest"] + fisher_df["wt_rest"])
    fisher_df["freq_diff"] = fisher_df["mut_cat_freq"] - fisher_df["mut_rest_freq"]
    return fisher_df.sort_values(["Cancer_Category", "pvalue"]).reset_index(drop=True)


def main():
    sample_info = load_tissue_sample_info()
    genes = common_tissue_panel_genes()
    maf = pd.read_csv(MAF_PATH, sep="\t")

    fisher_df = fisher_test_by_category(sample_info, maf, genes)
    fisher_df.to_csv(OUTPUT_PATH, index=False)


if __name__ == "__main__":
    main()
