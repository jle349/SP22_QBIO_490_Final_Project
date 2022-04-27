import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

if __name__ == '__main__':
    prot_data = pd.read_csv("~/qbio490/qbio_data_analysis_echo/analysis_data/cptac/cptac_protein.csv", index_col=0)
    rna_data = pd.read_csv("~/qbio490/qbio_data_analysis_echo/analysis_data/cptac/cptac_rna.csv", index_col=0)
    clinical_data = pd.read_csv("~/qbio490/qbio_data_analysis_echo/analysis_data/cptac/cptac_clinical.csv",
                                             index_col=0)

    name_intersects = [np.intersect1d(rna_data.index, prot_data.index)
        , np.intersect1d(prot_data.index, clinical_data.index)
        , np.intersect1d(rna_data.index, clinical_data.index)
    ]

    all_data_patients = name_intersects[0]
    intersect_clin = clinical_data.loc[all_data_patients, :]
    intersect_rna = rna_data.loc[all_data_patients, :]
    intersect_prot = prot_data.loc[all_data_patients, :]

    is_fem = intersect_clin["Gender"] == "Female" # boolean masks to determine male and female patients
    is_male = intersect_clin["Gender"] == "Male"

    fem_intersect_rna = intersect_rna.loc[is_fem, :]
    male_intersect_rna = intersect_prot.loc[is_male, :]
    fem_intersect_prot = intersect_prot.loc[is_fem, :]
    male_intersect_rna = intersect_prot.loc[is_male, :]
    male_intersect_prot = intersect_prot.loc[is_male, :]

    ncomparisons = 5  # define variable in case we want to change the number of correlations

    # all genes come from top 10 overexpressed genes in female patients and are present in cptac database
    fem_gene_names = ['LAS1L', "SLC25A5", "CFH", "CDKN2A", "FUCA2"]

    # initialize a correlation dataframe to store correlation coefficients
    fcorr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                           index=fem_gene_names,
                           columns=fem_gene_names)

    # grab correlation coefficients for each gene in respect to every other gene
    for g1 in fem_gene_names:
        for g2 in fem_gene_names:
            # calculate the correlations between protein and RNA and store in dataframe
            corr, pval = stats.spearmanr(fem_intersect_rna[g1], fem_intersect_prot[g2], nan_policy="omit")
            fcorr_df.loc[g1, g2] = corr

    # all genes come from top 10 overexpressed genes in male patients and are present in cptac database
    male_gene_names = ['PCDH1', "THBS4", "CALD1", "CDKN2A", "FUCA2"]

    # initialize a correlation dataframe to store correlation coefficients
    mcorr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                           index=male_gene_names,
                           columns=male_gene_names)
    for g1 in male_gene_names:
        for g2 in male_gene_names:
            # calculate the correlations between protein and RNA and store in dataframe
            corr, pval = stats.spearmanr(male_intersect_rna[g1], male_intersect_prot[g2], nan_policy="omit")
            mcorr_df.loc[g1, g2] = corr

    # generate side-by-side heat maps to compare between sexes
    fig, ax = plt.subplots(1,2)

    # generate heatmaps for female patients
    sns.heatmap(
        fcorr_df,
        cmap='mako', ax = ax[0]
    ).set(title = "Female Patients (Fig. 5a)")

    # generate heat map for male patients
    sns.heatmap(
        mcorr_df,
        cmap='mako', ax = ax[1]
    ).set(title = "Male Patients (Fig. 5b)")

    plt.show()
    plt.savefig('rna_protein_sex.png', bbox_inches='tight')