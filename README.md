# :cat: CRISPRa bioinformatics workflow

Bioinformatics pipeline based on limma-voom strategy for CRISPR activation screening data.

**:paw_prints: Install**

The pipeline was developed based on R version `R-4.4.3`. 
The required R packages are listed in `env/packages_requirements.txt`.\
To install the last version use the following command in your terminal:
```bash
mkdir CRISPRa_voom 
cd CRISPRa_voom
git clone https://github.com/seven1112233/CRISPRa_bioinformatics_workflow.git
cd CRISPRa_bioinformatics_workflow</pre>
```
**:star: Examples**

Run `Rscript bin/CRISPRa_workflow_T0med3mad_normbyNTC.R -h` in your terminal to get details.

A simple example :sun_with_face::
```bash
Rscript bin/CRISPRa_workflow_T0med3mad_normbyNTC.R \
                           -i data/250826_SP2_HCT116_mageck.count.txt \
                           -o 250826_SP2_HCT116_voom \
                           -p 250826_SP2_HCT116 \
                           -l SP2_HCT116 \
                           -c T17-T0
```
Output files including: \
Figures:
1. 250826_SP2_HCT116_Ref.aveLogCPM_dist.pdf\
   The distribution of sgRNA raw counts in T0 samples
2. 250826_SP2_HCT116_sgRNA_differential_abundance_limma_mean-variance.pdf\
   The voom mean variance plot and sample-specific weights
3. 250826_SP2_HCT116_log_norm_counts_dist.pdf\
   The distribution of sgRNA raw counts and normalized counts for each sample.
4. 250826_SP2_HCT116_sgRNA_differential_abundance_limma_log2normedC_corr.pdf\
   The correlation of normalized counts between replicates
5. 250826_SP2_HCT116_sgRNA_differential_abundance_limma_volcano.pdf\
   The volcano plot for the differential abundance analysis at the sgRNA level.\
   Threshold: |FC| >= 2 & pvalue < 0.05
6. 250826_SP2_HCT116_promoter_differential_abundance_limma_volcano.pdf\
   The volcano plot for the differential abundance analysis at the promoter level based on all of guides.\
   Threshold: |FC| >= 1.5 & pvalue < 0.05

Tables:
1. 250826_SP2_HCT116_T17_vs_T0_sgRNA_differential_abundance_limma.txt\
   The differential abundance analysis at the sgRNA level.\
   The key columns are sgRNA, T17_vs_T0(LFC), P.value, adj.P.Val(adjust pvalue).
2. 250826_SP2_HCT116_T17_vs_T0_promoter_level_fry_results.txt\
   The differential abundance analysis at the promoter level.\
   The key columns are promoter, Gene, LFC_median(the median of LFC of all guides), NGenes(the number of guides), PValue, FDR.

**:exclamation::exclamation::exclamation: Important notes:**

1. The input should be raw counts produced by MAGeCK count with columns separated by tab;
2. Two replicates for each time point are required;
3. The column names for input should be: sgRNA, Gene, SPn_Cellline_Tn_repn. For example: SP1_K562_T0_rep1, or SP1_H23_T24_DAR_rep1;
4. Each paramater should be specified.

:question: Please let Xiaozhen know if you have any questions!
