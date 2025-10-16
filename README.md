**CRISPRa bioinformatics workflow**

Bioinformatics pipeline based on limma-voom strategy for CRISPR activation screening data.

**Install**

To install the last version use the following command in your terminal:
<pre>mkdir CRISPRa_voom 
git clone https://github.com/seven1112233/CRISPRa_bioinformatics_workflow.git </pre>

**Examples**

Run `Rscript CRISPRa_workflow_T0med3mad_normbyNTC.R -h` in your terminal to get details.

A simple example:
<pre>Rscript CRISPRa_workflow_T0med3mad_normbyNTC.R\n

                           -i SP1_K562_mageck.count.txt\n

                           -o SP1_K562_voom\n

                           -p 25101010_SP1_K562\n

                           -l SP1_K562\n

                           -c T17-T0,T10-T0 </pre>


**Important notes:**

1. The input should be raw counts produced by MAGeCK count with columns separated by tab;
2. Two replicates for each time point are required;
3. The column names for input should be: sgRNA, Gene, SPn_Cellline_Tn_repn. For example: SP1_K562_T0_rep1, or SP1_H23_T24_DAR_rep1;
4. Each paramater should be specified.

Please let Xiaozhen know if you have any questions!
