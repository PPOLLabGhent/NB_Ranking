## Neuroblastoma Gene Ranking Analysis 

**Preparing Input Files**

In this study, we used the Neuroblastoma Research Consortium (NRC) dataset (GSE85047). 

*Copy Number Variation (CNV)* data in Circular Binary Segmentation (CBS) format:
  1.	Lift over from hg18 -> hg19 -> hg38 using chain files from the UCSC Genome Browser
  2.	Annotation of the regions using biomaRt 
  3.	Classify CNV regions based on Depuydt et al. (2018) thresholds
     
      - Gain: Mean CNV ratio ≥ 0.15
      - Loss: Mean CNV ratio ≤ -0.25
      - Normal: CNV ratios falling between these thresholds

The final results, including the CNV status (gain, loss, or normal) for each gene and sample, were compiled into a CSV file: *gene_level_cnvS_SBK_hg18_hg19_hg38.csv*. This first steps were done using *Transform_CNVData_hg18_to_hg38.R*. 

*Gene Expression* data mapped to the hg18 reference genome:
  1.	mRNA_expression_data.txt containing gene expression data and probe ID’s
     
     - Probe-level identifiers were mapped to gene symbols using mRNA info.txt
     - Probes that did not map to a known gene were excluded. 
     - If multiple probes mapped to the same gene, the first occurrence was retained
     - Expression data were normalized using L2 normalization (sklearn.preprocessing.normalize) across genes to standardize expression values.

  2.	20111216_NRC_samples.xlsx containing metadata information
     
     - Annotation was based on the INSS stage and the MYCN amplification status 
     - Stage 3 and 4 with a MYCN amplification – highstage_ampl
     - Stage 3 and 4 without a MYCN amplification – highstage_sc 
     - Stage 1, 2 and 4s with a MYCN amplification – lowstage_ampl
     - Stage 1, 2 and 4s without a MYCN amplification – lowstage_sc 

The intersection of sample IDs across metadata, expression, and CNV datasets identified 219 common samples for analysis.

**Linear Model**

For each gene, a linear model was fitted using the R stats package within a Python environment via rpy2. The model included normalized gene expression as the dependent variable and CNV status (normal, gain, or loss) alongside metadata categories as independent variables. CNV status was treated as a categorical factor, with normal set as the reference level. Depending on the presence of gain and/or loss categories, models were adapted accordingly. Genes with more than 44 missing values in CNV data (out of 219 patients) were excluded from analysis. The remaining missing values were handled by filtering out samples with NaN values for either expression or CNV status before model fitting.

The output of the linear model was saved as a CVS file: *gene_results_lm_SLB_V4_hg38.csv* containing the Estimate values for each metadata group. This was done using *Expression_LinearModel_Final_SLB.py*. 

**Ranking genes**

Based on the linear model, a ranking of genes per chromosomal arm is performed using the *NB_Ranking_SLB.py* script. 

The analyzeCNAspread function will determine how frequently a chromosomal arm is altered across patients.  This is based on a dataset combined of 556 high risk neuroblastoma patients from Depuydt et al. (2018). 

  1.	Segmentation analysis of the PDP, where thresholds classify regions as gains or losses based on their log2 ratios.
     
  2.	Each altered region is mapped to its corresponding genes using genomic coordinates (hg38)
     
  3.	The penetrance (frequency) of each CNA is calculated by dividing the number of patients exhibiting the alteration by the total patient count.

The next function rankNBaccordingToMYCNstatus will use this penetrance value and the linear model from the previous step to rank genes per chromosomal arm and according to their MYCN status. 

  1.	Calculate the frequency of chromosomal alterations (both gains and losses) within each MYCN group
     
  2.	Only chromosomal regions altered in at least 25% of cases are considered
     
  3.	Evaluate key properties for each altered chromosomal arm:
     
     - Gene Density: The number of genes per unit length of the chromosome.
     - Gain/Loss Ratio: The proportion of gains versus losses on a chromosomal arm.
     - Significance Testing: A binomial test checks whether gains or losses are more dominant, ensuring only statistically significant alterations are kept.

  4.	Aggregate CNAs for each patient and determines the frequency with which each gene is altered across all samples, producing a gene penetrance score
     
  6.	Genes within significantly altered chromosomal arms are ranked based on three biological metrics:
     
     - Copy Number Frequency (CNrank): Genes altered frequently across samples are ranked higher
     - Dosage Sensitivity (dosagerank): Derived from linear models predicting changes in gene expression due to CNAs, this metric assesses how sensitive gene expression is to copy number variations (higher sensitivity = better rank)
     - Survival Association/Risk Classification (riskrank): Based on statistical models linking gene alterations to patient survival outcomes

  7.	For each chromosomal arm, a final combined ranking is computed based on three key parameters:
     
     - Penetrance: This score reflects the fraction of patients in which that gene is altered.
     - Dosagelm: The correlation between copy number variation and gene expression.
     - Risklm: The association between gene alterations and survival risk classification.

To generate a comprehensive ranking score, the three parameters are first normalized using percentile ranking (percentile rank = rank / total count). Then, these ranks are multiplied across the three parameters to produce a final combined ranking score for each gene.

**Get started via Docker**

