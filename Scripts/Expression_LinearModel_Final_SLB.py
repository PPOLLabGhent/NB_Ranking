#!/usr/bin/env python3
# Linear model script utilizing expression and copy number variation data from the NRC Dataset (GSE85047)
# Author: Christophe Van Neste and adapted by Sarah-Lee Bekaert

import os
import zipfile
import pandas as pd
import numpy as np
import xarray as xr
from sklearn.preprocessing import normalize
from bidali.retro import ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, r # Added by SL
#from rpy2.rinterface import RRuntimeError - Adapted by SL

# Enable automatic conversion between pandas DataFrame and R data.frame
pandas2ri.activate()

# R functions
stats = importr('stats')

# Step 1: Load Expression and Metadata Files
expression_file = '/code/nb_ranking/InputData/mRNA_expression_data.txt'  # Update path
metadata_file = '/code/nb_ranking/InputData/20111216_NRC_samples.xlsx'   # Read file and skip first 4 lines because info 
cnvdata_file = '/code/nb_ranking/InputData/gene_level_cnvS_SBK_hg18_hg19_hg38.csv'    # Mapping based on lift over hg18 -> hg19 -> hg38


# Load the data (assuming tab-delimited)
expression = pd.read_csv(expression_file, sep='\t', index_col=0)
metadata = pd.read_excel(metadata_file, sheet_name='IDs NRC samples', index_col=0, skiprows=4)
cnv_data = pd.read_csv(cnvdata_file, index_col=0)
cnv_data = cnv_data.reset_index() # Sample is header now 

# Adding some extra lines because metadata is 319 and expression is 283 - 283 have all clinical info but only 255 CNV data -> in total only 219 match! 
# Make sure the NRC_ID column is stripped and in lowercase to ensure proper matching
metadata['NRC_ID'] = metadata['NRC_ID'].str.strip().str.lower()
# Get the column names from expression dataframe (strip and lowercase as well for matching)
expression_columns = expression.columns.str.strip().str.lower()
# Standardize sample names in CNV data
cnv_data['sample'] = cnv_data['sample'].str.strip().str.lower()

# Find the intersection of samples in all datasets
common_samples = set(metadata['NRC_ID']) & set(expression.columns) & set(cnv_data['sample'])
print(f"Number of common samples: {len(common_samples)}") #219 as exepcted 

# Filter metadata based on common samples
metadata = metadata[metadata['NRC_ID'].isin(common_samples)]
# Filter expression data (select only common columns)
expression = expression.loc[:, expression.columns.isin(common_samples)]
# Filter CNV data based on common samples
cnv_data = cnv_data[cnv_data['sample'].isin(common_samples)]

# Pivot the CNV data to have genes as rows and samples as columns 
cnv_pivot = cnv_data.pivot_table(index='gene', columns='sample', values='conugeal', aggfunc='first')
# Check the result - look good 
#print(cnv_pivot.head())

# Step 2: Rename the Probe ID's to Hugo Id's 
# Path to probe-gene mapping file
probe_to_gene_file = '/code/nb_ranking/InputData/mRNA info.txt'  # Update path
# Load the mapping file
probe_to_gene = pd.read_csv(probe_to_gene_file, sep='\t', header=None, names=['probe', 'gene'])
# Create a dictionary for mapping
probe_to_gene_dict = probe_to_gene.set_index('probe')['gene'].to_dict()

# Convert expression index and probe_to_gene keys to strings
expression.index = expression.index.astype(str)
probe_to_gene['probe'] = probe_to_gene['probe'].astype(str)
# Create the mapping dictionary
probe_to_gene_dict = probe_to_gene.set_index('probe')['gene'].to_dict()
# Map probe IDs to gene names
expression.index = expression.index.to_series().replace(probe_to_gene_dict)
# Remove unmapped probes (rows with NaN as index)
expression = expression[~expression.index.isna()]

# Extra Line to get common genes between expression and CNV Data 
# Standardize: remove spaces, convert to uppercase
expression.index = expression.index.str.strip().str.upper()
# If an index contains multiple genes, keep only the first gene
expression.index = expression.index.str.split(',').str[0]  

# Look for common genes -> 14170 gene
common_genes = expression.index.intersection(expression.index).intersection(cnv_pivot.index)

# Filtered files! 
expression_filtered = expression.loc[common_genes] #14899 in 219
expression_filtered = expression_filtered[~expression_filtered.index.duplicated(keep='first')] # remove duplicated genes -> 13529 in 219
cnv_pivot_filtered = cnv_pivot.loc[common_genes] #13529 in 219

# Normalize Filtered Expression Data
normalized_expression = pd.DataFrame(
    normalize(expression_filtered, axis=1),
    index=expression_filtered.index,
    columns=expression_filtered.columns
)

# Step 3: Create xarray Dataset
expression_filtered.index.name = 'gene' # Added this so Probset = gene 
data = xr.Dataset(
    {
        'expression': (['gene', 'sample'], expression_filtered.values),
        'normalized_expression': (['gene', 'sample'], normalized_expression.values),
        'conugeal': (['gene', 'sample'], cnv_pivot_filtered.values),  # Correct format now
    },
    coords={
        'gene': expression_filtered.index,  # Set 'gene' as a coordinate
        'sample': expression_filtered.columns  # Set 'sample' as a coordinate
    }
)

# Step 4: Process Metadata 
gene = 'BRIP1'                              # Example gene 
gene_data = data.loc[dict(gene=gene)]       # Example gene 

# Create Metadata in R format - based on INSS stage and MYCN amplification status 
# Create the concatenated labels in Python first
labels = (metadata.loc[:, 'nrc_inss'].apply(lambda x: 'highstage' if x in ('st3', 'st4') else 'lowstage') + '_' +
          metadata.loc[:, 'nrc_mycna'].apply(lambda x: 'amp' if x == 'yes' else 'sc')).tolist()

# Convert to an R character vector before creating a factor
labels_r = ro.StrVector(labels)

# Create the factor in R and relevel
metadata_r = ro.r.relevel(
    ro.r.factor(labels_r),  # Convert the list of strings to a factor
    ref='lowstage_sc'       # Set the reference level
)
#print(metadata_r)
metadata.set_index('NRC_ID', inplace=True) # NRC_ID needs to be index sample names

# Step 5: Linear Model Function
def calcLM(gene_data):
    # If a gene has a NaN value in any patient, exclude that sample from both expression and CNV data
    valid_samples = gene_data.dropna(dim="sample", subset=["normalized_expression", "conugeal"])
    # If no valid samples remain after filtering, return NaN
    if valid_samples.sizes["sample"] == 0:
        print(f"No valid samples available for gene {gene_data.gene}.")
    # Extract filtered values
    gexpr = valid_samples.normalized_expression.values
    cnfocus = valid_samples.conugeal.values
    filtered_sample_names = valid_samples.sample.values  # Get the remaining sample names
    # Ensure metadata is filtered to match valid_samples
    filtered_metadata = metadata.loc[filtered_sample_names]  # Keep only samples that remain
    # Recreate metadata labels **after filtering**
    labels = (filtered_metadata['nrc_inss'].apply(lambda x: 'highstage' if x in ('st3', 'st4') else 'lowstage') + '_' +
              filtered_metadata['nrc_mycna'].apply(lambda x: 'amp' if x == 'yes' else 'sc')).tolist()
    # Convert metadata to R format
    labels_r = ro.StrVector(labels)
    metadata_r = ro.r.relevel(ro.r.factor(labels_r), ref='lowstage_sc')  # Ensure correct reference level
    try:
        # Transform cnfocus into R series
        cnfocus_r = ro.r['relevel'](ro.r['factor'](cnfocus), ref='normal')
        gexpr_r = ro.FloatVector(gexpr)    
        # Check which categories are present in cnfocus_r
        cnfocus_levels = ro.r.levels(cnfocus_r)
        cnfocus_levels = list(cnfocus_levels)     
        # Create a DataFrame based on the presence of the levels (normal, gain, loss)
        if 'gain' in cnfocus_levels and 'loss' in cnfocus_levels:
            # Both gain and loss are present
            cnfocus_r = ro.r.relevel(cnfocus_r, ref='normal')  # Make sure 'normal' is the reference
            fmla = ro.Formula('gexpr ~ cnfocus + metadata')
            fmla.environment['gexpr'] = gexpr_r
            fmla.environment['cnfocus'] = cnfocus_r
            fmla.environment['metadata'] = metadata_r 
            fit_r = stats.lm(fmla)
            coefs = ro.r.summary(fit_r).rx2('coefficients')         
            # Convert coefficients to DataFrame with unique name for gain_loss case
            gain_loss_df = pd.DataFrame(
                np.array(coefs),
                columns=["Estimate", "Std. Error", "t value", "Pr(>|t|)"],
                index=["(Intercept)", "cnfocusgain", "cnfocusloss", "metadatahighstage_amp", "metadatahighstage_sc", "metadatalowstage_amp"]
            )         
            # Return DataFrame for gain and loss
            return xr.DataArray(gain_loss_df)    
        elif 'gain' in cnfocus_levels:
            # Only gain and normal present
            cnfocus_r = ro.r.relevel(cnfocus_r, ref='normal')
            fmla = ro.Formula('gexpr ~ cnfocus + metadata')
            fmla.environment['gexpr'] = gexpr_r
            fmla.environment['cnfocus'] = cnfocus_r
            fmla.environment['metadata'] = metadata_r 
            fit_r = stats.lm(fmla)
            coefs = ro.r.summary(fit_r).rx2('coefficients')        
            # Convert coefficients to DataFrame with unique name for gain case
            gain_df = pd.DataFrame(
                np.array(coefs),
                columns=["Estimate", "Std. Error", "t value", "Pr(>|t|)"],
                index=["(Intercept)", "cnfocusgain", "metadatahighstage_amp", "metadatahighstage_sc", "metadatalowstage_amp"]
            ) 
            # Return DataFrame for gain
            return xr.DataArray(gain_df) 
        elif 'loss' in cnfocus_levels:
            # Only loss and normal present
            cnfocus_r = ro.r.relevel(cnfocus_r, ref='normal')
            fmla = ro.Formula('gexpr ~ cnfocus + metadata')
            fmla.environment['gexpr'] = gexpr_r
            fmla.environment['cnfocus'] = cnfocus_r
            fmla.environment['metadata'] = metadata_r 
            fit_r = stats.lm(fmla)
            coefs = ro.r.summary(fit_r).rx2('coefficients')    
            # Convert coefficients to DataFrame with unique name for loss case
            loss_df = pd.DataFrame(
                np.array(coefs),
                columns=["Estimate", "Std. Error", "t value", "Pr(>|t|)"],
                index=["(Intercept)", "cnfocusloss", "metadatahighstage_amp", "metadatahighstage_sc", "metadatalowstage_amp"]
            )     
            # Return DataFrame for loss
            return xr.DataArray(loss_df)
    except Exception as e:
        print(f"Error processing gene data for {gene_data.gene}: {e}")
        return xr.DataArray(np.nan, dims=['gene', 'sample'], coords={'gene': [gene_data.gene], 'sample': gene_data.sample_names})

# There are some NAN values in the copy number data, count them for each gene 
nan_output_file = "/code/nb_ranking/OutputData/nan_counts_per_gene.csv"

# Open the file in write mode to start fresh
with open(nan_output_file, "w") as f:
    f.write("gene,NaN_count_conugeal\n")  # Write header

for gene_name in data.gene.values:
    try:
        gene_data = data.sel(gene=gene_name)
        nan_count = gene_data["conugeal"].isnull().sum().item()  # Use .isnull() instead of .isna()
        # Append results to file
        with open(nan_output_file, "a") as f:
            f.write(f"{gene_name},{nan_count}\n") 
    except Exception as e:
        print(f"Error processing gene {gene_name}: {e}")

# Now remove the genes that have more then 44 NAN's out of the 219 patients 
nan_counts_df = pd.read_csv('/code/nb_ranking/OutputData/nan_counts_per_gene.csv')

# Convert the 'gene' column to a list for easy lookup
nan_counts_dict = dict(zip(nan_counts_df['gene'], nan_counts_df['NaN_count_conugeal']))

# Loop over all the genes and apply calcLM function
all_results = []
for gene_name in data.gene.values:
    try:
        # Skip gene if its NaN count is greater than 44
        if nan_counts_dict.get(gene_name, 0) > 44:
            print(f"Skipping gene {gene_name} due to NaN count > 44.")
            continue
        # Select data for the current gene using .sel
        gene_data = data.sel(gene=gene_name) 
        # Apply the calcLM function to this gene data
        result = calcLM(gene_data)
        if result is not None:
            # Assign coordinates for the gene
            result = result.assign_coords({'gene': gene_name})
            all_results.append(result)  # Append result to the list    
        # print(f"Processed gene: {gene_name}")
    except Exception as e:
        print(f"Error processing gene {gene_name}: {e}")

# Concatenate all results into a single xarray.Dataset
final_results = xr.concat(all_results, dim='gene')

# Print the concatenated result to check
#print(final_results)

# Flatten the xarray dataset to a pandas DataFrame
final_results.name = "linear_model_results"  # Add name for results 
results_df = final_results.to_dataframe().reset_index()

# Filter the DataFrame to include only 'Estimate' values
results_estimate = results_df[results_df['dim_1'] == 'Estimate']

# Pivot the dataframe: `gene` becomes rows, `dim_0` (term) becomes columns
results_pivot = results_estimate.pivot_table(index='gene', columns='dim_0', values='linear_model_results')

# Flatten the multi-index columns and reset the index
results_pivot.columns = [col.strip() for col in results_pivot.columns.values]
results_pivot.reset_index(inplace=True)

# Save the formatted DataFrame to a CSV file
results_pivot.to_csv('/code/nb_ranking/OutputData/gene_resuls_lm_SLB_V4_hg38.csv', index=False)
