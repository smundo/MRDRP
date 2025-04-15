import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import psutil
import argparse
import concurrent.futures
warnings.filterwarnings('ignore')


def Rsquared(minor_allelefreq, beta):
    return 2*minor_allelefreq*(1-minor_allelefreq)*beta**2

def Fstat(rsquared, n):
    return rsquared*(n - 2)/(1 - rsquared)


def get_data(GWAS_dir, sep):

    """
    Collect full GWAS summary statistics for exposures of interest.

    Args:
        GWAS_dir (str): Directory containing GWAS subdirectories and files for
                        exposures of interest.
    
        sep (str): Delimiter used in individual GWAS files.

    Returns:
        List of dataframes for each exposure with full GWAS statistics.
    """

    
    exposure_list = sorted(glob.glob(os.path.join(GWAS_dir, '*')))
    
    exposure_chrom = []

    for i in exposure_list:
        files_list = glob.glob(os.path.join(i, '*'))
        exposure_chrom.append(files_list)

    full_gwas = []
    for i in range(len(exposure_chrom)):
    #for i in range(10):
        full_genome_gwas = pd.concat([pd.read_table(f, delimiter = sep) for f in exposure_chrom[i]], ignore_index = True)
        full_gwas.append(full_genome_gwas)

    return full_gwas


def map_rsID(gwas_list):

    """
    Use Olink mapping files to map rsIDs to SNPs

    Args:
        gwas_list: List of dataframes with summary statistics for each exposure.

    Returns:
        Updated list of dataframes with summary statistics for each exposure,
        each now containing a column for the rsIDs.
    """
    
    olink_directory = '/data/feng_lab/sergio/Olink_rsID_maps'
    olink_file_pattern = '*.tsv'
    olink_map_files = glob.glob(olink_directory + '/' + olink_file_pattern)
    
    olink_rsid_fullgenome_map = pd.concat([pd.read_table(f, delimiter = '\t') for f in olink_map_files], ignore_index = True)

    gwas_full_update = []
    for i in range(len(gwas_list)):
        gwas = gwas_list[i]
        gwas['LOG10P'] = 10**(np.negative(gwas['LOG10P']))
        gwas.rename(columns={'LOG10P': 'p-value'}, inplace = True)
        merged_withrsid = pd.merge(olink_rsid_fullgenome_map, gwas, how='inner', on='ID')
        gwas_full_update.append(merged_withrsid)

    return gwas_full_update


def select_sig_variants(GWAS_dir, output_dir, gwas_list, pval, POS, sig_threshold=5e-8, variant_type='cis', CHROM=None, genomic_coordinates=None, window=None):
    """
    Select significant variants for exposures of interest.

    Args:
        GWAS_dir (str): Directory containing GWAS subdirectories and files for
                        exposures of interest.

        output_dir (str): Directory to which to send output files.
    
        gwas_list: List containing GWAS statistics for exposures
    
        pval (str): Label used in GWAS statistics files to denote the p-value column

        POS (str): Label used in GWAS statistics files to denote the genomic position column.
    
        sig_threshold (float): p-value significance threshold
    
        variant_type (str): Type of the variants to analyze. Options are 'cis'
                            for cis-acting or 'cis+trans' for searching across genome.
    
        CHROM (str): Label used in GWAS statistics files to denote the chromosome column.
    
        genomic_coordinates : List of tuples with chromosome number and gene start and
                              stop coordinates for cis-variant analysis
    
        window: Amount (in bp) beyond genomic coordinates within which to search for
                cis-variants (e.g. 300000 will look for variants within 300000 bp of gene)

    Output:
        csv files for each exposure containing summary statistics for significant variants.
    """

    exposure_list = sorted(glob.glob(os.path.join(GWAS_dir, '*')))

    for i in range(len(gwas_list)):
        if variant_type == 'cis':
            try:
                gwas = gwas_list[i]
                chromosome = genomic_coordinates[i][0]
                lower_limit_pos = genomic_coordinates[i][1] - window
                upper_limit_pos = genomic_coordinates[i][2] + window
                gwas_cis_bychrom = gwas[gwas[CHROM] == chromosome]
                gwas_cis = gwas_cis_bychrom[(gwas_cis_bychrom[POS] > lower_limit_pos) & (gwas_cis_bychrom[POS] < upper_limit_pos)]
                gwas_cis = gwas_cis.sort_values(POS)
                gwas_genomesignificant = gwas_cis[gwas_cis[pval] < sig_threshold]

                # gwas_genomesignificant.loc[gwas_genomesignificant.A1FREQ > 0.5, 'MAF'] = 1 - gwas_genomesignificant['A1FREQ']
                # gwas_genomesignificant.loc[gwas_genomesignificant.A1FREQ < 0.5, 'MAF'] = gwas_genomesignificant['A1FREQ']

                # gwas_genomesignificant['R^2'] = gwas_genomesignificant.apply(lambda x: Rsquared(x['MAF'], x['BETA']), axis=1)
                # gwas_genomesignificant['F-statistic'] = gwas_genomesignificant.apply(lambda x: Fstat(x['R^2'], x['N']), axis=1)

                # gwas_genomesignificant_final = gwas_genomesignificant[gwas_genomesignificant['F-statistic'] > 10]
                
                dir_name = os.path.basename(exposure_list[i])
                if not gwas_genomesignificant.empty:
                    gwas_genomesignificant.to_csv(f'{output_dir}/{dir_name}_full_sig_variants.csv')
                elif gwas_genomesignificant.empty:
                    print(f"Dataframe for {dir_name} is empty, not writing to file.")
            except Exception as e:
                print(f"Error processing {dir_name}: {e}")
                continue
        elif variant_type == 'cis+trans':
            gwas = gwas_list[i]
            gwas_genomesignificant = gwas[gwas[pval] < sig_threshold]
            dir_name = os.path.basename(exposure_list[i])
            if not gwas_genomesignificant.empty:
                gwas_genomesignificant.to_csv(f'{output_dir}/{dir_name}_full_sig_variants.csv')
            elif gwas_genomesignificant.empty:
                print(f"Dataframe for {dir_name} is empty, not writing to file.")
                
