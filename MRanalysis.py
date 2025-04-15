import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import psutil
import re
warnings.filterwarnings('ignore')

import seaborn as sns
import rpy2
import rpy2.robjects as robjects

## To aid in printing HTML in notebooks
import rpy2.ipython.html
rpy2.ipython.html.init_printing()

## To see plots in an output cell
from rpy2.ipython.ggplot import image_png

from rpy2.robjects.packages import importr, data
from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib.embedded import RRuntimeError
pandas2ri.activate()

primary = "#404040"
accent = "#24B8A0"
pal = sns.color_palette([primary, accent])
style = {
    "axes.edgecolor": primary,
    "axes.labelcolor": primary,
    "text.color": primary,
    "xtick.color": primary,
    "ytick.color": primary,
}

sns.set_palette(pal)
sns.set_context("talk")
sns.set_style("ticks", style)

utils = importr('utils')
base = importr('base')
utils.chooseCRANmirror(ind=1)

utils.install_packages('TwoSampleMR')
TwoSampleMR = importr('TwoSampleMR')

import rpy2.robjects.lib.ggplot2 as ggplot2

plink_binary = robjects.r("genetics.binaRies::get_plink_binary()")

utils.install_packages('ieugwasr')
ieugwasr = importr('ieugwasr')

grdevices = importr('grDevices')

rprint = robjects.globalenv.find("print")

def LD_clump(sig_var_dir, rsid, beta, se, effect_allele, other_allele, eaf, pval, ref_genome_file=None, clump_threshold=0.01, population='EUR', samplesize=None, chrom=None, pos=None, delimiter=',', file_pattern='*.csv', local_clump=False):
    """
    Perform LD clumping for each exposure.

    Args:
        sig_var_dir (str): Directory containing csv files of significant variants
                           for each exposure (these files are the output from get_data())

        rsid (str): Label used for rsID column in exposure dataframe/csv file
    
        beta (str): Label used for effect size column
    
        se (str): Label used for standard error of effect size column

        effect_allele (str): Label used for effect allele column
    
        other_allele (str): Label used for other allele column
    
        eaf (str): Label used for effect allele frequency column
    
        pval (str): Label used for p-value column
    
        ref_genome_file (str): File for reference genome for local clumping
    
        clump_threshold (float): R^2 threshold for clumping

        population (str): Genetic ancestry of interest from reference genome

        samplesize (str): Label used for sample size column

        chrom (str): Label used for chromosome column

        pos (str): Label used for position column

        delimiter (str): Delimiter used in exposure dataframes

        file_pattern (str): Exposure dataframe file pattern

        local_clump (boolean): Whether local or remote clumping is desired

    Output:
        List of exposure IVs for each exposure
    """
    
    sig_var_files = sorted(glob.glob(sig_var_dir + '/' + file_pattern))
    #exp_names = [os.path.splitext(os.path.basename(sig_var_files[i]))[0] for i in range(len(sig_var_files))]

    exposure_IVs = []
    error_indices = []
    for i in range(len(sig_var_files)):
        
        if local_clump:
            try:
                exp_data = TwoSampleMR.read_exposure_data(filename = sig_var_files[i], clump = False, sep = delimiter, snp_col = rsid, beta_col = beta, se_col = se, effect_allele_col = effect_allele, other_allele_col = other_allele, eaf_col = eaf, pval_col = pval, samplesize_col = samplesize, chr_col = chrom, pos_col = pos)
                
                with (robjects.default_converter + pandas2ri.converter).context():
                    exp_data = robjects.conversion.get_conversion().rpy2py(exp_data)

                exp_data = exp_data.rename(columns={"SNP": "rsid", "pval.exposure":"pval", "id.exposure":"id"})

                with (robjects.default_converter + pandas2ri.converter).context():
                    exp_data = robjects.conversion.get_conversion().py2rpy(exp_data)

                IVs = ieugwasr.ld_clump(dat = exp_data, clump_r2 = clump_threshold, pop = population, plink_bin = plink_binary, bfile = ref_genome_file)

                with (robjects.default_converter + pandas2ri.converter).context():
                    IVs = robjects.conversion.get_conversion().rpy2py(IVs)

                IVs = IVs.rename(columns={"rsid": "SNP", "pval":"pval.exposure", "id":"id.exposure"})

                with (robjects.default_converter + pandas2ri.converter).context():
                    IVs = robjects.conversion.get_conversion().py2rpy(IVs)

                exposure_IVs.append(IVs)
            except RRuntimeError as e:
                print(f'Error processing {sig_var_files[i]}: {e}')
                error_indices.append(i)
                continue
        else:
            try:
                exp_data = TwoSampleMR.read_exposure_data(filename = i, clump = False, sep = delimiter, snp_col = rsid, beta_col = beta, se_col = se, effect_allele_col = effect_allele, other_allele_col = other_allele, eaf_col = eaf, pval_col = pval, samplesize_col = samplesize, chr_col = chrom, pos_col = pos)
                
                IVs = TwoSampleMR.clump_data(exp_data, clump_r2 = clump_threshold, pop = population)
                exposure_IVs.append(IVs)
            except RRuntimeError as e:
                print(f'Error processing {sig_var_files[i]}: {e}')
                error_indices.append(i)
                continue

    for i in sorted(error_indices, reverse=True):
        del sig_var_files[i]

    updated_exposures = [re.sub(r"_full_sig_variants", "", os.path.splitext(os.path.basename(sig_var_files[i]))[0]) for i in range(len(sig_var_files))]

    return exposure_IVs, updated_exposures


def MR_analysis(exposure_IVs, updated_exposures, outcome_gwas, delimiter, rsid_outcome, beta_outcome, se_outcome, effect_allele_outcome, other_allele_outcome, eaf_outcome, pval_outcome, res_out=None, data_out=None, pleiotropy_out=None, het_out=None, plot_out=None, singleSNP_out=None):

    """
    Perform MR analysis on each exposure with the desired outcome

    Args:
        exposure_IVs: List of exposure IVs for each exposure

        outcome_gwas (str): Outcome GWAS file
    
        delimiter (str): Delimiter used in outcome GWAS
    
        rsid_outcome (str): Label used for rsID column in outcome GWAS

        beta_outcome (str): Label used for effect size column in outcome GWAS
    
        se_outcome (str): Label used for standard error of effect size column in
                          outcome GWAS
    
        effect_allele_outcome (str): Label used for effect allele column in outcome GWAS
    
        other_allele_outcome (str): Label used for other allele column in outcome GWAS
    
        eaf_outcome (str): Label used for effect allele frequency column in outcome GWAS

        pval_outcome (str): Label used for p-value column in outcome GWAS

        res_out (str): Destination directory for MR results

        data_out (str): Destination directory for harmonized data

        pleiotropy_out (str): Destination directory for pleiotropy tests

        het_out (str): Destination directory for heterogeneity tests

        plot_out (str): Destination directory for MR scatter plots

        singleSNP_out (str): Destination directory for single-SNP plots

    Output:
        6 lists (each with length=number of exposures):
    
               (a) 4 lists of pandas dataframes for each exposure containing:
                   MR results, harmonized data, pleiotropy tests, heterogeneity tests
    
               (b) 2 lists of plots (R objects): MR scatter plots, single-SNP plots

        Also outputs these products (dataframes + plots) to desired directories
    """
    #outcome_dat_list = []
    MR_results = []
    MR_data_total = []
    pleiotropy_tests = []
    het_tests = []
    p_all = []
    psingle_all = []
    for i in range(len(exposure_IVs)):
        #outcome_file = outcome_gwas
        outcome_dat = TwoSampleMR.read_outcome_data(snps = exposure_IVs[i][0], filename=outcome_gwas, sep = delimiter, snp_col = rsid_outcome, beta_col = beta_outcome, se_col = se_outcome, effect_allele_col = effect_allele_outcome, other_allele_col = other_allele_outcome, eaf_col = eaf_outcome, pval_col = pval_outcome)
        #outcome_dat_list.append(outcome_dat)

        MR_data = TwoSampleMR.harmonise_data(exposure_IVs[i], outcome_dat)
        res = TwoSampleMR.mr(MR_data)
        pleiotropy = TwoSampleMR.mr_pleiotropy_test(MR_data)
        het = TwoSampleMR.mr_heterogeneity(MR_data)
        res_single = TwoSampleMR.mr_singlesnp(MR_data)
        psingle = TwoSampleMR.mr_forest_plot(res_single)
        p = TwoSampleMR.mr_scatter_plot(res, MR_data)
        
        with (robjects.default_converter + pandas2ri.converter).context():
            res = robjects.conversion.get_conversion().rpy2py(res)
            MR_data = robjects.conversion.get_conversion().rpy2py(MR_data)
            pleiotropy = robjects.conversion.get_conversion().rpy2py(pleiotropy)
            het = robjects.conversion.get_conversion().rpy2py(het)

        res.to_csv(f'{res_out}/{updated_exposures[i]}.csv')
        MR_data.to_csv(f'{data_out}/{updated_exposures[i]}.csv')
        pleiotropy.to_csv(f'{pleiotropy_out}/{updated_exposures[i]}.csv')
        het.to_csv(f'{het_out}/{updated_exposures[i]}.csv')

        grdevices.png(file=f"{plot_out}/{updated_exposures[i]}.png", width = 1000, height = 850)
        rprint(p)
        grdevices.dev_off()

        grdevices.png(file=f"{singleSNP_out}/{updated_exposures[i]}_singleSNP.png", width = 1000, height = 850)
        rprint(psingle)
        grdevices.dev_off()
            
        MR_results.append(res)
        MR_data_total.append(MR_data)
        pleiotropy_tests.append(pleiotropy)
        het_tests.append(het)
        p_all.append(p)
        psingle_all.append(psingle)

    return MR_results, MR_data_total, pleiotropy_tests, het_tests, p_all, psingle_all

        
