"""
From multiple sets of results from regressions, take average of beta weights and p values from common covariates.
Assumes that there is a root folder with one folder for each site and that inside each site folder is a folder called results,
which contains the results as nifti files. The files containing beta weights should be named as "beta_covariate.nii", and those containing p-values should be named "pval_covariate.nii."

The path to the root folder is expected as an optional parameter. If no path is given on the command line, 
then the current working directory is used.

Usage:
    
    $ python average_results /path/to/results

Expected path structure:

        /path/to/results/
            site1/
                results/
                    beta_cov1.nii
                    beta_cov2.nii
                    ...
                    pval_cov1.nii
                    pval_cov2.nii
            site2/
                results/
                    beta_cov1.nii
                    beta_cov2.nii
                    ...
                    pval_cov1.nii
                    pval_cov2.nii
                
"""
import sys
import os
import numpy as np
import nibabel as nib
import logging


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


FNAME_DELIM = '_'
BETA = 'beta'
PVAL = 'pval'
SUFFIX_NII = '.nii'
RESULTS_DIR = 'results'
OUTPUT_DIR = 'global_avg'


def get_covariates(path):
    filenames = [name for name in os.listdir(path)]
    beta_files = [f for f in filenames if is_nifti(f) and is_beta(f)]
    return [fname2covar(f) for f in beta_files]


def fname2covar(filename):
    root,_ = os.path.splitext(filename)
    parts = root.split(FNAME_DELIM)
    return FNAME_DELIM.join(parts[1:])


def is_nifti(filename):
    return SUFFIX_NII in filename


def is_beta(filename):
    return BETA in filename

        
def is_pval(filename):
    return BETA in filename
        

def get_common_covars(site_covars):
    covar_list = [covars for site, covars in site_covars.items()]
    common = set(covar_list[0])

    for covars in covar_list[1:]:
        common = common & set(covars)

    return list(common)


def average_results(covar, path, sites, kind):
    images = list()
    for site in sites:
        full_path = os.path.join(path, site, RESULTS_DIR, get_img_fname(kind, covar)) 
        images.append(nib.load(full_path))

    return average_images(images)
    

def average_images(images):
    data_arrays = [i.get_data() for i in images]
    avg_data = average_arrays(data_arrays)
    
    affine_arrays = [i.affine for i in images]
    avg_affine = average_arrays(affine_arrays)

    return nib.Nifti1Image(avg_data, avg_affine)


def average_arrays(arrays):
    return np.mean(np.array(arrays), axis=0)


def get_img_fname(kind, covar):
    return kind + FNAME_DELIM + covar + SUFFIX_NII


def main(root_folder):

    # Get list of sites from folders
    site_list = [name for name in os.listdir(root_folder) if os.path.isdir(os.path.join(root_folder, name))]
    if OUTPUT_DIR in site_list:
        site_list.remove(OUTPUT_DIR)

    log.info('Sites: ' + ', '.join(site_list))

    if site_list:

            # Get list of covariates from each site
            covars_by_site = dict()
            for site in site_list:
                covars_by_site[site] = get_covariates(os.path.join(root_folder, site, RESULTS_DIR))
                log.info('Covariates in site `{site}`: {covars}'.format(site=site, covars=', '.join(covars_by_site[site])))

            # Get list of common covariates
            common_covars = get_common_covars(covars_by_site)
            log.info('Common covariates: ' + ', '.join(common_covars))

            output_folder = os.path.join(root_folder, OUTPUT_DIR)
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

            # For each common covariate, find average betas and pvalues and save to file
            for covar in common_covars:
                avg_beta = average_results(covar, root_folder, site_list, BETA)
                output_path = os.path.join(output_folder, get_img_fname(BETA, covar))
                avg_beta.to_filename(output_path)
                log.info('Saved average beta weights for `{covar}` to {path}'.format(covar=covar, path=output_path)) 

                avg_pval = average_results(covar, root_folder, site_list, PVAL)
                output_path = os.path.join(output_folder, get_img_fname(PVAL, covar))
                avg_pval.to_filename(output_path)
                log.info('Saved average p-values for `{covar}` to {path}'.format(covar=covar, path=output_path)) 

    else:
        log.error('No sites found')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    else:
        folder = os.get_cwd()
    log.info('Root folder: ' + folder)
    main(folder)
