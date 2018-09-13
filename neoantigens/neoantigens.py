#!/usr/bin/env python
import os
import sys
from os.path import isfile, join, dirname, abspath
import click
import subprocess

from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import verify_file, safe_mkdir, verify_dir
from ngs_utils import logger
from ngs_utils.logger import warn, info
from ngs_utils.utils import set_locale; set_locale()
from ngs_utils import snakemake_utils
from neoantigens import package_path


@click.command()
@click.argument('dna_bcbio', type=click.Path(exists=True))
@click.argument('dna_sample')
@click.argument('target_rule', nargs=-1)
@click.option('-o', 'output_dir', type=click.Path(), help='Output directory (default is "umccrise")')
@click.option('-R', '--rna-bcbio', type=click.Path())
@click.option('-r', '--rna-sample')

@click.option('-j', '--jobs', 'jobs', default=1, help='Maximum number of cores to use at single time (works both for '
              'local and cluster runs)')
@click.option('-c', '--cluster-auto', 'cluster', is_flag=True, help='Submit jobs to cluster')
@click.option('--unlock', is_flag=True, help='Propagaded to snakemake')
@click.option('--restart-times', default=1, help='Propagaded to snakemake.')

def main(dna_bcbio, dna_sample, target_rule=list(), output_dir=None, rna_bcbio=None, rna_sample=None,
         jobs=None, cluster=False, unlock=False, restart_times=None,
    ):
    output_dir = output_dir or 'neoantigens'
    output_dir = safe_mkdir(abspath(output_dir))
    log_dir = safe_mkdir(join(output_dir, 'log'))
    logger.init(log_fpath_=join(log_dir, 'neoantigens.log'), save_previous=True)

    conf = dict()
    conf['dna_bcbio'] = verify_dir(dna_bcbio, is_critical=True)
    conf['dna_sample'] = dna_sample
    conf['rna_bcbio'] = rna_bcbio
    conf['rna_sample'] = rna_sample

    #########################
    #### Setting cluster ####
    #########################

    cluster_log_dir = ''
    if cluster:
        cluster_log_dir = safe_mkdir(join(log_dir, 'cluster'))
        cluster_param = snakemake_utils.make_cluster_cmdl(cluster_log_dir, 'neoantigens')

    ###############################
    #### Building command line ####
    ###############################

    target_rule = list(target_rule)
    snakefile = join(package_path(), 'Snakefile')
    cmd = (
        f'snakemake '
        f'{" ".join(target_rule)} '
        f'--snakefile {snakefile} '
        f'--printshellcmds '
        f'--directory {output_dir} '
        f'-j {jobs} '
        f'--rerun-incomplete ' 
        f'--restart-times {restart_times} '
        f'{cluster_param} '
        f'--config {" ".join(k + "=" + v for k, v in conf.items())} '
    )

    #################
    #### Running ####
    #################

    if unlock:
        print('* Unlocking previous run... *')
        run_simple(cmd + ' --unlock')
        print('* Now rerunning *')

    try:
        run_simple(cmd)
    except subprocess.CalledProcessError:
        logger.error('--------')
        logger.error(f'Error running Umccrise: snakemake returned a non-zero status. Working directory: {output_dir}')
        run_simple(f'chmod a+r {log_dir}')
        if cluster_log_dir:
            logger.error(f'Review cluster job logs in {cluster_log_dir}')
        sys.exit(1)
    else:
        logger.error('--------')
        run_simple(f'chmod a+r {log_dir}')
        logger.info(f'Finished. Output directory: {output_dir}')


if __name__ == '__main__':
    main()
