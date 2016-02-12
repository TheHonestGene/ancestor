"""
    ancestor
    ~~~~~~~~~~~~~
    The main module for running Risk Predictions
    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""
import os
import sys
import argparse
import logging, logging.config
import h5py
from .core import ancestry as an

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'default': {
            'format': '%(asctime)s %(levelname)s %(name)s %(message)s'
        },
    },
    'handlers': {
        'stdout': {
            'class': 'logging.StreamHandler',
            'stream': 'ext://sys.stdout',
            'formatter': 'default',
        },
        'stderr': {
            'class': 'logging.StreamHandler',
            'stream': 'ext://sys.stderr',
            'level': 'ERROR',
            'formatter': 'default',
        },
    },
    'root': {
        'handlers': ['stdout', 'stderr'],
        'level': 'INFO',
    },
}

logging.config.dictConfig(LOGGING)
log = logging.getLogger()


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot", dest="plot_file", metavar="plot file", help="If specified, PCs will saved in a plot")
    parser.add_argument('--hapmap', dest='hapmap_file', metavar="hapmap file",
                        help='Hapmap dataset file in HDF5 format. Only required if --pcs does not exist or is not specified')
    parser.add_argument('--pcs', dest='pcs_file', metavar="PCs file", help='The stored PCs file')
    parser.add_argument(dest='genotype_file', metavar="genotype-file",
                        help='Genotype file to check ancestry for (HDF5 format)')
    parser.add_argument(dest='weights_file', metavar="weights-file",
                        help='File that contains the PCs weights.')

    return parser


def main():
    # Process arguments
    parser = get_parser()
    args = vars(parser.parse_args())
    try:
        if args['pcs_file'] is None and args['hapmap_file'] is None:
            parser.error('At a minimum --pcs or --hapmap must be specified')
        if args['pcs_file'] is not None and not os.path.exists(args['pcs_file']):
            if args['hapmap_file'] is None:
                parser.error('If pcs_file does not exist you have to provide --hapmap')
        run(args)
        return 0
    except KeyboardInterrupt:
        return 0
    except Exception as e:
        log.exception(e)
        return 2


def run(args):
    genotype_file = args['genotype_file']
    weight_dict,stats = an.parse_pc_weights(args['weights_file'])
    if args['pcs_file'] and os.path.exists(args['pcs_file']):
        hapmap_pcs_dict = an.load_pcs_from_file(args['pcs_file'])
    else:
        snp_ids = _get_snp_ids_from_genotype(genotype_file)
        hapmap_pcs_dict = an.calculate_hapmap_pcs(args['hapmap_file'], weight_dict,snp_ids)
        hapmap_pcs_dict['populations'] = _get_populations_from_hapmap(args['hapmap_file'])
        if args['pcs_file'] is not None:
            an.save_hapmap_pcs(hapmap_pcs_dict['pcs'],hapmap_pcs_dict['populations'], args['pcs_file'])
    genotype_pcs = an.calc_genotype_pcs(genotype_file, weight_dict)

    pcs = hapmap_pcs_dict['pcs']
    populations = hapmap_pcs_dict['populations']
    eur_filter = populations['EUR']
    pc1 = genotype_pcs['pc1']
    pc2 = genotype_pcs['pc2']
    ancestry_dict = an.check_european(pcs[eur_filter], pc1, pc2)
    eur_mean_PC1, eur_mean_PC2 = ancestry_dict['eur_mean'].tolist()
    eur_std_PC1, eur_std_PC2 = ancestry_dict['eur_std'].tolist()
    eur_lim_PC1, eur_lim_PC2 = ancestry_dict['eur_lim'].tolist()
    ind_lim_PC1, ind_lim_PC2 = ancestry_dict['ind_lim'].tolist()

    log.debug('European mean: PC1: %.6f, PC2: %.6f and std: PC1 :%.6f, PC2: %.6f' % (
        eur_mean_PC1, eur_mean_PC2, eur_std_PC1, eur_std_PC2))
    log.debug('PC1: %.6f and PC2: %.6f' % (pc1, pc2))

    if ancestry_dict['is_non_european']:
        log.warn(
                'Sample appears to contain some non-European ancestry (%.6f,%.6f > %.6f,%.6f).  This means that current genetic predictions are probably not accurate.' % (
                    ind_lim_PC1, ind_lim_PC2, eur_lim_PC1, eur_lim_PC2))
    else:
        log.info('Sample appears to be of European ancestry (%.6f,%.6f < %.6f,%.6f).' % (
            ind_lim_PC1, ind_lim_PC2, eur_lim_PC1, eur_lim_PC2))

    if args['plot_file'] is not None:
        an.plot_pcs(args['plot_file'], pcs, populations, genotype_pcs)


def _get_snp_ids_from_genotype(genotype):
    f = h5py.File(genotype, 'r')
    snp_ids = {}
    for chr in f.keys():
        snp_ids[chr.lower()] = f[chr]['sids'][:]
    return snp_ids


def _get_populations_from_hapmap(hapmap_file):
    f = h5py.File(hapmap_file, 'r')
    populations = {}
    for population in ['EUR', 'SAS', 'AMR', 'AFR', 'EAS']:
        populations[population] = f['indivs']['continent'][...] == population
    f.close()
    return populations


if __name__ == '__main__':
    sys.exit(main())
