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
from os import path
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
    subparsers = parser.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')

    convert_parser = subparsers.add_parser('convert',help='Converting weights file from CSV format to HDF5 format')
    convert_parser.add_argument(dest="input_file", help="Original weights file (CSV format)")
    convert_parser.add_argument(dest="output_file", help="Name of the output file (must have .h5 or .hdf5 extension)")
    convert_parser.set_defaults(func=convert_snp_weights)

    generate_pcs_parser = subparsers.add_parser('prepare',help='Calculating PC projections and admixture decomposition information for a reference panel')
    generate_pcs_parser.add_argument(dest="input_file", help="Reference genotype panel file (in HDF5 format)")
    generate_pcs_parser.add_argument(dest="weights_file", help="Weights file (in HDF5 format or CSV format)")
    generate_pcs_parser.add_argument(dest="output_file", help="Name of the output file")
    generate_pcs_parser.add_argument("--ntmap", dest="nt_map_file", metavar="nt-map file", help="Name of the nucleotide map file to subset the SNPs")
    generate_pcs_parser.set_defaults(func=generate_pcs)


    ancestry_parser = subparsers.add_parser('pcs',help='Calculate the principal components for a given individual genotype using the specified weights')
    ancestry_parser.add_argument(dest="genotype_file", metavar="genotype-file", help="Genotype file to check ancestry for (HDF5 format)")
    ancestry_parser.add_argument(dest='weights_file', metavar="weights-file", help='Weights file (CSV or HDF5 format)')
    ancestry_parser.add_argument(dest='pcs_file', metavar="pcs-file", help='File with PCs for reference genotype panel (HDF5 format).')
    ancestry_parser.add_argument("--plot", dest="plot_file", metavar="plot file", help="If specified, PCs will saved in a plot")
    ancestry_parser.add_argument("--check",dest="check_population",help="Specify a population (i.e. EUR) to check if the individual genotype is part of")
    ancestry_parser.set_defaults(func=run)
    return parser


def main():
    # Process arguments
    parser = get_parser()
    args = vars(parser.parse_args())
    if 'func' not in args:
        parser.print_help()
        return 0
    try:
        args['func'](args)
        return 0
    except KeyboardInterrupt:
        return 0
    except Exception as e:
        log.exception(e)
        return 2


def run(args):
    genotype_file = args['genotype_file']
    ref_pcs_admix_file = args['pcs_file']
    check_population = args.get('check_population',None)
    weight_dict,stats = an.parse_pc_weights(args['weights_file'])
    pcs_dict = an.load_pcs_admixture_info(ref_pcs_admix_file)
    ancestry_results = an.ancestry_analysis(genotype_file, args['weights_file'], ref_pcs_admix_file, check_population=check_population)
    indiv_pcs = ancestry_results['indiv_pcs']
    pc1 = indiv_pcs[0]
    pc2 = indiv_pcs[1]
    log.info(ancestry_results)

    if args['plot_file'] is not None:
        an.plot_pcs(args['plot_file'], pcs_dict['pcs'], pcs_dict['pop_dict']['populations'], indiv_pcs=indiv_pcs)

def convert_snp_weights(args):
    input_pc_weights_file = args['input_file']
    output_pc_weights_file = args['output_file']
    if not path.basename(output_pc_weights_file).lower().endswith((".hdf5", ".h5")):
        raise Exception('% file must have extension .hdf5 or .h5' % output_pc_weights_file)
    if not path.basename(input_pc_weights_file).lower().endswith((".hdf5", ".h5")):
        sid_weights_map, stats_dict = an.parse_pc_weights(input_pc_weights_file)
        an.save_pc_weights(sid_weights_map, stats_dict, output_pc_weights_file)


def generate_pcs(args):
    nt_map_file = args.get('nt_map_file', None)
    ref_pcs_admix_file = args['output_file']
    pc_ref_genot_file = args['input_file']
    populations_to_use = ['TSI', 'FIN', 'IBS', 'GBR']
    pc_weights_file = args['weights_file']
    snp_filter = None
    if nt_map_file is not None:
        snps_filter = an.get_snps_filter(nt_map_file)

    sid_weights_map, stats_dict = an.parse_pc_weights(pc_weights_file)
    pcs_dict = an.calc_genot_pcs(pc_ref_genot_file, sid_weights_map, stats_dict, populations_to_use = populations_to_use, snps_filter=snps_filter)
    an.save_pcs_admixture_info(pcs_dict['pcs'], pcs_dict['pop_dict'], ref_pcs_admix_file)



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
        populations[population] = f['indivs']['continent'][...] == population.encode() # py3 fix
    f.close()
    return populations


if __name__ == '__main__':
    sys.exit(main())
