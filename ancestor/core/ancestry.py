"""
core functions to calculate ancestry
"""
import scipy as sp
import numpy as np
import h5py
import pylab
import logging
from os import path

log = logging.getLogger(__name__)


def parse_pc_weights(pc_weights_file):
    """
    Parses the weights from a file

    :param pc_weights_file: the file that contains the SNP weights
    """
    if path.basename(pc_weights_file).lower().endswith((".hdf5", ".h5")):
        return _parse_pc_weights_from_hdf5(pc_weights_file)
    return _parse_pc_weights_from_text(pc_weights_file)


def _parse_pc_weights_from_text(pc_weights_file):
    sid_dict = {}
    log.info('Parsing weights file %s' % pc_weights_file)
    with open(pc_weights_file) as f:
        shrinkage = map(float, f.readline().split())
        populations = {pop: {'count': count} for (pop, count) in
                       zip(f.readline().split(), map(int, f.readline().split()))}
        l = f.readline().split()
        for i in range(0, len(l), 3):
            populations[l[i]]['meanpc'] = map(float, [l[i + 1], l[i + 2]])
        linear_transform = map(float, f.readline().split())
        for line in f:
            l = line.split()
            sid = l[0]
            nts = [l[1].encode('latin-1'), l[2].encode('latin-1')]
            mean_g = np.float32(l[3])
            pc1w = np.float32(l[4])
            pc2w = np.float32(l[5])
            sid_dict[sid.encode('latin-1')] = {'pc1w': pc1w, 'pc2w': pc2w, 'mean_g': mean_g, 'nts': nts}

    return sid_dict, {'populations': populations, 'shrinkage': shrinkage, 'linear_transform': linear_transform}


def _parse_pc_weights_from_hdf5(pc_weights_file):
    gh5f = h5py.File(pc_weights_file, 'r')
    weights = {}
    attrs = gh5f.attrs;
    stats = {'linear_transform': attrs['linear_transform'].tolist(), 'shrinkage': attrs['shrinkage'].tolist()}
    populations = {}
    for key, value in attrs.items():
        if key in ('linear_transform', 'shrinkage'):
            continue
        key_s = key.split('_')
        if len(key_s) != 2:
            continue
        population = key_s[0]
        if population not in populations:
            populations[population] = {}
        populations[population][key_s[1]] = value.tolist()
    stats['populations'] = populations
    snpids = gh5f['snpid'][:]
    nts = gh5f['nts'][:]
    mean_g = gh5f['mean_g'][:]
    pc1w = gh5f['pc1w'][:]
    pc2w = gh5f['pc2w'][:]
    for i, snpid in enumerate(snpids):
        weights[snpid] = {'pc1w': pc1w[i], 'pc2w': pc2w[i], 'mean_g': mean_g[i], 'nts': nts[i].tolist()}
    gh5f.close()
    return weights, stats


def save_pc_weights(weights, stats, output_file):
    """
    Save the weights and stats in an HDF5 file

    :param weights: dictionary with weights
    :param stats: general statistics from the weights file
    :param output_file: hdf5 file to store the data in
    """
    gh5f = h5py.File(output_file, 'w')
    snpids = []
    nts = []
    mean_g = []
    pc1w = []
    pc2w = []
    gh5f.attrs['linear_transform'] = stats['linear_transform']
    gh5f.attrs['shrinkage'] = stats['shrinkage']
    for population, value in stats['populations'].items():
        gh5f.attrs['%s_count' % population] = value['count']
        gh5f.attrs['%s_meanpc' % population] = value['meanpc']

    for snpid, value in weights.items():
        snpids.append(snpid)
        nts.append(value['nts'])
        mean_g.append(value['mean_g'])
        pc1w.append(value['pc1w'])
        pc2w.append(value['pc2w'])
    gh5f.create_dataset('snpid', data=snpids, chunks=True, compression='lzf')
    gh5f.create_dataset('nts', data=np.asarray(nts), chunks=True, compression='lzf')
    gh5f.create_dataset('mean_g', data=np.asarray(mean_g, dtype='f4'), chunks=True, compression='lzf')
    gh5f.create_dataset('pc1w', data=np.asarray(pc1w, dtype='f4'), chunks=True, compression='lzf')
    gh5f.create_dataset('pc2w', data=np.asarray(pc2w, dtype='f4'), chunks=True, compression='lzf')
    gh5f.close()


def calc_genotype_pcs(genotype_file, weight_dict):
    """
    Calculate the principal components for a given genotype using the specified weights

    :param genotype_file: the genotype file for which the principal components should be calculated
    :param weight_dict: dictionary with SNP weights
    """
    gh5f = h5py.File(genotype_file, 'r')
    num_nt_issues = 0
    num_snps_used = 0
    log.info('Calculating Principal Components for genotype file %s' % genotype_file)
    pcs = sp.zeros((1, 2))
    for chrom in range(1, 23):
        log.info('Working on Chromosome %d' % chrom)

        chrom_str = 'Chr%d' % chrom
        g_cg = gh5f[chrom_str]

        # Loading data
        sids = g_cg['sids'][...]
        nts = g_cg['nts'][...]
        length = len(g_cg['snps'])
        snps = g_cg['snps'][...].reshape((length, 1))
        pcs_per_chr = _calc_pcs(weight_dict, sids, nts, snps)
        pcs += pcs_per_chr['pcs']
        num_snps_used += pcs_per_chr['num_snps_used']
        num_nt_issues += pcs_per_chr['num_nt_issues']
    gh5f.close()
    log.info('%d SNPs were excluded from the analysis due to nucleotide issues.' % (num_nt_issues))
    log.info('Number of SNPs uesd for PC projection: %d' % num_snps_used)
    return {'pc1': pcs[0][0], 'pc2': pcs[0][1], 'num_snps_used': num_snps_used}


def calculate_hapmap_pcs(hapmap_file, pc_weights_dict, snps_filter=None):
    """
    Calculates the principal components for the hapmap project

    :param hapmap_file: Hapmap file in HDF5 format
    :param pc_weights_dict: dictionary with SNP weights (key = snpid)
    :param snps_filter: list of snp-ids to subset (optional)
    :return: dictionary with pcs and number of snps that were used
    """
    assert not pc_weights_file or not pc_weights_file

    log.info('Calculating Principal components for Hapmap file %s' % hapmap_file)
    ok_sids = sp.array(pc_weights_dict.keys())
    log.info('Loaded PC weight for %d SNPs' % (len(ok_sids)))
    # Load genotypes
    log.info('Load Hapmap dataset')
    h5f = h5py.File(hapmap_file, 'r')
    num_indivs = len(h5f['indivs']['continent'][...])
    log.info('Found genotypes for %d individuals' % num_indivs)
    pcs = sp.zeros((num_indivs, 2))
    num_nt_issues = 0
    num_snps_used = 0
    log.info('Calculating PCs')
    for chrom in range(1, 23):
        log.info('Working on Chromosome %d' % chrom)
        chrom_str = 'chr%d' % chrom

        log.info('Identifying overlap')
        ok_snp_filter = sp.in1d(ok_sids, snps_filter[chrom_str])
        ok_chrom_sids = ok_sids.compress(ok_snp_filter, axis=0)
        sids = h5f[chrom_str]['variants']['ID'][...]
        ok_snp_filter = sp.in1d(sids, ok_chrom_sids)
        #         assert sids[ok_snp_filter]==ok_sids, 'WTF?'
        sids = sids.compress(ok_snp_filter, axis=0)

        log.info('Loading SNPs')
        snps = h5f[chrom_str]['calldata']['snps'][...]
        snps = snps.compress(ok_snp_filter, axis=0)
        length = len(h5f[chrom_str]['variants/REF'])
        nts = np.hstack((h5f[chrom_str]['variants/REF'][:].reshape(length, 1),
                         h5f[chrom_str]['variants/ALT'][:].reshape(length, 1)))
        nts = nts.compress(ok_snp_filter, axis=0)
        log.info('Updating PCs')
        pcs_per_chr = _calc_pcs(pc_weights_dict, sids, nts, snps)
        pcs += pcs_per_chr['pcs']
        num_nt_issues += pcs_per_chr['num_nt_issues']
        num_snps_used += pcs_per_chr['num_snps_used']

    h5f.close()
    log.info('%d SNPs were excluded from the analysis due to nucleotide issues.' % (num_nt_issues))
    log.info('%d SNPs were used for the analysis.' % (num_snps_used))
    return {'pcs': pcs, 'num_snps_used': num_snps_used}


def save_hapmap_pcs(pcs, populations, output_file):
    """
    Saves the hapmap pcs in an HDF5 file
    :param output_file: HDF5 file the pcs should be written to
    :param pcs: principal components
    :param populations: dictionary with various populations masks
    """
    log.info('Saving HapMap PCs in %s ' % output_file)
    # Store coordinates
    oh5f = h5py.File(output_file, 'w')
    oh5f.create_dataset('pcs', data=pcs)
    pop_g = oh5f.create_group('populations')
    for population, mask in populations.items():
        pop_g.create_dataset(population, data=mask)
    oh5f.close()


def load_pcs_from_file(input_file):
    """
    Loads pcs from an HDF5 file
    :param input_file: HDF5 file that contains the pcs
    :return: Dictionary with pcs and individual filters
    """
    log.info('Loading HapMap PCs from %s ' % input_file)
    ch5f = h5py.File(input_file)
    pop_g = ch5f['populations']
    populations = {}
    for key in pop_g.keys():
        populations[key] = pop_g[key][...]
    pcs = ch5f['pcs'][...]
    ch5f.close()
    return {'populations': populations,
            'pcs': pcs}


def check_european(eur_pcs, pc1, pc2):
    """
    Check if PC1 and PC2 are within the european PCs
    :param eur_pcs:
    :param pc1: PC1 of the individual
    :param pc2: PC2 of the individual
    :return: Dictionary with various statistics
    """
    # Report ancestry.
    eur_mean = sp.mean(eur_pcs, 0)
    eur_stds = sp.std(eur_pcs, 0)
    eur_lim = (3 * eur_stds) ** 2
    ind_pcs = sp.array([pc1, pc2])
    ind_lim = (ind_pcs - eur_mean) ** 2
    is_non_european = sp.any(ind_lim ** 2 > eur_lim)
    return {'eur_lim': eur_lim, 'eur_mean': eur_mean, 'eur_std': eur_stds, 'ind_lim': ind_lim,
            'is_non_european': is_non_european}


def plot_pcs(plot_file, pcs, populations, genotype_pcs_dict=None):
    """
    Plots the PCs of the hapmap and if provided of the genotype
    :param populations: dictionary with different population masks for coloring the individuals
    :param plot_file: Name of the file the plot should be saved (either .png or .pdf)
    :param pcs: principal components from the hapmap dataset
    :param genotype_pcs_dict: Dictionary with PCs of the individual (optional)
    """
    log.info('Plotting PCs of Hapmap')
    # Plot them
    for population, mask in populations.items():
        pylab.plot(pcs[mask][:, 0], pcs[mask][:, 1], label=population, ls='', marker='.', alpha=0.6)

    log.info('Plotting genome on plot')
    # Project genome on to plot.
    if genotype_pcs_dict is not None:
        pylab.plot(genotype_pcs_dict['pc1'], genotype_pcs_dict['pc2'], 'o', label='This is you')
    pylab.xlabel('PC 1')
    pylab.ylabel('PC 2')
    pylab.legend(loc=4, numpoints=1)
    pylab.savefig(plot_file)


def _calc_pcs(weight_dict, sids, nts, snps):
    num_nt_issues = 0
    num_snps_used = 0
    num_indivs = snps.shape[1]
    pcs = sp.zeros((num_indivs, 2))

    for snp_i, sid in enumerate(sids):
        try:
            d = weight_dict[sid]
        except:
            continue
        nt = nts[snp_i]
        snp = snps[snp_i]
        pc_weights = sp.array([d['pc1w'], d['pc2w']])
        pc_weights.shape = (1, 2)
        af = d['mean_g'] / 2.0
        if sp.all([nt[1], nt[0]] == d['nts']):
            # print 'Flip sign'
            pc_weights = -pc_weights
            af = 1 - af
        elif not sp.all(nt == d['nts']):
            num_nt_issues += 1
            continue
        mean_g = 2 * af
        sd_g = sp.sqrt(af * (1 - af))
        # "Normalizing" the SNPs with the given allele frequencies
        norm_snp = (snp - mean_g) / sd_g
        norm_snp.shape = (num_indivs, 1)

        # Project on the PCs
        pcs += sp.dot(norm_snp, pc_weights)
        num_snps_used += 1

    return {'num_snps_used': num_snps_used, 'num_nt_issues': num_nt_issues, 'pcs': pcs}
