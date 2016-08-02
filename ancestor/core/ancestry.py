"""
core functions to calculate ancestry
"""
import scipy as sp
from scipy import linalg
import numpy as np
import h5py
import pylab
import cPickle
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
        shrinkage = list(map(float, f.readline().split()))
        populations = {pop: {'count': count} for (pop, count) in
                       zip(f.readline().split(), list(map(int, f.readline().split())))}
        l = f.readline().split()
        num_pops = len(populations)
        for i in range(0, len(l), num_pops):
            populations[l[i]]['meanpc'] = list(map(float, l[i+1:i+num_pops]))
        linear_transform = list(map(float, f.readline().split()))
        for line in f:
            l = line.split()
            sid = l[0]
            nts = [l[1].encode('latin-1'), l[2].encode('latin-1')]
            mean_g = np.float32(l[3])
            pc_ws = list(map(float,(l[4:])))
            sid_dict[sid.encode('latin-1')] = {'pc_ws': pc_ws, 'mean_g': mean_g, 'nts': nts}
        num_pcs = len(pc_ws)
        assert num_pcs == num_pops-1, 'WTF?'
    return sid_dict, {'populations': populations, 'shrinkage': shrinkage, 'linear_transform': linear_transform, 'num_pcs':num_pcs}


def _parse_pc_weights_from_hdf5(pc_weights_file):
    """
    Loads the SNP weights dict!
    """
    gh5f = h5py.File(pc_weights_file, 'r')
    weights = {}
    attrs = gh5f.attrs;
    stats_dict = {'linear_transform': attrs['linear_transform'].tolist(), 'shrinkage': attrs['shrinkage'].tolist()}
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
    stats_dict['populations'] = populations
    snpids = gh5f['snpid'][...]
    nts = gh5f['nts'][...]
    mean_g = gh5f['mean_g'][...]
    pc_ws = gh5f['pc_ws'][...]
    stats_dict['num_pcs'] = len(pc_ws[0]) 
    for i, snpid in enumerate(snpids):
        weights[snpid] = {'pc_ws': pc_ws[i], 'mean_g': mean_g[i], 'nts': nts[i].tolist()}
    gh5f.close()
    return weights, stats_dict


def save_pc_weights(weights, stats, output_file):
    """
    Save the SNP weights and stats in an HDF5 file

    :param weights: dictionary with weights
    :param stats: general statistics from the weights file
    :param output_file: hdf5 file to store the data in
    """
    gh5f = h5py.File(output_file, 'w')
    snpids = []
    nts = []
    mean_g = []
    pc_ws = []
    gh5f.attrs['linear_transform'] = stats['linear_transform']
    gh5f.attrs['shrinkage'] = stats['shrinkage']
    for population, value in stats['populations'].items():
        gh5f.attrs['%s_count' % population] = value['count']
        gh5f.attrs['%s_meanpc' % population] = value['meanpc']

    for snpid, value in weights.items():
        snpids.append(snpid)
        nts.append(value['nts'])
        mean_g.append(value['mean_g'])
        pc_ws.append(value['pc_ws'])
    gh5f.create_dataset('snpid', data=snpids, chunks=True, compression='lzf')
    gh5f.create_dataset('nts', data=np.asarray(nts), chunks=True, compression='lzf')
    gh5f.create_dataset('mean_g', data=np.asarray(mean_g, dtype='f4'), chunks=True, compression='lzf')
    gh5f.create_dataset('pc_ws', data=np.asarray(pc_ws, dtype='f4'), chunks=True, compression='lzf')
    gh5f.close()


def calc_indiv_genot_pcs(genotype_file, weight_dict, num_pcs_to_uses,**kwargs):
    """
    Calculate the principal components for a given an individual genotype using the specified weights (
 
    :param genotype_file: the genotype file for which the principal components should be calculated
    :param weight_dict: dictionary with SNP weights
    """
    log_extra = kwargs.get('log_extra',{'progress':0})
    partial_progress_inc = (100-log_extra['progress'])/22
     
    gh5f = h5py.File(genotype_file, 'r')
    num_nt_issues = 0
    num_snps_used = 0
    log.info('Calculating Principal Components for genotype file %s' % genotype_file)
    pcs = sp.zeros((1, 2))
    for chrom in range(1, 23):
        log_extra['progress']+=partial_progress_inc
        log.info('Working on Chromosome %d' % chrom,extra=log_extra)
 
        chrom_str = 'Chr%d' % chrom
        g_cg = gh5f[chrom_str]
 
        # Loading data
        sids = g_cg['sids'][...]
        nts = g_cg['nts'][...]
        length = len(g_cg['snps'])
        snps = g_cg['snps'][...].reshape((length, 1))
        pcs_per_chr = _calc_pcs(weight_dict, sids, nts, snps, num_pcs_to_uses)
        pcs += pcs_per_chr['pcs']
        num_snps_used += pcs_per_chr['num_snps_used']
        num_nt_issues += pcs_per_chr['num_nt_issues']
        
    gh5f.close()
    log.info('%d SNPs were excluded from the analysis due to nucleotide issues.' % (num_nt_issues))
    log.info('Number of SNPs used for PC projection: %d' % num_snps_used)
    return {'pcs': pcs, 'num_snps_used': num_snps_used}




def calc_genot_pcs(genot_file, pc_weights_dict, pc_stats, populations_to_use = ['EUR','AFR','EAS'], snps_filter=None, 
                   verbose=False, debug_cutoff=None):
    """
    Calculates:
        - The principal components for the given genotype dataset.
        - The admixture decomposition matrix given the populations.

    :param genot_file: Genotype file in HDF5 format
    :param pc_weights_dict: dictionary with SNP weights (key = snpid)
    :param pc_stats: statistics about the PC SNP weights. 
    :param populations_to_use: Populations that we plan to use for the PC plot and the admixture decomposition. 
    :param snps_filter: list of snp-ids to subset (optional)
    :return: dictionary with pcs and number of snps that were used
    
    """
    assert len(populations_to_use)<=pc_stats['num_pcs']+1, 'There are not sufficiently many PCs to separate the populations.'
    
    log.info('Calculating Principal components for genotype file %s' % genot_file)
    ok_sids = np.asarray(list(pc_weights_dict.keys()))
    log.info('Loaded PC weight for %d SNPs' % (len(ok_sids)))
    # Load genotypes
    log.info('Loading genotypes')
    h5f = h5py.File(genot_file, 'r')

#     populations = h5f['indivs']['populations'][...]
    populations = h5f['indivs']['continent'][...]  #Currently this is using continents, but we can/should switch to populations.  
    #We should also consider to allow us to combine multiple populations in one group, or something like that. 
    
    indiv_filter = sp.in1d(populations, populations_to_use)
    filtered_populations = populations[indiv_filter]
    num_indivs = sp.sum(indiv_filter)
    log.info('Found genotypes for %d individuals' % num_indivs)
    num_pcs_to_use=len(populations_to_use)-1
    pcs = sp.zeros((num_indivs, num_pcs_to_use))
    num_nt_issues = 0
    num_snps_used = 0
    log.info('Calculating PCs')
    for chrom in range(1, 23):
        if verbose:
            print 'Working on Chromosome %d' % chrom
        log.info('Working on Chromosome %d' % chrom)
        chrom_str = 'chr%d' % chrom

        log.info('Identifying overlap')
        ok_snp_filter = sp.in1d(ok_sids, snps_filter[chrom_str])
        ok_chrom_sids = ok_sids.compress(ok_snp_filter, axis=0)
        if debug_cutoff is None:
            sids = h5f[chrom_str]['variants']['ID'][...]
        else:
            sids = h5f[chrom_str]['variants']['ID'][:debug_cutoff]  #A debugging hack
        ok_snp_filter = sp.in1d(sids, ok_chrom_sids)
        #         assert sids[ok_snp_filter]==ok_sids, 'WTF?'
        sids = sids.compress(ok_snp_filter, axis=0)

        log.info('Loading SNPs')
        if verbose:
            print 'Loading SNPs'
        if debug_cutoff is None:
            snps = h5f[chrom_str]['calldata']['snps'][...]
        else:
            snps = h5f[chrom_str]['calldata']['snps'][:debug_cutoff]  #A debugging hack
        if verbose:
            print 'Filtering SNPs'
        snps = snps.compress(ok_snp_filter, axis=0)
        if verbose:
            print 'Using %d SNPs'%sp.sum(ok_snp_filter)
            print 'Filtering individuals.'
        snps = snps[:,indiv_filter] #Filter individuals with certain ancestry for the analyses
        if verbose:
            print 'Using %d individuals'%sp.sum(indiv_filter)
        
        if debug_cutoff is None:
            length = len(h5f[chrom_str]['variants/REF'])
            nts = np.hstack((h5f[chrom_str]['variants/REF'][:].reshape(length, 1),
                             h5f[chrom_str]['variants/ALT'][:].reshape(length, 1)))
        else:
            length = debug_cutoff
            nts = np.hstack((h5f[chrom_str]['variants/REF'][:debug_cutoff].reshape(length, 1),
                             h5f[chrom_str]['variants/ALT'][:debug_cutoff].reshape(length, 1)))
            nts = nts.compress(ok_snp_filter, axis=0)
        log.info('Updating PCs')
        if verbose:
            print 'Calculating PC projections'
        pcs_per_chr = _calc_pcs(pc_weights_dict, sids, nts, snps, num_pcs_to_use)
        pcs += pcs_per_chr['pcs']
        num_nt_issues += pcs_per_chr['num_nt_issues']
        num_snps_used += pcs_per_chr['num_snps_used']

    
    log.info('%d Calculating average PC values for each population.' % (num_nt_issues))
    avg_pcs_list= []
    num_indiv_list = []
    num_pops_to_use = len(populations_to_use)
    E = sp.empty((num_pops_to_use,num_pops_to_use),dtype='float32')  #The avg. pop. PCs matrix, used for admixture decomposition.
    for i in range(num_pops_to_use):
        popul = populations_to_use[i]
        pop_filter = sp.in1d(filtered_populations, [popul])
        pop_num_indivs = sp.sum(pop_filter)
        avg_pcs = sp.mean(pcs[pop_filter],0)
        avg_pcs_list.append(avg_pcs)
        num_indiv_list.append(pop_num_indivs)
        E[i]=sp.concatenate((avg_pcs,[1.0]))
#     E = sp.transpose(E)
    print E 
    #For decomposition of admixture, we assume that the same set of SNPs are used.
    pop_dict = {'admix_decom_mat': linalg.inv(E), 'populations': filtered_populations, 'unique_populations':populations_to_use, 
                'avg_pcs':sp.array(avg_pcs_list), 'num_indivs':num_indiv_list}  
    h5f.close()
    log.info('%d SNPs were excluded from the analysis due to nucleotide issues.' % (num_nt_issues))
    log.info('%d SNPs were used for the analysis.' % (num_snps_used))
    
    return {'pcs': pcs, 'num_snps_used': num_snps_used, 'pop_dict':pop_dict}


def save_pcs_admixture_info(pcs, pop_dict, output_file):
    """
    Saves the PCs and admixture information to a HDF5 file
    :param output_file: HDF5 file the pcs should be written to
    :param pcs: principal components
    :param pop_dict: dictionary with various populations information
    """
    log.info('Saving genotype PCs in %s ' % output_file)
    # Store coordinates
    oh5f = h5py.File(output_file, 'w')
    oh5f.create_dataset('pcs', data=pcs)
    pop_g = oh5f.create_group('pop_group')
    for label, mask in pop_dict.items():
        pop_g.create_dataset(label, data=mask)
#     pop_g.create_dataset('admix_decom_mat', data=pop_dict['admix_decom_mat'])
#     pop_g.create_dataset('populations', data=pop_dict['populations'])
#     pop_g.create_dataset('unique_populations', data=pop_dict['unique_populations'])
#     pop_g.create_dataset('avg_pcs', data=pop_dict['avg_pcs'])
#     pop_g.create_dataset('num_indivs', data=pop_dict['num_indivs'])
    oh5f.close()


def load_pcs_admixture_info(input_file):
    """
    Loads pcs from an HDF5 file
    :param input_file: HDF5 file that contains the pcs
    :return: Dictionary with pcs and individual filters
    """
    log.info('Loading HapMap PCs from %s ' % input_file)
    ch5f = h5py.File(input_file)
    pop_g = ch5f['pop_group']
    pop_dict = {}
    for key in pop_g.keys():
        pop_dict[key] = pop_g[key][...]
    pcs = ch5f['pcs'][...]
    ch5f.close()
    return {'pop_dict': pop_dict, 'pcs': pcs}


def ancestry_analysis(genotype_file, weights_file, pcs_file, check_population='EUR', 
                      verbose=False, debug_cutoff=None, **kwargs):
    """
    Runs the ancestry analysis on a single genome.  It consists of three steps:
        1. Estimate the PC values for the individual.
        2. Estimate the population admixture proportion of the individual.
        3. Check if he is of the given ancestry, using only the top two PCs.
        
    :param genotype_file: The genotype file with the individual genotype.
    :param weights_file: The weights_file file with the PC SNP weights.
    :param pcs_file: The file with PCs and admixture information for reference genotypes (e.g. 1000 genomes).
    :param check_population: The population/ancestry to check.
        
    """
    weight_dict, stats_dict = parse_pc_weights(weights_file)
    pcs_admixture_dict = load_pcs_admixture_info(pcs_file)
    pcs = pcs_admixture_dict['pcs']
    pop_dict = pcs_admixture_dict['pop_dict']
    print pop_dict['unique_populations']
    genotype_d = calc_indiv_genot_pcs(genotype_file, weight_dict, len(pop_dict['unique_populations'])-1, **kwargs)
    genotype_pcs = genotype_d['pcs'][0]

    admixture = calc_admixture(genotype_pcs, pop_dict['admix_decom_mat'])
    
    pop_filter = sp.in1d(pop_dict['populations'],[check_population])
    check_population = check_in_population(pcs[pop_filter], genotype_pcs[0],genotype_pcs[1])
    ancestry_dict={}
    ancestry_dict = {'check_population': check_population, 'admixture': admixture}
    return ancestry_dict



def check_in_population(pcs, pc1, pc2):
    """
    Check if PC1 and PC2 are within the populations' PCs
    :param pcs: pcs of the population
    :param pc1: PC1 of the individual
    :param pc2: PC2 of the individual
    :return: Dictionary with various statistics
    """
    # Report ancestry.
    pop_mean = sp.mean(pcs, 0)
    pop_std = sp.std(pcs, 0)
    pop_lim = (3 * pop_std) ** 2
    ind_pcs = sp.array([pc1, pc2])
    ind_lim = (ind_pcs - pop_mean) ** 2
    is_in_population = sp.any(ind_lim ** 2 < pop_lim)
    return {'pop_lim': pop_lim, 'pop_mean': pop_mean, 'pop_std': pop_std, 'ind_lim': ind_lim,
            'is_in_population': is_in_population,'pc1':pc1,'pc2':pc2}


def plot_pcs(plot_file, pcs, populations, indiv_pcs=None):
    """
    Plots the PCs of the hapmap and if provided of the genotype
    :param populations: dictionary with different population masks for coloring the individuals
    :param plot_file: Name of the file the plot should be saved (either .png or .pdf)
    :param pcs: principal components from the hapmap dataset
    :param genotype_pcs_dict: Dictionary with PCs of the individual (optional)
    """
    log.info('Plotting PCs of Hapmap')
    # Plot them
    unique_pops = sp.unique(populations)
    for pop in unique_pops:
        pop_filter = sp.in1d(population, [pop])
        pylab.plot(pcs[pop_filter][0], pcs[pop_filter][1], label=pop, ls='', marker='.', alpha=0.6)

    log.info('Plotting genome on plot')
    # Project genome on to plot.
    if indiv_pcs is not None:
        pylab.plot(indiv_pcs[0], indiv_pcs[1], 'o', label='This is you')
    pylab.xlabel('PC 1')
    pylab.ylabel('PC 2')
    pylab.legend(loc=4, numpoints=1)
    pylab.savefig(plot_file)


def _calc_pcs(weight_dict, sids, nts, snps, num_pcs_to_use):
    num_nt_issues = 0
    num_snps_used = 0
    num_indivs = snps.shape[1]
    pcs = sp.zeros((num_indivs, num_pcs_to_use))

    for snp_i, sid in enumerate(sids):
        try:
            d = weight_dict[sid]
        except:
            continue
        nt = nts[snp_i]
        snp = snps[snp_i]
        pc_weights = sp.array(d['pc_ws'])
        pc_weights = pc_weights[:num_pcs_to_use]
        pc_weights.shape = (1, num_pcs_to_use)
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


def get_snps_filter(nt_map_file):
    nt_map = cPickle.load(open(nt_map_file,'r'))
    snps_filter = {}
    for chrom in range(1, 23):
        chrom_str = 'chr%d' % chrom
        snps_filter[chrom_str] = nt_map[chrom_str]['sids']
    return snps_filter

def calc_admixture(pred_pcs, admix_decomp_mat):    
    """
    Get admixture decomp.  Predicted PCs correspond to the admix_decomp_mat.
    """
    log.info('Decomposing individual admixture')

    v=sp.concatenate((pred_pcs,[1.0]))
    admixture = sp.dot(v, admix_decomp_mat)
    assert sp.sum(admixture)==1, "Admixture doesn't sum to 1: "+str(admixture)
    return admixture


#For debugging purposes
def _test_admixture_():
    pc_weights_file = '/faststorage/project/TheHonestGene/snpweights/snpwt.CO'
    pc_weights_hdf5_file = '/faststorage/project/TheHonestGene/snpweights/snpwt.CO.hdf5'
    nt_map_file = '/faststorage/project/TheHonestGene/data_for_pipeline/NT_DATA/23andme_v4_nt_map.pickled'
    pc_ref_genot_file = '/faststorage/project/TheHonestGene/data_for_pipeline/1k_genomes_hg.hdf5'
    indiv_genot_file = '/faststorage/project/TheHonestGene/prediction_data/23andme-genomes_imputed/3a9c0f27a91816e7.genome_imputed.hdf5'
    ref_pcs_admix_file = '/faststorage/project/TheHonestGene/test_data/1kg_CEPH_pcs_admix_data.hdf5'
    pcs_plot_file = '/faststorage/project/TheHonestGene/test_data/pc_plot.png'
    
    #Parse and save weights file
    print 'Parsing SNP weights from text file'
    sid_weights_map, stats_dict = _parse_pc_weights_from_text(pc_weights_file)
    print 'Storing SNP weights'
    save_pc_weights(sid_weights_map, stats_dict, pc_weights_hdf5_file)
#     sid_weights_map, stats_dict = _parse_pc_weights_from_hdf5(pc_weights_hdf5_file)
    
    #Generate a snps_filter based on an individual genotype??
    print 'Loading SNP filter'
    snps_filter = get_snps_filter(nt_map_file)
    
    #Generate and save PC/admixture info file for 1000 genomes.
    print 'Calculating PC projections and admixture decomposition information'
    pcs_dict = calc_genot_pcs(pc_ref_genot_file, sid_weights_map, stats_dict, populations_to_use = ['EUR','AFR','EAS'], 
                              snps_filter=snps_filter, verbose=True)
    print 'Save projected PCs and admixture decomposition to file'
    save_pcs_admixture_info(pcs_dict['pcs'], pcs_dict['pop_dict'], ref_pcs_admix_file)

    print 'Loading pre-calculated projected PCs and admixture decomposition to file'
    pcs_dict = load_pcs_admixture_info(ref_pcs_admix_file)
    
    print "Plot PC projection for the genotypes (for debugging purposes)."
    plot_pcs(pcs_plot_file, pcs_dict['pcs'], pcs_dict['pop_dict']['populations'])

    # Calculate admixture for an individual
    print 'Calculate admixture for an individual.'
    ancestry_results =  ancestry_analysis(indiv_genot_file, pc_weights_hdf5_file, ref_pcs_admix_file, check_population='EUR',)
    print ancestry_results
    
    #Plot PCs..
    
    