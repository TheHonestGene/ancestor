"""

1. Parse PC weights for the SNPs in question (or all SNPs)

2. Predict ancestry for 1K genomes individuals
     - Plot individuals.
     - Store plot coordinates.
     - Identify ancestry clusters and boundaries.
     - Store ancestry boundaries.
     
3. Predict ancestry for target individuals  

4.

"""
from itertools import izip
import scipy as sp
import h5py
import pylab

cloud_dir = '/Users/bjv/Dropbox/Cloud_folder/'
repos_dir = '/Users/bjv/REPOS/'
Kg_nt_decoder = {1:'A', 2:'T', 3:'C', 4:'G', }


def parse_PC_weights(pc_weights_file='/Users/bjv/REPOS/ancestor/tests/snpwt.CO'):
    """
Parses a 0.988256270989 0.975872942881 
YRI CEU ASI 
113 112 169 
YRI -0.0776451327434 -0.0167681415929 CEU 0.0157571428571 0.0783553571429 ASI 0.041475739645 -0.0407165680473 
-8.77591739034 -1.89552894559 0.286808231529 1.76506268675 8.77952492895 0.28426484374 7.01085470359 -6.88399598336 0.428926924731 
rs11240767 C T 0.94797 1.9086e-07 8.4583e-08
....

    """
    sid_dict = {}
    with open(pc_weights_file) as f:
        print f.next()
        print f.next()
        print f.next()
        print f.next()
        print f.next()
        for line in f:
            l = line.split()
            sid = l[0]
            nts = [l[1],l[2]]
            mean_g = float(l[3])
            pc1w = float(l[4])
            pc2w = float(l[5])
            sid_dict[sid]={'pc1w':pc1w, 'pc2w':pc2w, 'mean_g':mean_g, 'nts':nts}
            
    return sid_dict



def calc_pc_vals(genotype_file, weight_dict=None, weight_file=None):
    """
    """
    gh5f = h5py.File(genotype_file)
    pc1 = 0.0
    pc2 = 0.0
    
    num_snps_used = 0

    for chrom in range(1,23):
        print 'Working on Chromosome %d'%chrom

        chrom_str = 'Chr%d'%chrom            
        g_cg = gh5f[chrom_str]
        
        #Loading data 
        sids = g_cg['sids'][...]
        nts = g_cg['nts'][...]
        snps = g_cg['snps'][...]
        print sids
        
        num_snps_found = 0
        num_nt_issues = 0
        for sid, nt, snp in izip(sids,nts,snps):
            try:
                d = weight_dict[sid]
            except:
                
                continue
            
            num_snps_found +=1
            if sp.all([nt[1],nt[0]]==d['nts']):
                #print 'Flip sign'
                if snp==0:
                    snp=2
                elif snp==2:
                    snp=0                
            elif not sp.all(nt == d['nts']):
                num_nt_issues+=1
                continue
            pc1w = d['pc1w']
            pc2w = d['pc2w']
            mean_g = d['mean_g']
            af = mean_g/2.0
            sd_g = sp.sqrt(af*(1-af))
            norm_snp = (snp-mean_g)/sd_g
            pc1+=norm_snp*pc1w
            pc2+=norm_snp*pc2w
            num_snps_used+=1
            
        print 'Number of SNP weights found: %d'%num_snps_found
        print 'Number of SNP weight discarded for nucleotide issues: %d'%num_nt_issues
        
    gh5f.close()

    print 'Number of SNPs uesd for PC projection: %d'%num_snps_used
    return {'pc1':pc1, 'pc2':pc2}



def plot_1KG_PCs(kg_file=cloud_dir+'Data/1Kgenomes/1K_genomes_v3.hdf5', plot_file = cloud_dir+'tmp/PC_plot.png', 
                 out_file=repos_dir+'ancestor/tests/data/1kg_pc_coord.hdf5'):
    
    #Load weights to identify which SNPs to use.
    sid_dict = parse_PC_weights()
    ok_sids = sp.array(sid_dict.keys())
    print 'Loaded PC weight for %d SNPs'%(len(ok_sids))
    
    #Load genotypes
    h5f = h5py.File(kg_file)
    eur_filter = h5f['indivs']['continent'][...]=='EUR'
    amr_filter = h5f['indivs']['continent'][...]=='AMR'
    asn_filter = h5f['indivs']['continent'][...]=='ASN'
    afr_filter = h5f['indivs']['continent'][...]=='AFR'
    num_indivs = len(h5f['indivs']['continent'][...])
    print 'Found genotypes for %d individuals'%num_indivs
    
    pcs = sp.zeros((num_indivs,2))
    
    num_nt_issues = 0
    num_snps_used = 0
    for chrom in range(1,3):
        print 'Working on Chromosome %d'%chrom
        chrom_str = 'chr%d'%chrom
        print 'Identifying overlap'
        sids = h5f[chrom_str]['snp_ids'][...]
        ok_snp_filter = sp.in1d(sids, ok_sids)
#         assert sids[ok_snp_filter]==ok_sids, 'WTF?'
        sids = sids[ok_snp_filter]
        
        print 'Loading SNPs'
        snps = h5f[chrom_str]['raw_snps'][...]
        snps = snps[ok_snp_filter]
        

        nts = h5f[chrom_str]['nts'][...]
        nts = nts[ok_snp_filter]
        
        print 'updating PCs'
        for snp_i, sid in enumerate(sids):
            d = sid_dict[sid]
            nt = nts[snp_i]
            nt = [Kg_nt_decoder[nt[0]],Kg_nt_decoder[nt[1]]]
            snp = snps[snp_i]
            pc_weights = sp.array([d['pc1w'],d['pc2w']])
            pc_weights.shape = (1,2)
            af = d['mean_g']/2.0
            if sp.all([nt[1],nt[0]]==d['nts']):
                #print 'Flip sign'
                pc_weights = -pc_weights
                af = 1-af
            elif not sp.all(nt == d['nts']):
                num_nt_issues+=1
                continue
            
            mean_g = 2*af
            sd_g = sp.sqrt(af*(1-af))

            #"Normalizing" the SNPs with the given allele frequencies
            norm_snp = (snp-mean_g)/sd_g
            norm_snp.shape = (num_indivs,1)

            #Project on the PCs
            pcs += sp.dot(norm_snp,pc_weights)
            num_snps_used+=1

    h5f.close()
    print pcs[eur_filter]
    print pcs[asn_filter]
    print pcs[afr_filter]
    print pcs[amr_filter]
    print '%d SNPs were excluded from the analysis due to nucleotide issues.'%(num_nt_issues)
    print '%d SNPs were used for the analysis.'%(num_snps_used)
            
        
    #Plot them
    pylab.plot(pcs[eur_filter][:,0],pcs[eur_filter][:,1], label='EUR', ls='', marker='.',alpha=0.6)
    pylab.plot(pcs[asn_filter][:,0],pcs[asn_filter][:,1], label='ASN', ls='',  marker='.',alpha=0.6)
    pylab.plot(pcs[afr_filter][:,0],pcs[afr_filter][:,1], label='AFR', ls='', marker='.',alpha=0.6)
    pylab.plot(pcs[amr_filter][:,0],pcs[amr_filter][:,1], label='AMR', ls='', marker='.',alpha=0.6)
    pylab.xlabel('PC 1')
    pylab.ylabel('PC 2')
    pylab.legend(loc=4,numpoints=1)
    pylab.savefig(plot_file)
        
    #Store coordinates
    oh5f = h5py.File(out_file)
    oh5f.create_dataset('pcs',data=pcs)
    oh5f.create_dataset('asn_filter',data=asn_filter)
    oh5f.create_dataset('afr_filter',data=afr_filter)
    oh5f.create_dataset('amr_filter',data=amr_filter)
    oh5f.create_dataset('eur_filter',data=eur_filter)
    oh5f.close()
    
    
    
def plot_genome_pcs(genotype_file=repos_dir+'imputor/tests/data/test_out_genotype_imputed.hdf5', 
                    kgenomes_pc_coord_file = repos_dir+'ancestor/tests/data/1kg_pc_coord.hdf5', 
                    plot_file=cloud_dir+'tmp/PC_plot.png'):
    #Load 1K genomes coordinates and plot.
    print 'Loading 1000 genomes coordinates'
    ch5f = h5py.File(kgenomes_pc_coord_file)

    eur_filter = ch5f['eur_filter'][...]
    afr_filter = ch5f['afr_filter'][...]
    amr_filter = ch5f['amr_filter'][...]
    asn_filter = ch5f['asn_filter'][...]
    pcs = ch5f['pcs'][...]
    ch5f.close()
    
    pylab.plot(pcs[eur_filter][:,0],pcs[eur_filter][:,1], label='EUR', ls='', marker='.',alpha=0.6)
    pylab.plot(pcs[asn_filter][:,0],pcs[asn_filter][:,1], label='ASN', ls='', marker='.',alpha=0.6)
    pylab.plot(pcs[afr_filter][:,0],pcs[afr_filter][:,1], label='AFR', ls='', marker='.',alpha=0.6)
    pylab.plot(pcs[amr_filter][:,0],pcs[amr_filter][:,1], label='AMR', ls='', marker='.',alpha=0.6)
    
    #Project genome on to plot.
    sid_dict = parse_PC_weights()    
    ret_dict = calc_pc_vals(genotype_file,sid_dict)
    pylab.plot(ret_dict['pc1'], ret_dict['pc2'], 'o', label='This is you')
    
    pylab.xlabel('PC 1')
    pylab.ylabel('PC 2')
    pylab.legend(loc=4,numpoints=1)
    pylab.savefig(plot_file)
    
    #Report ancestry. 
    

    
