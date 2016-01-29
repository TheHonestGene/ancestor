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
            mean_g = l[3]
            pcw1 = l[4]
            pcw2 = l[5]
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
        snps = g_cg['sids'][...]
        
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



def plot_1KG_PCs(kg_file='...'):
    #Load genotypes
    #Project on the PCs
    #Plot them
    #Store coordinates
    pass

def plot_genome_pcs():
    #Load 1K genomes coordinates and plot.
    #Project genome on to plot.
    #Report ancestry. 
    pass

    
