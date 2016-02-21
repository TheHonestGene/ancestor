from os import path
from ancestor.core import ancestry as an
import numpy as np
from numpy.testing import assert_almost_equal

data_path = path.join(path.dirname(__file__), 'data')

stats_to_assert = {'shrinkage': [0.988256270989, 0.975872942881],
                   'linear_transform': [-8.77591739034, -1.89552894559, 0.286808231529, 1.76506268675, 8.77952492895,
                                        0.28426484374, 7.01085470359, -6.88399598336, 0.428926924731],
                   'populations': {'YRI': {'count': 113, 'meanpc': [-0.0776451327434, -0.0167681415929]},
                                   'CEU': {'count': 112, 'meanpc': [0.0157571428571, 0.0783553571429]},
                                   'ASI': {'count': 169, 'meanpc': [0.041475739645, -0.0407165680473]}
                                   }
                   }
                   
ancestry_to_assert = {'pop_std': np.asarray([ 0.00059839,  0.00111145]), 
    'pc2': -0.01174512492642139, 'is_in_population': True,
    'ind_lim': np.asarray([  3.77975763e-08,   4.77972171e-07]),
    'pop_lim': np.asarray([  3.22262193e-06,   1.11179861e-05]),
    'pc1': -0.0070095622547185791,
    'pop_mean':np.asarray([-0.00681515, -0.01105377]),
    'population':'EUR'}



def test_parse_weights_from_hdf5(weights_dict_hdf5):
    assert_weights(weights_dict_hdf5, 712835)


def test_parse_weights_from_text(weights_dict):
    assert_weights(weights_dict, 15)


def test_load_pcs_from_file():
    pcs_dict = an.load_pcs_from_file(path.join(data_path, 'pcs.hdf5'))
    assert type(pcs_dict) == dict
    populations = pcs_dict['populations']
    assert type(populations) == dict
    assert len(populations.keys()) == 5
    assert pcs_dict['pcs'].shape == (2504, 2)


def test_calc_genotype_pcs_from_hdf5_weights(weights_dict_hdf5):
    weights = weights_dict_hdf5[0]
    pcs_dict = an.calc_genotype_pcs(path.join(data_path, 'genotype.hdf5'), weights)
    assert type(pcs_dict) == dict
    assert pcs_dict == {'pc1': -0.0070095622547185791, 'num_snps_used': 237399, 'pc2': -0.01174512492642139}

def test_ancestry_analysis():
    genotype_file = path.join(data_path, 'genotype.hdf5')
    pcs_file = path.join(data_path, 'pcs.hdf5')
    weights_file = path.join(data_path,'weights.hdf5')
    ancestry_dict = an.ancestry_analysis(genotype_file,weights_file,pcs_file)
    assert_ancestry_dict(ancestry_dict)
     


def assert_weights(weights_dict, length):
    assert type(weights_dict) == tuple
    weights, stats = weights_dict
    assert len(weights) == length
    weight = weights[b'rs4040604']
    assert weight == {'mean_g': np.float32(0.58883), 'nts': [b'T', b'G'], 'pc1w': np.float32(1.6099e-07),
                      'pc2w': np.float32(4.9054e-07)}
    assert_stats(stats)

def assert_ancestry_dict(ancestry_dict):
    assert_almost_equal(ancestry_dict['pop_std'],ancestry_to_assert['pop_std'])
    assert_almost_equal(ancestry_dict['ind_lim'],ancestry_to_assert['ind_lim'])
    assert_almost_equal(ancestry_dict['pop_lim'],ancestry_to_assert['pop_lim'])
    assert_almost_equal(ancestry_dict['pop_mean'],ancestry_to_assert['pop_mean'])
    assert ancestry_dict['pc1'] == ancestry_to_assert['pc1']
    assert ancestry_dict['pc2'] == ancestry_to_assert['pc2']
    assert ancestry_dict['is_in_population'] == ancestry_to_assert['is_in_population']
    assert ancestry_dict['population'] == ancestry_to_assert['population']

def assert_stats(stats):
    assert stats['shrinkage'] == stats_to_assert['shrinkage']
    assert stats['linear_transform'] == stats_to_assert['linear_transform']
    for population in ('YRI', 'CEU', 'ASI'):
        assert stats['populations'][population]['meanpc'] == stats_to_assert['populations'][population]['meanpc']
        assert int(stats['populations'][population]['count']) == int(
                stats_to_assert['populations'][population]['count'])
