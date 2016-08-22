from os import path
from ancestor.core import ancestry as an
import numpy as np
from numpy.testing import assert_almost_equal

data_path = path.join(path.dirname(__file__), 'data')

stats_to_assert = {'shrinkage': [0.869045547838, 0.769657143423, 0.731475178282, 0.707642249858,
                                 0.655952609988, 0.646263302717, 0.641821371264, 0.625572399068,
                                 0.61926635542, 0.606826428827],
                   'linear_transform': [-5.76876414694, -8.03949196612, 6.07540134769, 1.2481264795, 2.93971486278, -5.04224636121, 2.65646034124, 5.85445717023, 11.004353167, -3.16204430647, 0.154123326805, -4.23309525956, -0.168369912964, 8.92767144022, 1.1195689241, 7.37548672737, -9.22367965627, 16.4590902363, 7.88148057715, 29.121896897, -4.76478153779, 0.154158760754, -4.76171760857, -6.68961957283, -1.73440542166, 0.406757339145, 9.00788503093, -13.4228861676, 18.1586313146, 1.0944738413, 23.6640819765, -4.34631699487, 0.167239433183, -13.084541036, -50.8225343421, 20.383086267, 4.13062817327, 96.8584940546, -71.6177177904, 69.6092326526, 135.150354996, 209.233789529, -38.8646837848, 0.0563236129631, 2.615736314, -0.608369710473, -14.3164494857, -2.50835820036, 16.5527841997, -22.8793507543, -8.36718001031, -69.0177171871, -4.44324006314, -10.2501158653, 0.0917284057361, 126.321645274, 419.795128973, -261.846693135, -27.1486636822, -767.548209973, 319.103916214, -1065.35820143, -1057.60492954, -1919.27771579, 283.487280286, 0.0904965315739, -56.4286699375, -191.307609814, 118.449971985, 18.1470014918, 313.797621467, -180.259567667, 443.715234025, 343.193618884, 682.977739944, -130.112125791, 0.0539801545111, 5.96870803415, 10.8397159673, -12.8055032755, 0.32843385167, 11.5915662282, 23.7971472681, 16.3115545734, -6.17178374893, 70.8464492181, 4.42254845392, 0.0820991431138, 55.8393796328, 157.859656584, -125.885938292, -10.9917118328, -163.017858994, 189.770767912, -279.57981735, -407.588172501, -667.841236553, 144.823249929, 0.0245848374495, -79.6149570212, -231.769921095, 207.286562856, 10.9247851334, 275.30794351, -123.490535047, 578.172937741, 764.254689503, 1089.55878692, -185.803268294, 0.0453500054135, -26.8537242456, -99.0885851111, 55.4662957138, 4.34343232251, 197.134572886, -106.73584795, 208.222057902, 282.953528002, 475.155094749, -55.4297420951, 0.079915788497],
                   'populations': {'French': {'count': 29, 'meanpc': [0.0559620689655, -0.0142068965517, 0.00764137931034, 0.00938620689655, -0.00535862068966, -0.0044275862069, 0.00198620689655, 0.000972413793103, -0.000648275862069, 0.00444482758621]}, 'Sardinian': {'count': 28, 'meanpc': [0.0333214285714, -0.0108035714286, 0.016375, -0.131114285714, 0.00819285714286, 0.0062, -0.00438214285714, -0.00646428571429, -0.000246428571429, -0.02365]}, 'Tuscan': {'count': 8, 'meanpc': [0.0282375, -0.002975, 0.0022875, -0.016075, 0.006275, 0.0043625, 0.0023625, 0.0036125, -0.0049, -0.0066625]}, 'Bedouin': {'count': 47, 'meanpc': [-0.081414893617, -0.0603680851064, 0.0277191489362, 0.0132319148936, 0.0131042553191, 0.040170212766, -0.0355127659574, -0.0163276595745, 0.0100106382979, 0.0158893617021]}, 'Russian': {'count': 25, 'meanpc': [0.070876, -0.01818, -0.003232, 0.104568, -0.011388, 0.003784, -0.001404, -0.004228, 0.008872, 0.005972]}, 'Druze': {'count': 47, 'meanpc': [-0.0337638297872, 0.0987872340426, 0.0441255319149, 0.00546170212766, 0.0245042553191, -0.0104361702128, 0.00306170212766, 0.00269361702128, -0.00494680851064, -0.0241234042553]}, 'Palestinian': {'count': 51, 'meanpc': [-0.0508058823529, 0.00212549019608, -0.0853882352941, -0.00789215686275, -0.030968627451, -0.0264725490196, 0.0307078431373, 0.0133019607843, -0.00541764705882, 0.0138450980392]}, 'Orcadian': {'count': 16, 'meanpc': [0.07079375, -0.01940625, 0.00970625, 0.04715625, -0.00541875, -0.00460625, 0.0049625, -0.000625, -0.001725, 0.00550625]}, 'Adygei': {'count': 17, 'meanpc': [0.0145705882353, 0.0179352941176, -0.0219941176471, 0.0730294117647, 0.0141705882353, 0.00723529411765, -0.00432352941176, 0.0127411764706, -0.00914705882353, -0.0471117647059]}},
                   'num_pcs': 10}

ancestry_to_assert = {'indiv_pcs': np.asarray([-0.03520987, -0.00157889, -0.02793205]),
                      'check_population': {'ref_pop_std':np.asarray([ 0.00125915,  0.00123291,  0.00124285]),
                                           'ref_pop_mean_pcs':np.asarray([-0.03426199, -0.00177332, -0.02924874]),
                                           'is_in_population':True,
                                           'check_pop':b'GBR'},
                      'admixture': {'unadjusted_admixture': np.asarray([ 0.09676908,  0.57078818, -0.02539509,  0.35783837]),
                                    'admixture':np.asarray([ 0.09437244,  0.55665166,  0.        ,  0.34897591]),
                                    'confidence': 'Very good',
                                    'confidence_score':0.0002322598308503443}}




def test_parse_weights_from_hdf5(weights_dict_hdf5):
    assert_weights(weights_dict_hdf5, 574874)


def test_parse_weights_from_text(weights_dict):
    assert_weights(weights_dict, 15)


def test_load_pcs_from_file():
    pcs_dict = an.load_pcs_admixture_info(path.join(data_path, 'pcs.hdf5'))
    assert isinstance(pcs_dict, dict)
    pop_dict = pcs_dict['pop_dict']
    assert isinstance(pop_dict, dict)
    assert pop_dict['unique_populations'].tolist() == [b'TSI', b'FIN', b'IBS', b'GBR']
    assert_almost_equal(pop_dict['avg_pcs'],np.array([[-0.01927088, -0.00627075, -0.02725555],
       [-0.03808748, -0.00068919, -0.02726736],
       [-0.02579543, -0.0021997 , -0.02896758],
       [-0.03426199, -0.00177332, -0.02924874]]))
    assert pop_dict['num_indivs'].tolist() == [107,  99, 107,  91]
    assert_almost_equal(pop_dict['admix_decom_mat'], np.array([[ -14.72152901,   -5.30980682,  141.77890015, -121.74756622],
       [-228.58602905,  162.11912537,  477.99157715, -411.52468872],
       [  96.64972687,  405.74151611,   12.19804382, -514.58929443],
       [   1.91713834,   11.97299385,    6.06203556,  -18.95216751]], dtype=np.float32))
    assert pop_dict['populations'].shape == (404,)


def test_calc_genotype_pcs_from_hdf5_weights(weights_dict_hdf5):
    weights = weights_dict_hdf5[0]
    num_of_pcs = 3
    pcs_dict = an.calc_indiv_genot_pcs(path.join(data_path, 'genotype.hdf5'), weights,num_of_pcs)
    assert isinstance(pcs_dict, dict)
    assert pcs_dict['num_snps_used'] == 318808
    assert_almost_equal(pcs_dict['pcs'], np.asarray([[-0.03520987, -0.00157889, -0.02793205]]))

def test_ancestry_analysis():
    genotype_file = path.join(data_path, 'genotype.hdf5')
    pcs_file = path.join(data_path, 'pcs.hdf5')
    weights_file = path.join(data_path, 'weights.hdf5')
    ancestry_dict = an.ancestry_analysis(genotype_file, weights_file, pcs_file,check_population=b'GBR')
    assert_ancestry_dict(ancestry_dict)



def assert_weights(weights_dict, length):
    assert isinstance(weights_dict, tuple)
    weights, stats = weights_dict
    assert len(weights) == length
    weight = weights[b'rs11260588']

    assert weight['mean_g'] == np.float32(0.99835998)
    assert weight['nts'] == [b'G', b'A']
    assert_almost_equal(weight['pc_ws'],
                        np.array([5.48959974e-07, 1.05260006e-06, -1.21640005e-06,
                                  -2.98470013e-07, 1.83859996e-07, 5.10689972e-07,
                                  -1.60859997e-07, 5.83359991e-08, 2.54170001e-08,
                                  1.73710006e-07], dtype=np.float32))
    assert_stats(stats)

def assert_ancestry_dict(ancestry_dict):
    pop_check = ancestry_dict['check_population']
    pop_chec_assert = ancestry_to_assert['check_population']
    admixture = ancestry_dict['admixture']
    admixture_assert = ancestry_to_assert['admixture']

    assert_almost_equal(ancestry_dict['indiv_pcs'],ancestry_to_assert['indiv_pcs'])
    assert_almost_equal(pop_check['ref_pop_std'], pop_chec_assert['ref_pop_std'])
    assert_almost_equal(pop_check['ref_pop_mean_pcs'], pop_chec_assert['ref_pop_mean_pcs'])
    assert pop_check['is_in_population'] == pop_chec_assert['is_in_population']
    assert pop_check['check_pop'] == pop_chec_assert['check_pop']
    assert_almost_equal(admixture['unadjusted_admixture'],admixture_assert['unadjusted_admixture'])
    assert_almost_equal(admixture['admixture'],admixture_assert['admixture'])
    assert admixture['confidence'] == admixture_assert['confidence']
    assert_almost_equal(admixture['confidence_score'], admixture_assert['confidence_score'])



def assert_stats(stats):
    assert stats['shrinkage'] == stats_to_assert['shrinkage']
    assert stats['linear_transform'] == stats_to_assert['linear_transform']
    for population in ('French', 'Sardinian', 'Tuscan', 'Bedouin', 'Russian', 'Palestinian', 'Orcadian', 'Druze'):
        assert stats['populations'][population]['meanpc'] == stats_to_assert['populations'][population]['meanpc']
        assert int(stats['populations'][population]['count']) == int(stats_to_assert['populations'][population]['count'])
