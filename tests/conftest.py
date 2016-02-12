import pytest
from os import path
from ancestor.core import ancestry as an

data_path = path.join(path.dirname(__file__),'data')

@pytest.fixture(scope="module")
def weights_dict():
    return an.parse_pc_weights(path.join(data_path,'weights_subset.txt'))

@pytest.fixture(scope="module")
def weights_dict_hdf5():
    return an.parse_pc_weights(path.join(data_path,'weights.hdf5'))