import pytest
import numpy as np
from life123.units import K, C, F, KJ_PER_MOL, J_PER_MOL, standardize_units, convert



def test_standardize_units():

    result = standardize_units(value=3)
    assert result == 3

    result = standardize_units(value=(3, K))
    assert np.allclose(result, 3)

    result = standardize_units(value=(0, C))
    assert np.allclose(result, 273.15)

    result = standardize_units(value=(32, F))
    assert np.allclose(result, 273.15)


    result = standardize_units(value=(5, KJ_PER_MOL))
    assert np.allclose(result, 5)

    result = standardize_units(value=(13000, J_PER_MOL))
    assert np.allclose(result, 13)



def test_convert():

    result = convert(quantity=212, from_unit=F, to_unit=C)
    assert np.allclose(result, 100)

    result = convert(quantity=37, from_unit=C, to_unit=C)
    assert np.allclose(result, 37)

    result = convert(quantity=13000, from_unit=J_PER_MOL, to_unit=KJ_PER_MOL)
    assert np.allclose(result, 13)

    with pytest.raises(Exception):
        convert(quantity=37, from_unit=C, to_unit=J_PER_MOL)    # Nonsensical
