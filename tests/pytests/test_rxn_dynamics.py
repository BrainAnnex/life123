import pytest
import numpy as np

from life123 import RxnDynamics



def test_norm_A():
    delta_conc = np.array([1, 4])
    result = RxnDynamics.norm_A(delta_conc)
    assert np.allclose(result, 4.25)

    delta_conc = np.array([.5, 2])
    result = RxnDynamics.norm_A(delta_conc)
    assert np.allclose(result, 1.0625)

    delta_conc = np.array([.5, 2, 1])
    result = RxnDynamics.norm_A(delta_conc)
    assert np.allclose(result, 0.5833333333)

    delta_conc = np.array([0.5, 1, 4, -2, -5])
    result = RxnDynamics.norm_A(delta_conc)
    assert np.allclose(result, 1.85)



def test_norm_B():
    delta = np.array([4, -1])
    base = np.array([10, 2])
    result = RxnDynamics.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.5)

    base = np.array([10, 0])
    result = RxnDynamics.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.4)     # The zero baseline concentration is disregarded

    delta = np.array([300,        6, 0, -1])
    base = np.array([0.00000001, 10, 0,  2])
    result = RxnDynamics.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.6)

    delta = np.array([300,   6, 0, -1])
    base = np.array([0.001, 10, 0,  2])
    result = RxnDynamics.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 300000)

    base = np.array([2,   5, 5, 14, 14])
    delta = np.array([0.5, 1, 4, -2, -5])
    result = RxnDynamics.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.8)

    with pytest.raises(Exception):
        RxnDynamics.norm_B(baseline_conc=np.array([1, 2, 3]), delta_conc=delta) # Too many entries in array

    with pytest.raises(Exception):
        RxnDynamics.norm_B(baseline_conc=base, delta_conc=np.array([1]))        # Too few entries in array



def test_norm_C():
    result = RxnDynamics.norm_C(prev_conc=np.array([1]), baseline_conc=np.array([2]), delta_conc=np.array([0.5]))
    assert result == 0

    result = RxnDynamics.norm_C(prev_conc=np.array([8]), baseline_conc=np.array([5]), delta_conc=np.array([1]))
    assert result == 0

    result = RxnDynamics.norm_C(prev_conc=np.array([8]), baseline_conc=np.array([5]), delta_conc=np.array([4]))
    assert np.isclose(result, 4/3)

    result = RxnDynamics.norm_C(prev_conc=np.array([10]), baseline_conc=np.array([14]), delta_conc=np.array([-2]))
    assert result == 0

    result = RxnDynamics.norm_C(prev_conc=np.array([10]), baseline_conc=np.array([14]), delta_conc=np.array([-5]))
    assert np.isclose(result, 5/4)

    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])
    result = RxnDynamics.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3)

    # A scenario where the 'prev' and 'baseline' values are almost identical
    prev = np.append(prev, 3)
    baseline = np.append(baseline, 2.999999999)
    delta = np.append(delta, 8)
    result = RxnDynamics.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3)

    # A scenario where the 'delta' dwarfs the change between 'prev' and 'baseline'
    prev = np.append(prev, 10)
    baseline = np.append(baseline, 10.05)
    delta = np.append(delta, -9)
    result = RxnDynamics.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3)

    prev = np.append(prev, 10)
    baseline = np.append(baseline, 10.2)
    delta = np.append(delta, -9)
    result = RxnDynamics.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 45)



def test_norm_D():
    prev =     np.array([ 12.96672432,  31.10067726,  55.93259842,  44.72389482, 955.27610518])
    baseline = np.array([ 12.99244738,  31.04428765,  55.96326497,  43.91117372, 956.08882628])
    delta =    np.array([-2.56160549,   -0.03542113,   2.59702662,   1.36089356,  -1.36089356])
    result = RxnDynamics.norm_D(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 37.64942285399873)
