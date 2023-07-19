import pytest
import numpy as np
from src.modules.reactions.reaction import ThermoDynamics



def test_K_from_delta_G():
    assert np.allclose(1., ThermoDynamics.K_from_delta_G(delta_G = 0, temp = 100))
    assert np.allclose(1., ThermoDynamics.K_from_delta_G(delta_G = 0, temp = 2000))

    assert np.allclose(1.00120344711, ThermoDynamics.K_from_delta_G(delta_G = -1, temp = 100))
    assert np.allclose(0.99879799943, ThermoDynamics.K_from_delta_G(delta_G =  1, temp = 100))

    assert np.allclose(3.32917175365, ThermoDynamics.K_from_delta_G(delta_G = -1000, temp = 100))
    assert np.allclose(0.30037501036, ThermoDynamics.K_from_delta_G(delta_G =  1000, temp = 100))

    assert np.allclose(11.0833845653, ThermoDynamics.K_from_delta_G(delta_G = -1000, temp = 50))
    assert np.allclose(0.09022514685, ThermoDynamics.K_from_delta_G(delta_G =  1000, temp = 50))

    assert np.allclose(1.27194180118, ThermoDynamics.K_from_delta_G(delta_G = -1000, temp = 500))
    assert np.allclose(0.78619949361, ThermoDynamics.K_from_delta_G(delta_G =  1000, temp = 500))

    delta_G = 8000.
    K = ThermoDynamics.K_from_delta_G(delta_G = delta_G, temp = 300)
    assert np.allclose(delta_G , ThermoDynamics.delta_G_from_K(K = K, temp = 300))   # Reverse ops



def test_delta_G_from_K():
    assert np.allclose(0, ThermoDynamics.delta_G_from_K(K = 1, temp = 100))
    assert np.allclose(0, ThermoDynamics.delta_G_from_K(K = 1, temp = 50))
    assert np.allclose(0, ThermoDynamics.delta_G_from_K(K = 1, temp = 1000))

    assert np.allclose(-831.446261815324, ThermoDynamics.delta_G_from_K(K = 2.71828182846, temp = 100))
    assert np.allclose(831.446261815324, ThermoDynamics.delta_G_from_K(K = 1/2.71828182846, temp = 100))

    assert np.allclose(-8314.46261815324, ThermoDynamics.delta_G_from_K(K = 2.71828182846, temp = 1000))
    assert np.allclose(8314.46261815324, ThermoDynamics.delta_G_from_K(K = 1/2.71828182846, temp = 1000))

    assert np.allclose(-83.1446261815324, ThermoDynamics.delta_G_from_K(K = 2.71828182846, temp = 10))
    assert np.allclose(83.1446261815324, ThermoDynamics.delta_G_from_K(K = 1/2.71828182846, temp = 10))

    K = 4.2
    delta_G = ThermoDynamics.delta_G_from_K(K = K, temp = 273.)
    assert np.allclose(K , ThermoDynamics.K_from_delta_G(delta_G = delta_G, temp = 273.))   # Reverse ops



def test_delta_G_from_enthalpy():
    # When there's no Entropy change, the Delta_G equals the change in Enthalpy
    assert np.allclose(800. , ThermoDynamics.delta_G_from_enthalpy(delta_H=800., delta_S=0., temp=100))
    assert np.allclose(800. , ThermoDynamics.delta_G_from_enthalpy(delta_H=800., delta_S=0., temp=200))

    # Likewise, at absolute zero temperature, the Delta_G would equal the change in Enthalpy
    assert np.allclose(800. , ThermoDynamics.delta_G_from_enthalpy(delta_H=800., delta_S=-15., temp=0))
    assert np.allclose(800. , ThermoDynamics.delta_G_from_enthalpy(delta_H=800., delta_S=50.,  temp=0))

    # But as soon as the temperature deviates from absolute zero, they're no longer the same
    assert not np.allclose(800. , ThermoDynamics.delta_G_from_enthalpy(delta_H=800., delta_S=50.,  temp=1))

    # Progressively larger changes in Entropy will make the reaction increasingly favored energetically
    assert np.allclose(4000. , ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=-20.,temp=100))
    assert np.allclose(1900. , ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=1.,  temp=100))
    assert np.allclose(   0. , ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=20., temp=100))
    assert np.allclose(-2000. ,ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=40., temp=100))

    # When the Entropy change is positive, higher temps will make the reaction more favored energetically...
    assert np.allclose(1000. , ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=10., temp=100))
    assert np.allclose(0. ,    ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=10., temp=200))
    assert np.allclose(-1000. ,ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=10., temp=300))

    # but when the Entropy change is negative, higher temps will make the reaction LESS favored energetically
    assert np.allclose(3000. , ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=-10., temp=100))
    assert np.allclose(4000. , ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=-10., temp=200))
    assert np.allclose(5000. , ThermoDynamics.delta_G_from_enthalpy(delta_H=2000., delta_S=-10., temp=300))



def test_delta_H_from_gibbs():
    # When there's no Entropy change, the Delta_H equals the change in Gibbs free energy
    assert np.allclose(800. , ThermoDynamics.delta_H_from_gibbs(delta_G=800., delta_S=0., temp=100))
    assert np.allclose(800. , ThermoDynamics.delta_H_from_gibbs(delta_G=800., delta_S=0., temp=200))

    # Likewise, at absolute zero temperature, the Delta_H would equal to the change in Gibbs free energy
    assert np.allclose(800. , ThermoDynamics.delta_H_from_gibbs(delta_G=800., delta_S=-15., temp=0))
    assert np.allclose(800. , ThermoDynamics.delta_H_from_gibbs(delta_G=800., delta_S=50.,  temp=0))

    # But as soon as the temperature deviates from absolute zero, they're no longer the same
    assert not np.allclose(800. , ThermoDynamics.delta_H_from_gibbs(delta_G=800., delta_S=50.,  temp=1))

    delta_H = ThermoDynamics.delta_H_from_gibbs(delta_G=3000., delta_S=20.,temp=100)
    assert np.allclose(delta_H, 5000.)
    # Get back the original delta_G
    assert np.allclose(3000. , ThermoDynamics.delta_G_from_enthalpy(delta_H=delta_H, delta_S=20., temp=100))



def test_delta_S_from_gibbs():
    with pytest.raises(Exception):
        # Division by zero
        ThermoDynamics.delta_S_from_gibbs(delta_G=3000, delta_H=3000, temp=0)

    # If delta_G and delta_H are the same, then there's no Entropy change
    assert np.allclose(0. , ThermoDynamics.delta_S_from_gibbs(delta_G=2000, delta_H=2000, temp=100))

    delta_S = ThermoDynamics.delta_S_from_gibbs(delta_G=3000, delta_H=5000, temp=100)
    assert np.allclose(delta_S, 20.)
    # Get back the original delta_G
    assert np.allclose(3000. , ThermoDynamics.delta_G_from_enthalpy(delta_H=5000, delta_S=delta_S, temp=100))
