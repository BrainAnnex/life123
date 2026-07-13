import pytest
import numpy as np
from life123 import ThermoDynamics
from life123.units import K, C, J_PER_MOL



def test_equilibrium_constant_from_gibbs_energy():
    assert np.allclose(1., ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = 0, temp = 100))
    assert np.allclose(1., ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = 0, temp = 2000))

    assert np.allclose(1.00120344711, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = -0.001, temp = 100))
    assert np.allclose(0.99879799943, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G =  0.001, temp = 100))

    assert np.allclose(3.32917175365, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = -1, temp = 100))
    assert np.allclose(0.30037501036, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G =  1, temp = 100))

    assert np.allclose(11.0833845653, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = -1, temp = 50))
    assert np.allclose(0.09022514685, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G =  1, temp = 50))

    assert np.allclose(1.27194180118, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = -1, temp = 500))
    assert np.allclose(0.78619949361, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G =  1, temp = 500))

    assert np.allclose(0.78619949361, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G =  1, temp = (226.85, C)))

    delta_G = 8.
    K = ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = delta_G, temp = 300)
    assert np.allclose(delta_G, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = K, temp = 300))   # Reverse ops



def test_gibbs_energy_from_equilibrium_constant():
    assert np.allclose(0, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 1, temp = 100))
    assert np.allclose(0, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 1, temp = 50))
    assert np.allclose(0, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 1, temp = 1000))

    assert np.allclose(-0.831446261815324, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 2.71828182846, temp = 100))
    assert np.allclose(0.831446261815324, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 1 / 2.71828182846, temp = 100))

    assert np.allclose(-8.31446261815324, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 2.71828182846, temp = 1000))
    assert np.allclose(8.31446261815324, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 1 / 2.71828182846, temp = 1000))

    assert np.allclose(-0.0831446261815324, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 2.71828182846, temp = 10))
    assert np.allclose(0.0831446261815324, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 1 / 2.71828182846, temp = 10))

    assert np.allclose(-0.0831446261815324, ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = 2.71828182846, temp = (-263.15, C)))

    K = 4.2
    delta_G = ThermoDynamics.gibbs_energy_from_equilibrium_constant(K = K, temp = 273.)
    assert np.allclose(K, ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = delta_G, temp = 273.))   # Reverse ops



def test_gibbs_energy_from_enthalpy_entropy():
    # When there's no Entropy change, the Delta_G equals the change in Enthalpy
    assert np.allclose(800., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=800., delta_S=0., temp=100))
    assert np.allclose(800., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=800., delta_S=0., temp=200))

    # Likewise, at absolute zero temperature, the Delta_G would equal the change in Enthalpy
    assert np.allclose(800., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=800., delta_S=-15., temp=0))
    assert np.allclose(800., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=800., delta_S=50., temp=0))

    # But as soon as the temperature deviates from absolute zero, they're no longer the same
    result = ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=800., delta_S=50., temp=10)
    assert not np.allclose(800, result)
    assert np.allclose(799.5, result)

    # Progressively larger changes in Entropy will make the reaction increasingly favored energetically
    assert np.allclose(4000., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000, delta_S=-20000, temp=100))
    assert np.allclose(1900., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000, delta_S=1000, temp=100))
    assert np.allclose(0., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000, delta_S=20000, temp=100))
    assert np.allclose(-2000., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000, delta_S=40000, temp=100))

    # When the Entropy change is positive, higher temps will make the reaction more favored energetically...
    assert np.allclose(1000., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000., delta_S=10000, temp=100))
    assert np.allclose(0., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000., delta_S=10000, temp=200))
    assert np.allclose(-1000., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000., delta_S=10000, temp=300))

    # but when the Entropy change is negative, higher temps will make the reaction LESS favored energetically
    assert np.allclose(3000., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000., delta_S=-10000, temp=100))
    assert np.allclose(4000., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000., delta_S=-10000, temp=200))
    assert np.allclose(5000., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=2000., delta_S=-10000, temp=300))



def test_enthalpy_from_gibbs_energy():
    # When there's no Entropy change, the Delta_H equals the change in Gibbs free energy
    assert np.allclose(800., ThermoDynamics.enthalpy_from_gibbs_energy(delta_G=800., delta_S=0., temp=100))
    assert np.allclose(800., ThermoDynamics.enthalpy_from_gibbs_energy(delta_G=800., delta_S=0., temp=200))

    # Likewise, at absolute zero temperature, the Delta_H would equal to the change in Gibbs free energy
    assert np.allclose(800., ThermoDynamics.enthalpy_from_gibbs_energy(delta_G=800., delta_S=-15., temp=0))
    assert np.allclose(800., ThermoDynamics.enthalpy_from_gibbs_energy(delta_G=800., delta_S=50., temp=0))

    # But as soon as the temperature deviates from absolute zero, they're no longer the same
    assert not np.allclose(800., ThermoDynamics.enthalpy_from_gibbs_energy(delta_G=800., delta_S=50., temp=1))

    delta_H = ThermoDynamics.enthalpy_from_gibbs_energy(delta_G=3000, delta_S=20000, temp=100)
    assert np.allclose(delta_H, 5000)
    # Get back the original delta_G
    assert np.allclose(3000., ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=delta_H, delta_S=20000, temp=100))



def test_entropy_from_gibbs_energy():
    with pytest.raises(Exception):
        # Division by zero
        ThermoDynamics.entropy_from_gibbs_energy(delta_G=3, delta_H=3, temp=0)

    # If delta_G and delta_H are the same, then there's no Entropy change
    assert np.allclose(0., ThermoDynamics.entropy_from_gibbs_energy(delta_G=2, delta_H=2, temp=100))

    delta_S = ThermoDynamics.entropy_from_gibbs_energy(delta_G=3, delta_H=5, temp=100)
    assert np.allclose(delta_S, 20)
    # Get back the original delta_G
    assert np.allclose(3, ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=5, delta_S=delta_S, temp=100))



def test_extract_thermodynamic_data():
    with pytest.raises(Exception):
        ThermoDynamics.extract_thermodynamic_data(delta_H=1, delta_S=2, temp=None)  # Missing temp


    result = ThermoDynamics.extract_thermodynamic_data(temp=100)
    assert result == {"K": None, "delta_H": None, "delta_S": None, "delta_G": None}


    result = ThermoDynamics.extract_thermodynamic_data(temp=100, K=3)
    assert result["K"] == 3
    assert result["delta_H"] is None
    assert result["delta_S"] is None
    assert np.allclose(result["delta_G"], -0.913437080597)


    result = ThermoDynamics.extract_thermodynamic_data(temp=100, delta_G=0.8)
    assert np.allclose(result["K"], 0.38205953171)
    assert result["delta_H"] is None
    assert result["delta_S"] is None
    assert np.allclose(result["delta_G"], 0.8)


    with pytest.raises(Exception):
        # Inconsistent thermodynamic data
        ThermoDynamics.extract_thermodynamic_data(temp=100, K=3, delta_G=666.)


    result = ThermoDynamics.extract_thermodynamic_data(temp=300, delta_H=10, delta_S=-40)
    assert np.allclose(result["K"], 0.000147752393)
    assert result["delta_H"] == 10
    assert result["delta_S"] == -40
    assert result["delta_G"] == 22


    with pytest.raises(Exception):
        # Inconsistent thermodynamic data
        ThermoDynamics.extract_thermodynamic_data(temp=100, delta_H=500, delta_S=-3, delta_G=999)


    result = ThermoDynamics.extract_thermodynamic_data(temp=100, delta_H=500, K=0.38205953171)
    assert np.allclose(result["K"], 0.38205953171)
    assert result["delta_H"] == 500
    assert np.allclose(result["delta_S"], 4992)
    assert np.allclose(result["delta_G"], 0.8)


    with pytest.raises(Exception):
        # Inconsistent thermodynamic data
        ThermoDynamics.extract_thermodynamic_data(temp=100, K=2, delta_H=500, delta_S=-3)


    result = ThermoDynamics.extract_thermodynamic_data(temp=300, delta_H=30, delta_G=20)
    assert np.allclose(result["K"], 0.00032942659596407563)
    assert result["delta_H"] == 30
    assert np.allclose(result["delta_S"], 33.3333333)
    assert result["delta_G"] == 20


    result = ThermoDynamics.extract_thermodynamic_data(temp=350, delta_S=-50, delta_G=10)
    assert np.allclose(result["K"], 0.03218183869535201)
    assert np.allclose(result["delta_H"], -7.5)
    assert result["delta_S"] == -50
    assert result["delta_G"] == 10



def test_relative_population_states():
    # Example data about the conformations of methylcyclohexane molecules
    # State 1 : axial ; state 2 : equatorial
    # From Atkins Physical Chemistry (12th edn, 2023), p. xxxvii
    result = ThermoDynamics.relative_population_states(delta_molar_energy=-6, temp=(26.85, C))
    assert np.allclose(result, 0.090225146850218)

    result = ThermoDynamics.relative_population_states(delta_molar_energy=-6, temp=(-270, C))
    assert np.allclose(result, 0)

    result = ThermoDynamics.relative_population_states(delta_molar_energy=-50, temp=(25, C))
    assert np.allclose(result, 0)

    result = ThermoDynamics.relative_population_states(delta_molar_energy=10, temp=(25, C))
    assert np.allclose(result, 56.48383858609265)

    result = ThermoDynamics.relative_population_states(delta_molar_energy=10, temp=(298.15, K))
    assert np.allclose(result, 56.48383858609265)

    result = ThermoDynamics.relative_population_states(delta_molar_energy=10, temp=298.15)
    assert np.allclose(result, 56.48383858609265)

    result = ThermoDynamics.relative_population_states(delta_molar_energy=(10000, J_PER_MOL), temp=298.15)
    assert np.allclose(result, 56.48383858609265)

    with pytest.raises(Exception):
         # Nonsensical unit for molar energy
        ThermoDynamics.relative_population_states(delta_molar_energy=(666, C), temp=298.15)
