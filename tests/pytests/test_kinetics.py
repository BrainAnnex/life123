import numpy as np
from life123.reactions import Stoichiometry
from life123.kinetics import MassAction_Model, Custom_Model



#####################  Class MassAction_Model  #####################

def rate_MassAction_Model():

    # Reaction  P <-> R
    st = Stoichiometry({"P": -1, "R": 1})

    ma = MassAction_Model(stoichiometry=st)
    ma.set_parameters({"kF": 7., "kR": 3.})

    result = ma.rate(conc_dict={"P": 6, "R": 11.})
    assert np.allclose(result, 7. * 6 - 3. * 11.)   # 9.


    # Reaction  2A <-> B
    st = Stoichiometry({"A": -2, "B": 1})

    ma = MassAction_Model(stoichiometry=st)
    ma.set_parameters({"kF": 5., "kR": 2.})

    result = ma.rate(conc_dict={"A": 4.5, "B": 6.})
    assert np.allclose(result, 5. * 4.5 **2 - 2. * 6.)      # 89.25

    ma.overwrite_parameters({"kF": 5.})     # No reverse reaction
    result = ma.rate(conc_dict={"A": 4.5, "B": 6.})
    assert np.allclose(result, 5. * 4.5 **2)                # 101.25   (irreversible reaction)

    ma.overwrite_parameters({"kF": 3., "kR": 2.})       # A different kR from the first run
    result = ma.rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 59.)     # 3. * 5 **2 - 2. * 8.


    # Reaction  B <-> 2C
    st = Stoichiometry({"B": -1, "C": 2})

    ma = MassAction_Model(stoichiometry=st)
    ma.set_parameters({"kF": 4., "kR": 2.})

    result = ma.rate(conc_dict={"B": 5., "C": 4})
    assert np.allclose(result, 4. * 5. - 2. * 4. **2)       # -12.0


    # Reaction A + B <-> C + D , HYPOTHETICALLY with mass-action kinetics (1st-order kinetics for each species)
    st = Stoichiometry({"A": -1, "B": -1, "C": 1, "D": 1})

    ma = MassAction_Model(stoichiometry=st)
    ma.set_parameters({"kF": 5., "kR": 2.})

    result = ma.rate( conc_dict={"A": 5., "B": 8., "C": 15., "D": 7.})
    assert np.allclose(result, -10.)





#####################  Class Custom_Model  #####################


def test_kinetic_rate_first_order():
    # All reactions below are hypothetically modeled
    # as having 1st-order kinetics with respect to each of the involved species


    # Reaction A <-> B
    st = Stoichiometry({"A": -1, "B": 1})
    result = Custom_Model.kinetic_rate_first_order(stoichiometry=st,
                                                   kinetic_parameters={"kF": 20., "kR": 2.},
                                                   conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5. - 2. * 8.)  # 84.0

    result = Custom_Model.kinetic_rate_first_order(stoichiometry=st,
                                                   kinetic_parameters={"kF": 20.},  # Irreversible
                                                   conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5)             # 100.0


    # Reaction A + B <-> C + D , with 1st-order kinetics for each species
    st = Stoichiometry({"A": -1, "B": -1, "C": 1, "D": 1})

    result = Custom_Model.kinetic_rate_first_order(stoichiometry=st,
                                                   kinetic_parameters={"kF": 10.},  # Irreversible
                                                   conc_dict={"A": 2, "B": 4, "C": 0, "D": 3})
    assert np.allclose(result,  10. * 2 * 4 )  # 80

    result = Custom_Model.kinetic_rate_first_order(stoichiometry=st,
                                                   kinetic_parameters={"kF": 5., "kR": 2.},
                                                   conc_dict={"A": 3.5, "B": 9., "C": 11., "D": 7.})
    assert np.allclose(result,  5. * 3.5 * 9. - 2. * 11. * 7.)  # 3.5

    result = Custom_Model.kinetic_rate_first_order(stoichiometry=st,
                                                   kinetic_parameters={"kF": 5., "kR": 2.},
                                                   conc_dict={"A": 5., "B": 8., "C": 15., "D": 7.})
    assert np.allclose(result,  -10.)
