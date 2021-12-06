from bio_sim import BioSim


def test_initialize_universe():
    BioSim.initialize_universe(n_cells=10, n_species=1)

    assert BioSim.n_cells == 10
    assert BioSim.n_species == 1