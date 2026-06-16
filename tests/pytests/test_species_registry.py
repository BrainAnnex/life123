import pytest
import pandas as pd
from life123.species_registry import Species, SpeciesRegistry
from tests.utilities.comparisons import *


def test_constructor_Species():

    s = Species(id="test")
    assert s == Species(id='test', name='test', label='test', sort_order=0, diffusion_rate=None, charge=0, formula='',
                        ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')

    s = Species(id="metK", name="methionine adenosyltransferase", label="MK", molecular_weight=43000, categories=["protein", "enzyme"],
                ec_number="2.5.1.6", diffusion_rate=123,
                metadata={"compartment": "cytosol"})
    assert s == Species(id='metK', name='methionine adenosyltransferase', label='MK', sort_order=0, diffusion_rate=123,
                        charge=0, formula='', ec_number='2.5.1.6',
                        cas_number='', chebi='', locus_tag='', molecular_weight=43000, categories=['protein', 'enzyme'],
                        metadata={'compartment': 'cytosol'}, annotation='', plot_color='')

    s = Species(id="test", diffusion_rate=0)
    assert s == Species(id='test', name='test', label='test', sort_order=0, diffusion_rate=0, charge=0, formula='',
                        ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')

    with pytest.raises(Exception):
        Species(id="test", unknown_arg="I shouldn't be here")   # Non-existent arg

    with pytest.raises(Exception):
        Species(id="test", diffusion_rate=-1)       # Negative diffusion

    with pytest.raises(Exception):
        Species(id="test", diffusion_rate="do I look like a rate??")

    with pytest.raises(Exception):
        Species(id="test", categories=123)          # Bad value





##################  class SpeciesRegistry  ##################

def test_constructor_SpeciesRegistry():

    sr = SpeciesRegistry()
    assert sr.by_id == {}

    with pytest.raises(Exception):
        SpeciesRegistry(ids=["A", "B"], n_species=2)    # Passing too many args


    # `ids` argument
    with pytest.raises(Exception):
        SpeciesRegistry(ids=123)    # Bad value


    sr = SpeciesRegistry(ids="X")

    assert sr.number_of_species() == 1
    assert sr.count == 1
    assert len(sr.by_id) == 1
    assert sr.by_id['X'] == Species(id='X', name='X', label='X', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')


    sr = SpeciesRegistry(ids=["X"])

    assert sr.number_of_species() == 1
    assert sr.count == 1
    assert len(sr.by_id) == 1
    assert sr.by_id['X'] == Species(id='X', name='X', label='X', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')


    sr = SpeciesRegistry(ids=("X", "Y"))
    assert sr.number_of_species() == 2
    assert sr.count == 2
    assert len(sr.by_id) == 2
    assert list(sr.by_id.keys()) == ["X", "Y"]
    assert sr.by_id['X'] == Species(id='X', name='X', label='X', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')
    assert sr.by_id['Y'] == Species(id='Y', name="Y", label='Y', sort_order=1,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')


    # `species` argument
    s1 = Species(id="X")
    s2 = Species(id="Y")

    with pytest.raises(Exception):
        SpeciesRegistry(species=123)    # Bad value

    sr = SpeciesRegistry(species=s1)
    assert sr.number_of_species() == 1
    assert sr.count == 1
    assert len(sr.by_id) == 1
    assert sr.by_id['X'] == Species(id='X', name='X', label='X', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')

    sr = SpeciesRegistry(species=[s1])
    assert sr.number_of_species() == 1
    assert sr.count == 1
    assert len(sr.by_id) == 1
    assert sr.by_id['X'] == Species(id='X', name='X', label='X', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')

    sr = SpeciesRegistry(species=(s1, s2))
    assert sr.number_of_species() == 2
    assert sr.count == 2
    assert len(sr.by_id) == 2
    assert list(sr.by_id.keys()) == ["X", "Y"]
    assert sr.by_id['X'] == Species(id='X', name='X', label='X', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')
    assert sr.by_id['Y'] == Species(id='Y', name="Y", label='Y', sort_order=1,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')


    # `n_species` argument
    with pytest.raises(Exception):
        SpeciesRegistry(n_species="X")

    with pytest.raises(Exception):
        SpeciesRegistry(n_species=0)

    sr = SpeciesRegistry(n_species=1)

    assert sr.number_of_species() == 1
    assert sr.count == 1
    assert len(sr.by_id) == 1
    assert sr.by_id['A'] == Species(id='A', name='A', label='A', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')

    sr = SpeciesRegistry(n_species=2)
    assert sr.number_of_species() == 2
    assert sr.count == 2
    assert len(sr.by_id) == 2
    assert list(sr.by_id.keys()) == ["A", "B"]
    assert sr.by_id['A'] == Species(id='A', name='A', label='A', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')
    assert sr.by_id['B'] == Species(id='B', name="B", label='B', sort_order=1,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')



def test_get_species():

    sr = SpeciesRegistry()

    with pytest.raises(Exception):
        sr.get_species("A")         # Doesn't exist

    with pytest.raises(Exception):
        sr.get_species(id=123)      # id isn't a string

    result = sr.add_species(id="A")
    assert result == Species(id='A', name='A', label='A', sort_order=0,
                            diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                            locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')

    result = sr.add_species(id="B", name="B's long name", plot_color='red')
    assert result == Species(id='B', name="B's long name", label='B', sort_order=1,
                        diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                        locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='red')



def test_get_species_id():

    sr = SpeciesRegistry()

    with pytest.raises(Exception):
        sr.get_species_id(0)

    sr.add_species(id="A")
    assert sr.get_species_id(0) == "A"

    with pytest.raises(Exception):
        sr.get_species_id(1)        # Doesn't exist

    sr.add_species(id="B")
    assert sr.get_species_id(1) == "B"  # Now it exists
    assert sr.get_species_id(0) == "A"



def test_number_of_species():

    sr = SpeciesRegistry()
    assert sr.number_of_species() == 0

    sr.add_species(id="A")
    assert sr.number_of_species() == 1

    sr.add_species(id="B")
    assert sr.number_of_species() == 2



def test_add_species():

    sr = SpeciesRegistry()

    with pytest.raises(Exception):
        sr.add_species(id=123)      # Bad value

    result = sr.add_species(id="A")
    assert sr.number_of_species() == 1
    assert sr.count == 1
    assert len(sr.by_id) == 1
    species_A = Species(id='A', name='A', label='A', sort_order=0,
                        diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                        locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')
    assert sr.by_id['A'] == species_A
    assert result == species_A


    result = sr.add_species(id="B", name="B's long name")
    assert sr.number_of_species() == 2
    assert sr.count == 2
    assert len(sr.by_id) == 2
    species_B = Species(id='B', name="B's long name", label='B', sort_order=1,
                        diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                        locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')
    assert sr.by_id['B'] == species_B
    assert result == species_B


    with pytest.raises(Exception):
        sr.add_species(id="B")      # Already exists

    assert sr.add_species(id="B", skip_duplicates=True) is None     # Already exists



def test_as_recordset():

    sr = SpeciesRegistry()
    assert sr.as_recordset() == []

    sr.add_species(id="A")
    result = sr.as_recordset()
    assert result == [
                        {   'id': 'A', 'name': 'A', 'label': 'A', 'sort_order': 0, 'diffusion_rate': None, 'charge': 0,
                            'formula': '', 'ec_number': '', 'cas_number': '', 'chebi': '', 'locus_tag': '',
                            'molecular_weight': None, 'categories': [], 'metadata': {}, 'annotation': '', 'plot_color': ''
                        }
                     ]

    sr.add_species(id="NAD", cas_number="53-84-9")
    result = sr.as_recordset()
    assert result == [
                        {   'id': 'A', 'name': 'A', 'label': 'A', 'sort_order': 0, 'diffusion_rate': None, 'charge': 0,
                            'formula': '', 'ec_number': '', 'cas_number': '', 'chebi': '', 'locus_tag': '',
                            'molecular_weight': None, 'categories': [], 'metadata': {}, 'annotation': '', 'plot_color': ''
                        },
                        {   'id': 'NAD', 'name': 'NAD', 'label': 'NAD', 'sort_order': 1, 'diffusion_rate': None, 'charge': 0,
                            'formula': '', 'ec_number': '', 'cas_number': '53-84-9', 'chebi': '', 'locus_tag': '',
                            'molecular_weight': None, 'categories': [], 'metadata': {}, 'annotation': '', 'plot_color': ''
                        }
                     ]



def test_as_dataframe():

    sr = SpeciesRegistry()
    result = sr.as_dataframe()
    assert result.empty

    sr.add_species(id="A")
    result = sr.as_dataframe()

    expected_recordset = [
                            {   'id': 'A', 'name': 'A', 'label': 'A', 'sort_order': 0, 'diffusion_rate': None, 'charge': 0,
                                'formula': '', 'ec_number': '', 'cas_number': '', 'chebi': '', 'locus_tag': '',
                                'molecular_weight': None, 'categories': [], 'metadata': {}, 'annotation': '', 'plot_color': ''
                            }
                         ]
    expected_df = pd.DataFrame(expected_recordset)
    assert expected_df.equals(result)


    sr.add_species(id="NAD", cas_number="53-84-9", categories=["enzyme"])
    result = sr.as_dataframe()

    expected_recordset = [
                        {   'id': 'A', 'name': 'A', 'label': 'A', 'sort_order': 0, 'diffusion_rate': None, 'charge': 0,
                            'formula': '', 'ec_number': '', 'cas_number': '', 'chebi': '', 'locus_tag': '',
                            'molecular_weight': None, 'categories': [], 'metadata': {}, 'annotation': '', 'plot_color': ''
                        },
                        {   'id': 'NAD', 'name': 'NAD', 'label': 'NAD', 'sort_order': 1, 'diffusion_rate': None, 'charge': 0,
                            'formula': '', 'ec_number': '', 'cas_number': '53-84-9', 'chebi': '', 'locus_tag': '',
                            'molecular_weight': None, 'categories': ["enzyme"], 'metadata': {}, 'annotation': '', 'plot_color': ''
                        }
                     ]
    expected_df = pd.DataFrame(expected_recordset)
    assert expected_df.equals(result)

    sr.by_id["A"].sort_order = 3

    result = sr.as_dataframe(sort=True)
    expected_recordset = [
                        {   'id': 'NAD', 'name': 'NAD', 'label': 'NAD', 'sort_order': 1, 'diffusion_rate': None, 'charge': 0,
                            'formula': '', 'ec_number': '', 'cas_number': '53-84-9', 'chebi': '', 'locus_tag': '',
                            'molecular_weight': None, 'categories': ["enzyme"], 'metadata': {}, 'annotation': '', 'plot_color': ''
                        },
                        {   'id': 'A', 'name': 'A', 'label': 'A', 'sort_order': 3, 'diffusion_rate': None, 'charge': 0,
                            'formula': '', 'ec_number': '', 'cas_number': '', 'chebi': '', 'locus_tag': '',
                            'molecular_weight': None, 'categories': [], 'metadata': {}, 'annotation': '', 'plot_color': ''
                        }
                     ]     # Reversed order, with updated the `sort_order` of 'A`
    expected_df = pd.DataFrame(expected_recordset)
    assert expected_df.equals(result)

    #print()
    #print(result)




##################  PRIVATE METHODS  ##################


def test__generate_generic_names():
    sr = SpeciesRegistry()

    assert sr._generate_generic_names(1) == ["A"]
    assert sr._generate_generic_names(2) == ["A", "B"]
    assert sr._generate_generic_names(26) == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                                                     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    assert sr._generate_generic_names(27) == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                                                     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                                                     'Z2']
    assert sr._generate_generic_names(28) == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                                                     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                                                     'Z2', 'Z3']
