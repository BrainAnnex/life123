import pytest
import pandas as pd
from life123.species_registry import Species, SpeciesRegistry, Macromol
from life123.visualization.colors import Colors
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


    s = Species(id="test")
    s.diffusion_rate = 10
    assert s == Species(id='test', name='test', label='test', sort_order=0, diffusion_rate=10, charge=0, formula='',
                        ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')

    with pytest.raises(Exception):
        s.diffusion_rate = -1           # Negative diffusion

    with pytest.raises(Exception):
        Species(id="test", name=123)    # Bad value

    with pytest.raises(Exception):
        Species(id="test", label=["confused"])      # Bad value


    s = Species(id="test", name="       ")          # `id` is used for `name` when `name` is missing
    assert s.name == "test"
    s = Species(id="test", name=None)               # `id` is used for `name` when `name` is missing
    assert s.name == "test"

    s = Species(id="test", label="       ")         # `id` is used for `label` when `label` is missing
    assert s.label == "test"
    s = Species(id="test", label=None)              # `id` is used for `label` when `label` is missing
    assert s.label == "test"





########################  class SpeciesRegistry  ########################

def test_constructor_SpeciesRegistry():

    sr = SpeciesRegistry()
    assert sr.by_id == {}

    with pytest.raises(Exception):
        SpeciesRegistry(id=["A", "B"], n_species=2)    # Passing too many args


    # `ids` argument
    with pytest.raises(Exception):
        SpeciesRegistry(id=123)    # Bad value


    sr = SpeciesRegistry(id="X")

    assert sr.number_of_species() == 1
    assert sr.count == 1
    assert len(sr.by_id) == 1
    assert sr.by_id['X'] == Species(id='X', name='X', label='X', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')


    sr = SpeciesRegistry(id=["X"])

    assert sr.number_of_species() == 1
    assert sr.count == 1
    assert len(sr.by_id) == 1
    assert sr.by_id['X'] == Species(id='X', name='X', label='X', sort_order=0,
                                    diffusion_rate=None, charge=0, formula='', ec_number='', cas_number='', chebi='',
                                    locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='')


    sr = SpeciesRegistry(id=("X", "Y"))
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


    # Test special arguments
    with pytest.raises(Exception):
        SpeciesRegistry(diffusion_rate=11)      # Mismatched size


    sr = SpeciesRegistry(n_species=1, diffusion_rate=123)
    assert sr.by_id['A'] == Species(id='A', name='A', label='A', sort_order=0, diffusion_rate=123, charge=0, formula='',
                                    ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                    categories=[], metadata={}, annotation='', plot_color='')

    sr = SpeciesRegistry(n_species=1, plot_color=['pink'])
    assert sr.by_id['A'] == Species(id='A', name='A', label='A', sort_order=0, diffusion_rate=None, charge=0, formula='',
                                    ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                    categories=[], metadata={}, annotation='', plot_color='pink')

    with pytest.raises(Exception):
        SpeciesRegistry(id="X", diffusion_rate=[1, 2])      # Mismatched size

    sr = SpeciesRegistry(id=["X", "Y"], diffusion_rate=[1, 2])
    assert sr.get_all_values(field="diffusion_rate") == [1, 2]

    sr = SpeciesRegistry(id=("X", "Y"), diffusion_rate=[1, 2], plot_color=("blue", "turquoise"))
    assert sr.get_all_values(field="diffusion_rate") == [1, 2]
    assert sr.get_all_values(field="plot_color") == ["blue", "turquoise"]



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



def test_get_all_species_ids():
    sr = SpeciesRegistry()

    assert sr.get_all_species_ids() == []   # There's nothing yet

    sr.add_species(id="A")
    assert sr.get_all_species_ids() == ["A"]

    sr.add_species(id="B")
    assert sr.get_all_species_ids() == ["A", "B"]



def test_get_species_index():
    sr = SpeciesRegistry()

    with pytest.raises(Exception):
        sr.get_species_index("A")   # Doesn't exist

    sr.add_species(id="A")
    assert sr.get_species_index("A") == 0

    with pytest.raises(Exception):
        sr.get_species_index("B")   # Doesn't exist

    sr.add_species(id="B")
    assert sr.get_species_index("B") == 1
    assert sr.get_species_index("A") == 0



def test_get_value():
    sr = SpeciesRegistry()

    with pytest.raises(Exception):
        sr.get_value(species_id="A", field="xyz")  # No such species exist

    sr.add_species(id="A", ec_number="2.5.1.6")

    with pytest.raises(Exception):
        sr.get_value(species_id="A", field="xyz")  # No such field exist

    assert sr.get_value(species_id="A", field="ec_number") == "2.5.1.6"
    assert sr.get_value(species_id="A", field="charge") == 0
    assert sr.get_value(species_id="A", field="diffusion_rate") is None



def test_get_all_values():
    sr = SpeciesRegistry()

    assert sr.get_all_values(field="charge") == []

    sr.add_species(id="A")
    assert sr.get_all_values(field="charge") == [0]

    sr.add_species(id="B", ec_number="2.5.1.6", diffusion_rate=123)
    assert sr.get_all_values(field="ec_number") == ["", "2.5.1.6"]
    assert sr.get_all_values(field="charge") == [0, 0]
    assert sr.get_all_values(field="diffusion_rate") == [None, 123]



def test_get_max_value():
    sr = SpeciesRegistry()
    
    with pytest.raises(Exception):
        sr.max_value(field="diffusion_rate")       # Empty registry


    sr.add_species(id="A")
    sr.set_value(species_id="A", field="diffusion_rate", value=10)
    result = sr.max_value(field="diffusion_rate")
    assert result == 10

    sr.add_species(id="B")
    sr.set_value(species_id="B", field="diffusion_rate", value=3)
    result = sr.max_value(field="diffusion_rate")
    assert result == 10

    sr.set_value(species_id="B", field="diffusion_rate", value=11)
    result = sr.max_value(field="diffusion_rate")
    assert result == 11

    sr.add_species(id="X")      # This won't get the `diffusion_rate` set; its value will be None
    with pytest.raises(Exception):
        sr.max_value(field="diffusion_rate")

    

def test_missing_values_present():
    sr = SpeciesRegistry()

    assert not sr.has_missing_values(field="diffusion_rate")

    sr.add_species(id="A")
    sr.set_value(species_id="A", field="diffusion_rate", value=10)
    assert not sr.has_missing_values(field="diffusion_rate")

    sr.add_species(id="B")
    sr.set_value(species_id="B", field="diffusion_rate", value=3)
    assert not sr.has_missing_values(field="diffusion_rate")

    sr.set_value(species_id="B", field="diffusion_rate", value=11)
    assert not sr.has_missing_values(field="diffusion_rate")

    sr.add_species(id="X")      # This won't get the `diffusion_rate` set; its value will be None
    assert sr.has_missing_values(field="diffusion_rate")



def test_number_of_species():
    sr = SpeciesRegistry()
    assert sr.number_of_species() == 0

    sr.add_species(id="A")
    assert sr.number_of_species() == 1

    sr.add_species(id="B")
    assert sr.number_of_species() == 2



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



##################  SET VALUES  ##################

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



def test_set_update():
    sr = SpeciesRegistry()

    sr.add_species(id="A")
    sr.update(species_id="A", diffusion_rate=10)
    assert  sr.get_value(species_id="A", field="diffusion_rate") == 10

    with pytest.raises(Exception):
        sr.update(species_id="A", diffusion_rate=-4)     # Bad value

    with pytest.raises(Exception):
        sr.update(species_id="A", categories="I'm not a list")    # Bad value


    sr.add_species(id="B", plot_color="magenta")
    sr.update(species_id="B", categories=["protein", "enzyme"], plot_color="plum")

    assert sr.get_species("A") == Species(id='A', name='A', label='A', sort_order=0, diffusion_rate=10, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                          categories=[], metadata={}, annotation='', plot_color='')

    assert sr.get_species("B") == Species(id='B', name='B', label='B', sort_order=1, diffusion_rate=None, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                          categories=["protein", "enzyme"], metadata={}, annotation='', plot_color='plum')

    with pytest.raises(Exception):
        sr.update(species_id="B", name={"a": 1})     # Bad value

    sr.update(species_id="B", name="      ")    # Missing value; the `id` is used instead
    assert sr.get_species("B") == Species(id='B', name='B', label='B', sort_order=1, diffusion_rate=None, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                          categories=["protein", "enzyme"], metadata={}, annotation='', plot_color='plum')

    with pytest.raises(Exception):
        sr.update(species_id="B", label=1234)     # Bad value

    sr.update(species_id="B", label="      ")    # Missing value; the `id` is used instead
    assert sr.get_species("B") == Species(id='B', name='B', label='B', sort_order=1, diffusion_rate=None, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                          categories=["protein", "enzyme"], metadata={}, annotation='', plot_color='plum')

    with pytest.raises(Exception):
        sr.update(species_id="B")   # Missing arguments



def test_set_value():
    sr = SpeciesRegistry()

    sr.add_species(id="A")

    sr.set_value(species_id="A", field="diffusion_rate", value=10)
    assert  sr.get_value(species_id="A", field="diffusion_rate") == 10

    with pytest.raises(Exception):
        sr.set_value(species_id="A", field="diffusion_rate", value=-4)     # Bad value

    with pytest.raises(Exception):
        sr.set_value(species_id="A", field="categories", value="I'm not a list")    # Bad value

    sr.set_value(species_id="A", field="categories", value=["protein", "enzyme"])
    result = sr.get_value(species_id="A", field="categories")
    assert result == ["protein", "enzyme"]
    result = sr.get_value(species_id="A", field="diffusion_rate")
    assert result == 10


    sr.add_species(id="B", plot_color="magenta")
    sr.set_value(species_id="B", field="plot_color", value="plum")
    result = sr.get_value(species_id="B", field="plot_color")
    assert result == "plum"

    assert sr.get_species("A") == Species(id='A', name='A', label='A', sort_order=0, diffusion_rate=10, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                          categories=["protein", "enzyme"], metadata={}, annotation='', plot_color='')

    assert sr.get_species("B") == Species(id='B', name='B', label='B', sort_order=1, diffusion_rate=None, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                          categories=[], metadata={}, annotation='', plot_color='plum')

    with pytest.raises(Exception):
        sr.set_value(species_id="B", field="name", value={"a": 1})     # Bad value

    sr.set_value(species_id="B", field="name", value="      ")    # Missing value; the `id` is used instead
    assert sr.get_species("B") == Species(id='B', name='B', label='B', sort_order=1, diffusion_rate=None, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                          categories=[], metadata={}, annotation='', plot_color='plum')

    with pytest.raises(Exception):
        sr.set_value(species_id="B", field="label", value=1234)     # Bad value

    sr.set_value(species_id="B", field="label", value="      ")    # Missing value; the `id` is used instead
    assert sr.get_species("B") == Species(id='B', name='B', label='B', sort_order=1, diffusion_rate=None, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None,
                                          categories=[], metadata={}, annotation='', plot_color='plum')


    sr.set_value(species_id="A", field="diffusion_rate", value=None)    # None is a legit value to use
    assert sr.get_value(species_id="A", field="diffusion_rate") is None

    with pytest.raises(Exception):
        sr.set_value(species_id="A", field="Impostor", value=666)           # Field doesn't exist in Species




def test_set_all_values():
    sr = SpeciesRegistry()

    with pytest.raises(Exception):
        sr.set_all_values(field="diffusion_rate", values=[5])  # Mismatch with registry size

    sr.add_species(id="A")
    assert sr.get_value(species_id="A", field="diffusion_rate") is None    # Not yet set
    sr.set_all_values(field="diffusion_rate", values=[5])
    assert sr.get_value(species_id="A", field="diffusion_rate") == 5       # Now set

    sr.add_species(id="B")
    assert sr.get_all_values(field="diffusion_rate") == [5, None]

    with pytest.raises(Exception):
        sr.set_all_values(field="diffusion_rate", values=[5])  # Mismatch with registry size

    sr.set_all_values(field="diffusion_rate", values=[8, 12])
    assert sr.get_all_values(field="diffusion_rate") == [8, 12]

    sr.set_all_values(field="diffusion_rate", values=[None, 666])
    assert sr.get_all_values(field="diffusion_rate") == [None, 666]

    sr.set_all_values(field="plot_color", values=["red", "pink"])
    assert sr.get_all_values(field="plot_color") == ["red", "pink"]

    assert sr.get_species("A") == Species(id='A', name='A', label='A', sort_order=0, diffusion_rate=None, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='red')

    assert sr.get_species("B") == Species(id='B', name='B', label='B', sort_order=1, diffusion_rate=666, charge=0, formula='',
                                          ec_number='', cas_number='', chebi='', locus_tag='', molecular_weight=None, categories=[], metadata={}, annotation='', plot_color='pink')



def test_assign_colors():
    sr = SpeciesRegistry()

    with pytest.raises(Exception):
        sr.assign_colors(123)

    with pytest.raises(Exception):
        sr.assign_colors("Unknown")


    # Make a note of the first 3 default colors
    default_col_1, default_col_2, default_col_3 = Colors.assign_default_colors(3)

    sr = SpeciesRegistry(id="A")
    sr.set_value(species_id="A", field="plot_color", value="red")
    result = sr.assign_colors("A")
    assert result == ["red"]
    assert sr.get_value(species_id="A", field="plot_color") == "red"

    sr = SpeciesRegistry(id="A")
    result = sr.assign_colors(["A"])
    assert result == [default_col_1]
    assert sr.get_value(species_id="A", field="plot_color") == default_col_1   # Was made into a permanent assignment


    sr = SpeciesRegistry(id=["A", "B"])
    sr.set_all_values(field="plot_color", values=["red", "blue"])
    result = sr.assign_colors(["A", "B"])
    assert result == ["red", "blue"]
    assert sr.get_all_values(field="plot_color") == ["red", "blue"]

    sr = SpeciesRegistry(id="A")
    sr.set_value(species_id="A", field="plot_color", value="red")
    sr.add_species(id="B")
    result = sr.assign_colors(["A", "B"])
    assert result == ["red", default_col_1]
    assert sr.get_all_values(field="plot_color") == ["red", default_col_1]

    sr = SpeciesRegistry(id=["A", "B"])
    sr.set_value(species_id="B", field="plot_color", value="yellow")
    result = sr.assign_colors(["A", "B"])
    assert result == [default_col_1, "yellow"]
    assert sr.get_all_values(field="plot_color") == [default_col_1, "yellow"]

    sr = SpeciesRegistry(id=["A", "B"])
    result = sr.assign_colors(["A", "B"])
    assert result == [default_col_1, default_col_2]
    assert sr.get_all_values(field="plot_color") == [default_col_1, default_col_2]

    sr = SpeciesRegistry(id=["A", "B"])
    result = sr.assign_colors()
    assert result == [default_col_1, default_col_2]
    assert sr.get_all_values(field="plot_color") == [default_col_1, default_col_2]

    sr = SpeciesRegistry(id=["A", "B", "C"])
    sr.set_value(species_id="A", field="plot_color", value="red")
    sr.set_value(species_id="C", field="plot_color", value="yellow")
    result = sr.assign_colors(["A", "B", "C"])
    assert result == ["red", default_col_1, "yellow"]
    assert sr.get_all_values(field="plot_color") == ["red", default_col_1, "yellow"]

    sr = SpeciesRegistry(id=["A", "B", "C", "D", "E"])
    result = sr.assign_colors(["B", "C", "D"])
    assert result == [default_col_1, default_col_2, default_col_3]
    assert sr.get_value(species_id="B", field="plot_color") == default_col_1
    assert sr.get_value(species_id="C", field="plot_color") == default_col_2
    assert sr.get_value(species_id="D", field="plot_color") == default_col_3

    sr = SpeciesRegistry(id=["A", "B", "C", "D", "E"])
    sr.set_value(species_id="B", field="plot_color", value="red")
    sr.set_value(species_id="D", field="plot_color", value="yellow")
    result = sr.assign_colors(["B", "C", "D"])
    assert result == ["red", default_col_1, "yellow"]
    assert sr.get_value(species_id="B", field="plot_color") == "red"
    assert sr.get_value(species_id="C", field="plot_color") == default_col_1
    assert sr.get_value(species_id="D", field="plot_color") == "yellow"



def test_get_color_mapping_by_index():
    sr = SpeciesRegistry()
    assert sr.get_color_mapping_by_index() == {}

    sr = SpeciesRegistry(n_species=3)
    sr.set_all_values(field="plot_color", values=['red', 'green', 'blue'])
    assert sr.get_color_mapping_by_index() == {0: "red", 1: "green", 2: "blue"}





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
                                                     
                                       
                                                     



########################  MACRO-MOLECULES  ########################


def test_add_macromolecules():
    chem_data = Macromol()

    chem_data.add_macromolecules(["M1", "M2"])
    assert chem_data.get_macromolecules() == ["M1", "M2"]

    chem_data.add_macromolecules("M3")
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3"]

    chem_data.add_macromolecules(["M4", "M5"])
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3", "M4", "M5"]

    chem_data.add_macromolecules("M2")          # Redundant
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3", "M4", "M5"]

    chem_data.add_macromolecules(["M1", "M6"])  # Partially redundant
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3", "M4", "M5", "M6"]



def test_get_macromolecules():
    chem_data = Macromol()

    assert chem_data.get_macromolecules() == []
    chem_data.add_macromolecules(["M1", "M2"])
    assert chem_data.get_macromolecules() == ["M1", "M2"]
    chem_data.add_macromolecules(["M3", "M1"])
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3"]



def test_set_binding_site_affinity():
    sr = SpeciesRegistry(id=["A", "B", "ZZZ"])

    chem_data = Macromol()          # TODO: probably pass the SpeciesRegistry as argument
    chem_data.add_macromolecules(["M1", "M2"])

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=3)
    result = chem_data.get_binding_site_affinity(macromolecule="M1", site_number=1)
    assert result == ("A", 3)
    assert result.chemical == "A"
    assert result.Kd == 3

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=5)
    assert chem_data.get_binding_site_affinity("M1", site_number=1) == ("A", 3)
    assert chem_data.get_binding_site_affinity("M1", site_number=2) == ("B", 5)

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="B", Kd=11)
    assert chem_data.get_binding_site_affinity("M1", site_number=1) == ("A", 3)
    assert chem_data.get_binding_site_affinity("M1", site_number=2) == ("B", 5)
    assert chem_data.get_binding_site_affinity("M2", site_number=1) == ("B",11)

    # Over-write previous value of affinity of "A" to site 1 of macromolecule "M1"
    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=999)
    assert chem_data.get_binding_site_affinity("M1", site_number=1) == ("A", 999)
    assert chem_data.get_binding_site_affinity("M1", site_number=2) == ("B", 5)
    assert chem_data.get_binding_site_affinity("M2", site_number=1) == ("B",11)

    # Attempting to associate another chemical to site 1 of macromolecule "M1", will result in error
    with pytest.raises(Exception):
        chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="ZZZ", Kd=999)

    chem_data.set_binding_site_affinity(macromolecule="M3", site_number=10, ligand="A", Kd= 4)
    assert chem_data.get_binding_site_affinity("M3", site_number=10) == ("A", 4)
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3"]     # "M3" got automatically added
    assert chem_data.get_binding_site_affinity("M1", site_number=1) == ("A", 999)
    assert chem_data.get_binding_site_affinity("M1", site_number=2) == ("B", 5)
    assert chem_data.get_binding_site_affinity("M2", site_number=1) == ("B",11)

    # Unknown chemical "X"
    with pytest.raises(Exception):
        chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="X", Kd=100)



def test_get_binding_site_affinity():
    pass  # TODO



def test_get_binding_sites():
    sr = SpeciesRegistry(id=["A", "B"])

    chem_data = Macromol()
    chem_data.add_macromolecules(["M1", "M2"])

    with pytest.raises(Exception):
        chem_data.get_binding_sites("M9999")    # Unknown macromolecule

    assert chem_data.get_binding_sites("M1") == []
    assert chem_data.get_binding_sites("M2") == []

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=3)
    assert chem_data.get_binding_sites("M1") == [1]
    assert chem_data.get_binding_sites("M2") == []

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=5)
    assert chem_data.get_binding_sites("M1") == [1, 2]
    assert chem_data.get_binding_sites("M2") == []

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="B", Kd=11)
    assert chem_data.get_binding_sites("M1") == [1, 2]
    assert chem_data.get_binding_sites("M2") == [1]

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=3, ligand="B", Kd=102)
    assert chem_data.get_binding_sites("M1") == [1, 2]
    assert chem_data.get_binding_sites("M2") == [1, 3]



def test_get_binding_sites_and_ligands():
    sr = SpeciesRegistry(id=["A", "B"])

    chem_data = Macromol()
    chem_data.add_macromolecules(["M1", "M2"])

    with pytest.raises(Exception):
        chem_data.get_binding_sites_and_ligands("M9999")    # Unknown macromolecule

    assert chem_data.get_binding_sites_and_ligands("M1") == {}
    assert chem_data.get_binding_sites_and_ligands("M2") == {}

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=3)
    assert chem_data.get_binding_sites_and_ligands("M1") == {1: "A"}
    assert chem_data.get_binding_sites_and_ligands("M2") == {}

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=5)
    assert chem_data.get_binding_sites_and_ligands("M1") == {1: "A", 2: "B"}
    assert chem_data.get_binding_sites_and_ligands("M2") == {}

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="B", Kd=11)
    assert chem_data.get_binding_sites_and_ligands("M1") == {1: "A", 2: "B"}
    assert chem_data.get_binding_sites_and_ligands("M2") == {1: "B"}

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=3, ligand="B", Kd=102)
    assert chem_data.get_binding_sites_and_ligands("M1") == {1: "A", 2: "B"}
    assert chem_data.get_binding_sites_and_ligands("M2") == {1: "B", 3: "B"}



def test_get_ligand_name():
    sr = SpeciesRegistry(id=["A", "B"])

    chem_data = Macromol()
    chem_data.add_macromolecules(["M1", "M2"])

    with pytest.raises(Exception):
        chem_data.get_ligand_name(macromolecule="M9999", site_number=123)    # Unknown macromolecule

    with pytest.raises(Exception):
        chem_data.get_ligand_name(macromolecule="M1", site_number=1)    # M1 has no sites associated to it

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="A", Kd=3.1)

    with pytest.raises(Exception):
        chem_data.get_ligand_name(macromolecule="M1", site_number=1)    # No site number 1 on M1

    assert chem_data.get_ligand_name(macromolecule="M1", site_number=2) == "A"

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="B", Kd=6.)
    assert chem_data.get_ligand_name(macromolecule="M1", site_number=1) == "B"

    sr.add_species("C")
    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=5, ligand="C", Kd=11.4)
    assert chem_data.get_ligand_name(macromolecule="M2", site_number=5) == "C"



def test_show_binding_affinities():
    pass



def test_reset_macromolecule():
    pass   # TODO


def test_clear_macromolecules():
    chem_data = Macromol()
    chem_data.add_macromolecules(["M1", "M2"])
    chem_data.clear_macromolecules()
    assert chem_data.macromolecules == []
    assert chem_data.binding_sites == {}
