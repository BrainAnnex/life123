import pytest
import numpy as np
from life123 import ChemData, Membranes2D



#########   TESTS OF INITIALIZATION, SETTING AND VIEWING    #########

def test_constructor():
    with pytest.raises(Exception):
        Membranes2D()                  # Missing required arguments

    with pytest.raises(Exception):
        Membranes2D(x_bins="I'm not an integer", y_bins=9)

    with pytest.raises(Exception):
        Membranes2D(x_bins=9, y_bins="I'm not an integer")


    with pytest.raises(Exception):
        Membranes2D(x_bins=0, y_bins=9)      # Number of bins must be at least 1

    with pytest.raises(Exception):
        Membranes2D(x_bins=6, y_bins=0)      # Number of bins must be at least 1


    bio = Membranes2D(x_bins=6, y_bins=4)

    assert bio.n_bins_x == 6
    assert bio.n_bins_y == 4
    assert bio.membrane_list == []
    assert bio.permeability == {}


def test_describe_state():
    pass



def test_uses_membranes():
    bio = Membranes2D(x_bins=6, y_bins=4)

    assert not bio.uses_membranes()
    '''
    bio.add_membrane(membranes=[(2, 4)])
    assert bio.uses_membranes()
    
    
    bio = Membranes2D(n_bins=123)
    assert not bio.uses_membranes()

    bio.add_membrane(membranes=[(0, 4)])
    assert bio.uses_membranes()

    bio.add_membrane(membranes=[])
    assert not bio.uses_membranes()
    '''



def test_add_membrane_1():
    bio = Membranes2D(x_bins=15, y_bins=20) # 15 bins in x-dir, numbered 0 thru 14
                                            # 20 bins in x-dir, numbered 0 thru 19
    with pytest.raises(Exception):
        bio.add_membrane()     # Missing required arg

    with pytest.raises(Exception):
        bio.add_membrane(vertex_list=123)           # Wrong arg type

    with pytest.raises(Exception):
        bio.add_membrane(vertex_list=[1, 2, 3])     # Too few elements

    with pytest.raises(Exception):
        bio.add_membrane([1, 2, 3, 4])              # Wrong elements in the list
        
    with pytest.raises(Exception):
        bio.add_membrane([(7, 8), (7, 8, 9), (2, 7), (3,7)])    # Element isn't a pair

    with pytest.raises(Exception):
        bio.add_membrane([(7, 8), (7, "some_string"), (2, 7), (3,7), (3,7)])  # Element isn't a pair of integers

    with pytest.raises(Exception):
        bio.add_membrane([(-1,5), (10,5), (10,8), (3,7)])       # x-value out of bounds

    with pytest.raises(Exception):
        bio.add_membrane([(7,5), (10,5), (15,5), (3,7)])        # x-value out of bounds

    with pytest.raises(Exception):
        bio.add_membrane([(7,5), (10,5), (10,-1), (3,7)])       # y-value out of bounds

    with pytest.raises(Exception):
        bio.add_membrane([(7,5), (10,5), (10,20), (3,7)])       # y-value out of bounds

    with pytest.raises(Exception):
        bio.add_membrane([(7,5), (10,5), (10,5), (10,8)])       # Repeated point

    with pytest.raises(Exception):
        bio.add_membrane([(7,5), (10,5), (11,8), (3,7)])        # Diagonal segment (not axes-aligned)

    with pytest.raises(Exception):
        bio.add_membrane([(7,5), (10,5), (10,8), (12,8)])       # Diagonal segment closing the polygon


    with pytest.raises(Exception):
        bio.add_membrane([(7,5), (12,5), (10,5), (10,12), (7,12)])  # Collinear adjacent segments along x-axis

    with pytest.raises(Exception):
        bio.add_membrane([(7,5), (10,5), (10,8), (10,12), (7,12)])  # Collinear adjacent segments along y-axis


    result = bio.add_membrane([(7,5), (10,5), (10,12), (7,12)])
    assert result == 0
    assert bio.membrane_list == [[(7,5), (10,5), (10,12), (7,12)]]

    result = bio.add_membrane([(14,8), (14,10), (7,10), (7,8)])
    assert result == 1
    assert bio.membrane_list == [
                                    [(7,5), (10,5), (10,12), (7,12)],
                                    [(14,8), (14,10), (7,10), (7,8)]
                                ]


def test_add_membrane_2():
    bio = Membranes2D(x_bins=15, y_bins=20) # 15 bins in x-dir, numbered 0 thru 14
                                            # 20 bins in x-dir, numbered 0 thru 19

    result = bio.add_membrane([(7,5), (10,5), (10,8), (14,8), (14,10), (7,10), (7,5)])
    assert result == 0
    assert bio.membrane_list == [[(7,5), (10,5), (10,8), (14,8), (14,10), (7,10)]]  # The redundant last vertex was dropped



def test_add_membrane_3():
    bio = Membranes2D(x_bins=99, y_bins=99)

    with pytest.raises(Exception):
        # Edges 1 and 4 intersect
        bio._assert_simple_polygon([(0,0), (2,0), (2,7), (10,7), (10,4), (0,4)])

    with pytest.raises(Exception):
        # Edges 3 and 7 overlap
        bio._assert_simple_polygon([(10,25),(20,25),(20,60),(15,60),(15,40),(18,40),(18,30),(15,30),(15,50),(10,50)])

    with pytest.raises(Exception):
        # Edges 3 and 7 have a point in common
        bio._assert_simple_polygon([(10,25),(20,25),(20,60),(15,60),(15,40),(18,40),(18,30),(15,30),(15,40),(10,40)])



def test__assert_simple_polygon():
    bio = Membranes2D(x_bins=99, y_bins=99)

    bio._assert_simple_polygon([(0,0), (10,0), (10,10), (0,10), (0,0)])

    with pytest.raises(Exception):
        # Edges 1 and 4 intersect
        bio._assert_simple_polygon([(0,0), (2,0), (2,7), (10,7), (10,4), (0,4), (0,0)])

    with pytest.raises(Exception):
        # Edges 3 and 7 overlap
        bio._assert_simple_polygon([(10,25),(20,25),(20,60),(15,60),(15,40),(18,40),(18,30),(15,30),(15,50),(10,50),(10,25)])

    with pytest.raises(Exception):
        # Edges 3 and 7 have a point in common
        bio._assert_simple_polygon([(10,25),(20,25),(20,60),(15,60),(15,40),(18,40),(18,30),(15,30),(15,40),(10,40),(10,25)])



def test__segments_intersect():
    bio = Membranes2D(x_bins=99, y_bins=99)

    # 2 horizontal segments
    assert not bio._segments_intersect([(0,0), (10,0)] , [(11,0), (12,0)])
    assert     bio._segments_intersect([(0,0), (10,0)] , [(10,0), (12,0)])
    assert     bio._segments_intersect([(0,0), (10,0)] , [(9,0), (12,0)])
    assert     bio._segments_intersect([(0,0), (10,0)] , [(2,0), (7,0)])

    # 2 vertical segments
    assert not bio._segments_intersect([(8,30), (8,40)], [(8,10), (8,29)])
    assert     bio._segments_intersect([(8,30), (8,40)], [(8,10), (8,30)])
    assert     bio._segments_intersect([(8,30), (8,40)], [(8,10), (8,31)])
    assert     bio._segments_intersect([(8,30), (8,40)], [(8,31), (8,39)])

    # One horizontal segments and a vertical one
    assert not bio._segments_intersect([(30,8), (40,8)], [(35,0), (35,7)])
    assert     bio._segments_intersect([(30,8), (40,8)], [(35,0), (35,8)])
    assert     bio._segments_intersect([(30,8), (40,8)], [(35,0), (35,9)])
    assert     bio._segments_intersect([(0,10), (0,0)] , [(0,0), (10,0)])




def test_set_permeability():
    pass

def test_change_permeability():
    pass
