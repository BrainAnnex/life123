import pytest
import numpy as np
from modules.numerical.numerical import Numerical as num
import pandas as pd


def test_curve_intersect_interpolate():
    df = pd.DataFrame({"TIME": [1, 2, 3], "A": [2,4,10], "B": [12,4,0]})
    result = num.curve_intersect_interpolate(df, row_index=1, x="TIME", var1="A", var2="B")
    assert np.allclose(result, (2, 4))

    df = pd.DataFrame({"SYSTEM TIME": [10, 20, 30], "p": [2,4,10], "q": [12,4,0]})
    result = num.curve_intersect_interpolate(df, row_index=1, x="SYSTEM TIME", var1="p", var2="q")
    assert np.allclose(result, (20, 4))

    df = pd.DataFrame({"SYSTEM TIME": [10, 20, 30], "p": [2,4,10], "q": [12,5,0]})
    result = num.curve_intersect_interpolate(df, row_index=1, x="SYSTEM TIME", var1="p", var2="q")
    assert np.allclose(result, (20.90909090909091, 4.545454545454546))
    alternate_result = num.line_intersect((20, 4), (30, 10), (20, 5), (30, 0))  # Worked out by hand
    assert np.allclose(result, alternate_result)

    df = pd.DataFrame({"SYSTEM TIME": [10, 20, 30], "p": [2,4,10], "q": [12,1,0]})
    result = num.curve_intersect_interpolate(df, row_index=1, x="SYSTEM TIME", var1="p", var2="q")
    assert np.allclose(result, (17.692307692307693, 3.5384615384615383))
    alternate_result = num.line_intersect((10, 2), (20, 4), (10, 12), (20, 1))
    assert np.allclose(result, alternate_result)

    df = pd.DataFrame({"x-coord": [5, 6, 10], "line 1": [-2,-2,-2], "line 2": [2,0,-5]})
    result = num.curve_intersect_interpolate(df, row_index=1, x="x-coord", var1="line 1", var2="line 2")
    assert np.allclose(result, (7.6, -2.0))
    alternate_result = num.line_intersect((6, -2), (10, -2), (6, 0), (10, -5))
    assert np.allclose(result, alternate_result)

    df = pd.DataFrame({"x coord": [5, 7, 15], "L1": [-6,-3,3], "L2": [1,-2,-0]})
    result = num.curve_intersect_interpolate(df, row_index=1, x="x coord", var1="L1", var2="L2")
    assert np.allclose(result, (9.0, -1.5))
    alternate_result = num.line_intersect((7, -3), (15, 3), (7, -2), (15, 0))
    assert np.allclose(result, alternate_result)

    # Extra row at end of dataframe
    df = pd.DataFrame({"x coord": [5, 7, 15, 99], "L1": [-6,-3,3,99], "L2": [1,-2,-0,99]})
    result = num.curve_intersect_interpolate(df, row_index=1, x="x coord", var1="L1", var2="L2")
    assert np.allclose(result, (9.0, -1.5))

    # Extra row at start of dataframe
    df = pd.DataFrame({"x coord": [-123, 5, 7, 15], "L1": [-123,-6,-3,3], "L2": [-123,1,-2,-0]})
    result = num.curve_intersect_interpolate(df, row_index=2, x="x coord", var1="L1", var2="L2")
    assert np.allclose(result, (9.0, -1.5))

    # No actual intersection
    df = pd.DataFrame({"x coord": [8, 10, 15], "L1": [-3,-1,-6], "L2": [3,1,30]})
    result = num.curve_intersect_interpolate(df, row_index=1, x="x coord", var1="L1", var2="L2")
    assert np.allclose(result, (10, 0.))


    # Error cases
    df = pd.DataFrame({"TIME": [1, 2, 3], "A": [4,4,4], "B": [4,4,4]})
    with pytest.raises(Exception):
        num.curve_intersect_interpolate(df, row_index=1, x="TIME", var1="A", var2="B")  # Extended overlap on both sides

    df = pd.DataFrame({"TIME": [1, 2, 3], "A": [10,12,30], "B": [10,12,40]})
    with pytest.raises(Exception):
        num.curve_intersect_interpolate(df, row_index=1, x="TIME", var1="A", var2="B")  # Extended overlap on left segments

    df = pd.DataFrame({"TIME": [1, 2, 3], "A": [5,12,30], "B": [-8,12,30]})
    with pytest.raises(Exception):
        num.curve_intersect_interpolate(df, row_index=1, x="TIME", var1="A", var2="B")  # Extended overlap on right segments


    df = pd.DataFrame({"SYSTEM TIME": [5, 6, 10], "A": [-2,-2,-2], "B": [2,0,-5]})
    # Unknown column names
    with pytest.raises(Exception):
        num.curve_intersect_interpolate(df, row_index=1, x="unknown", var1="A", var2="B")
    with pytest.raises(Exception):
        num.curve_intersect_interpolate(df, row_index=1, x="SYSTEM TIME", var1="unknown", var2="B")
    with pytest.raises(Exception):
        num.curve_intersect_interpolate(df, row_index=1, x="SYSTEM TIME", var1="A", var2="unknown")

    df = pd.DataFrame({"TIME": [1, 2], "A": [10,12], "B": [80,100]})
    with pytest.raises(Exception):
        num.curve_intersect_interpolate(df, row_index=1, x="TIME", var1="A", var2="B")  # Too few rows


    df = pd.DataFrame({"x coord": [-2, 0, 5], "L1": [-3,1,-6], "L2": [3,-1,30]})
    # Multiple intersections
    with pytest.raises(Exception):
        num.curve_intersect_interpolate(df, row_index=1, x="x coord", var1="L1", var2="L2")



def test_intersecting_curves():
    with pytest.raises(Exception):
        num.intersecting_curves_NOT_USED_TO_DITCH(times=(0, 2, 2), y1=(1, 2, 3), y2=(4, 5, 6)) # repeated times
        num.intersecting_curves_NOT_USED_TO_DITCH(times=(1, 2, 0), y1=(1, 2, 3), y2=(4, 5, 6)) # non-consecutive times



def test_segment_intersect():
    # Same points as in test_line_intersect()

    # 45-degree angles
    result = num.segment_intersect(t=(-2, 2), y1=(-1, 1), y2=(1, -1))
    assert np.allclose(result, [0., 0.])
    result = num.segment_intersect(t=(-2, 2), y2=(-1, 1), y1=(1, -1))
    assert np.allclose(result, [0., 0.])

    # x-axis (t value) stretch by 2 and shift by 14
    result = num.segment_intersect(t=(10, 18), y1=(-1, 1), y2=(1, -1))
    assert np.allclose(result, [14., 0.])
    result = num.segment_intersect(t=(10, 18), y2=(-1, 1), y1=(1, -1))
    assert np.allclose(result, [14., 0.])

    # y-axis shift by 10
    result = num.segment_intersect(t=(10, 18), y1=(9, 11), y2=(11, 9))
    assert np.allclose(result, [14., 10.])
    result = num.segment_intersect(t=(10, 18), y2=(9, 11), y1=(11, 9))
    assert np.allclose(result, [14., 10.])

    # "3,4,5" right triangles
    result = num.segment_intersect(t=(0, 8), y1=(0, 0), y2=(3, -3))
    assert np.allclose(result, [4., 0.])

    # x-axis (t value) shift by -10
    result = num.segment_intersect(t=(-10, -2), y1=(0, 0), y2=(3, -3))
    assert np.allclose(result, [-6., 0.])

    # y-axis shift by 100
    result = num.segment_intersect(t=(-10, -2), y1=(100, 100), y2=(103, 97))
    assert np.allclose(result, [-6., 100.])

    # Other examples
    result = num.segment_intersect(t=(2, 3), y1=(10, 11), y2=(12, 10))
    assert np.allclose(result, [2.666666, 10.666666])

    result = num.segment_intersect(t=(20,30), y1=(4,10), y2=(5,0))
    assert np.allclose(result, (20.90909090909091, 4.545454545454546))

    result = num.segment_intersect(t=(10,20), y1=(2,4), y2=(12,1))
    assert np.allclose(result, (17.692307692307693, 3.5384615384615383))

    with pytest.raises(Exception):
        num.segment_intersect(t=(2, 3), y1=(10, 10), y2=(12, 12))   # Parallel segments
        num.segment_intersect(t=(2, 3), y2=(10, 10), y1=(12, 12))   # Parallel segments



def test_segment_intersect_OLD():       # TODO: ditch
    # Same points as in test_line_intersect()

    # 45-degree angles
    result = num.segment_intersect_OLD(t=(-2, 2), y_rise=(-1, 1), y_drop=(1, -1))
    assert np.allclose(result, [0., 0.])

    # x-axis (t value) stretch by 2 and shift by 14
    result = num.segment_intersect_OLD(t=(10, 18), y_rise=(-1, 1), y_drop=(1, -1))
    assert np.allclose(result, [14., 0.])

    # y-axis shift by 10
    result = num.segment_intersect_OLD(t=(10, 18), y_rise=(9, 11), y_drop=(11, 9))
    assert np.allclose(result, [14., 10.])

    # "3,4,5" right triangles
    result = num.segment_intersect_OLD(t=(0, 8), y_rise=(0, 0), y_drop=(3, -3))
    assert np.allclose(result, [4., 0.])

    # x-axis (t value) shift by -10
    result = num.segment_intersect_OLD(t=(-10, -2), y_rise=(0, 0), y_drop=(3, -3))
    assert np.allclose(result, [-6., 0.])

    # y-axis shift by 100
    result = num.segment_intersect_OLD(t=(-10, -2), y_rise=(100, 100), y_drop=(103, 97))
    assert np.allclose(result, [-6., 100.])

    result = num.segment_intersect_OLD(t=(2, 3), y_rise=(10, 11), y_drop=(12, 10))
    assert np.allclose(result, [2.666666, 10.666666])

    with pytest.raises(Exception):
        num.segment_intersect_OLD(t=(2, 3), y_rise=(10, 10), y_drop=(12, 12))   # Parallel segments



def test_line_intersect():
    # Same points as in test_segment_intersect()

    # 45-degree angles
    result = num.line_intersect((-2,-1), (2, 1), (-2,1), (2,-1))
    assert np.allclose(result, [0., 0.])

    # x-axis (t value) stretch by 2 and shift by 14
    result = num.line_intersect((10,-1), (18, 1), (10,1), (18,-1))
    assert np.allclose(result, [14., 0.])

    # y-axis shift by 10
    result = num.line_intersect((10,9), (18, 11), (10,11), (18,9))
    assert np.allclose(result, [14., 10.])

    # "3,4,5" right triangles
    result = num.line_intersect((0,0), (8, 0), (0,3), (8,-3))
    assert np.allclose(result, [4., 0.])

    # x-axis (t value) shift by -10
    result = num.line_intersect((-10,0), (-2, 0), (-10,3), (-2,-3))
    assert np.allclose(result, [-6., 0.])

    # y-axis shift by 100
    result = num.line_intersect((-10,100), (-2, 100), (-10,103), (-2,97))
    assert np.allclose(result, [-6., 100.])

    # Other examples
    result = num.line_intersect((2,10), (3, 11), (2,12), (3,10))
    assert np.allclose(result, [2.666666, 10.666666])

    result = num.line_intersect((20, 4), (30, 10), (20, 5), (30, 0))
    assert np.allclose(result, (20.90909090909091, 4.545454545454546))

    result = num.line_intersect((10, 2), (20, 4), (10, 12), (20, 1))
    assert np.allclose(result, (17.692307692307693, 3.5384615384615383))

    assert num.line_intersect((2, 10), (3, 10), (2, 12), (3, 12)) is None   # Parallel segments



def test_deep_flatten():
    assert num.deep_flatten(5) == [5]
    assert num.deep_flatten([1,2]) == [1,2]
    assert num.deep_flatten( (1,2) ) == [1,2]
    assert num.deep_flatten([(1, 2), (3,4)])  == [1, 2, 3, 4]
    assert num.deep_flatten( ((1, 2), (3,4)) )  == [1, 2, 3, 4]
    assert num.deep_flatten([[1,2,3], [4,5,6], [7], [8,9]]) == [1, 2, 3, 4, 5, 6, 7, 8, 9]
    ragged = [[1, 2], [[[[3, 4, 5], 6]]], 7, [8, [9, [10, 11], 12, [13, 14, [15, [[16, 17], 18]]]]]]
    assert num.deep_flatten(ragged) == list(range(1, 19))
    assert num.deep_flatten("hello") == ["hello"]



def test_compare_vectors():
    a = np.array([4.3, 123, 5.2, 532.5, -123.42, 0, 2.3])
    res = num.compare_vectors(a, a)
    assert np.allclose(res, 0.)

    a = np.array([25.])
    b = np.array([26])
    res = num.compare_vectors(a, b)
    assert np.allclose(res, 1.)

    a = np.array([15., 28])
    b = np.array([18, 24.])
    res = num.compare_vectors(a, b)
    assert np.allclose(res, 5.)


    a = np.array([4.3, 123, 5.2, 532.5, -123.42, 0, 2.3])
    b = np.array([52.4, 25.3, -2.4, 5, 5, 9.1, -11])
    res1 = num.compare_vectors(a, b)
    res2 = num.compare_vectors(b, a)
    assert np.allclose(res1, res2)


    a = np.array([15., 28, 4])
    b = np.array([18, 24., 10])
    res = num.compare_vectors(a, b)
    assert np.allclose(res, 7.810249675906654)
    res2 = num.compare_vectors(a, b, trim_edges=1)
    assert np.allclose(res2, 4)



def manual_test_compare_states():     # MANUAL TEST
    state1 = np.array([10, 20, 30])
    state2 = np.array([10.3, 19.9, 30.2])

    print()
    num.compare_states(state1, state2, verbose=True)
    """
    Max of unsigned absolute differences:  0.3000000000000007
    L2 norm of differences (vectors) / Frobenius norm (matrices):  0.3741657386773947
    Relative differences:  [-0.03        0.005      -0.00666667]
    Max of unsigned relative differences:  0.030000000000000072
    Mean of relative differences:  -0.010555555555555547
    Median of relative differences:  -0.006666666666666643
    Standard deviation of relative differences:  0.014550889837454275
    np.allclose with lax tolerance?  (rtol=1e-01, atol=1e-01) :  True
    np.allclose with mid tolerance?  (rtol=1e-02, atol=1e-03) :  False
    np.allclose with tight tolerance?  (rtol=1e-03, atol=1e-05) :  False
    np.allclose with extra-tight tolerance?  (rtol=1e-05, atol=1e-08) :  False
    """


    state1 = np.array([[10, 20, 30],
                       [100, 200, 300]])
    state2 = np.array([[10.3, 19.9, 30.2],
                       [103, 199, 302]])

    print()
    num.compare_states(state1, state2, verbose=True)
    """
    Max of unsigned absolute differences:  3.0
    L2 norm of differences (vectors) / Frobenius norm (matrices):  3.760319135392633
    Relative differences:  [[-0.03        0.005      -0.00666667]
     [-0.03        0.005      -0.00666667]]
    Max of unsigned relative differences:  0.030000000000000072
    Mean of relative differences:  -0.010555555555555552
    Median of relative differences:  -0.006666666666666655
    Standard deviation of relative differences:  0.014550889837454246
    np.allclose with lax tolerance?  (rtol=1e-01, atol=1e-01) :  True
    np.allclose with mid tolerance?  (rtol=1e-02, atol=1e-03) :  False
    np.allclose with tight tolerance?  (rtol=1e-03, atol=1e-05) :  False
    np.allclose with extra-tight tolerance?  (rtol=1e-05, atol=1e-08) :  False
    """


    state1 = np.array([[10, 20, 30],
                       [100, 200, 300]])
    state2 = np.array([[10.0001, 19.9999, 30.0001],
                       [100.0004, 199.9985, 300.0003]])

    print()
    num.compare_states(state1, state2, verbose=True)
    """
    Max of unsigned absolute differences:  0.0014999999999929514
    L2 norm of differences (vectors) / Frobenius norm (matrices):  0.001590597372048279
    Relative differences:  [[-1.00000000e-05  5.00000000e-06 -3.33333333e-06]
     [-4.00000000e-06  7.50000000e-06 -1.00000000e-06]]
    Max of unsigned relative differences:  9.999999999976694e-06
    Mean of relative differences:  -9.722222222130482e-07
    Median of relative differences:  -2.166666666632011e-06
    Standard deviation of relative differences:  5.82651718172416e-06
    np.allclose with lax tolerance?  (rtol=1e-01, atol=1e-01) :  True
    np.allclose with mid tolerance?  (rtol=1e-02, atol=1e-03) :  True
    np.allclose with tight tolerance?  (rtol=1e-03, atol=1e-05) :  True
    np.allclose with extra-tight tolerance?  (rtol=1e-05, atol=1e-08) :  True
    """



def test_gradient_order4_1d():
    result = num.gradient_order4_1d([10, 20, 30, 40, 50, 60], 1)
    assert np.allclose(result, np.full(6, 10., dtype=float))

    result = num.gradient_order4_1d([10, 80, 30, 50, 10, 90], dx=2)
    #print(result)
    assert np.allclose(result, [136.66666667, -24.16666667, -10. , -7.08333333, -17.91666667, 138.75])



def test_gradient_order4():
    result = num.gradient_order4(np.array([10, 20, 30, 40, 50, 60]), 1)
    assert np.allclose(result, np.full(6, 10., dtype=float))

    dx = 2.
    arr = np.array([10, 80, 30, 50, 10, 90])
    result = num.gradient_order4(arr, dx)
    print(result)
    assert np.allclose(result, [35., -25., -10. , -7.08333333, -20., 40.])  # Note how the edges and next-to-edge values
                                                                            # differ from those checked in test_gradient_order4_1d

    # At the boundary values, gradient_order4() performs simple 2-point forward and backward differences (accuracy order 1),
    # just like np.gradient() with the argument edge_order=1
    first_order_boundary = np.gradient(arr, dx, edge_order=1)
    assert np.allclose(result[0] , first_order_boundary[0])
    assert np.allclose(result[-1], first_order_boundary[-1])

    # At points 1 away from the edges, gradient_order4() performs simple 2-point forward and backward differences (accuracy order 1)
    assert np.allclose(result[1] ,  (arr[2]-arr[1])/dx)     # Forward difference quotient, at 2nd point
    assert np.allclose(result[-2] , (arr[-2]-arr[-3])/dx)   # Backward difference quotient, at 2nd-to-last point

    # At the inner points (at least 2 away from the edges), both gradient_order4() and gradient_order4_1d() give the same results
    result_1d = num.gradient_order4_1d(arr, dx=dx)
    assert np.allclose(result[2], result_1d[2])
    assert np.allclose(result[3], result_1d[3])



def test_expand_matrix():
    m = np.array([[1, 2],
                  [3, 4]])
    result = num.expand_matrix_boundary(m)
    assert np.allclose(result ,
                       [[1, 1, 2, 2],
                        [1, 1, 2, 2],
                        [3, 3, 4, 4],
                        [3, 3, 4, 4]])

    m = np.array([[1]])
    result = num.expand_matrix_boundary(m)
    assert np.allclose(result ,
                       [[1, 1, 1],
                        [1, 1, 1],
                        [1, 1, 1]])

    m = np.array([[1, 2, 3]])
    result = num.expand_matrix_boundary(m)
    assert np.allclose(result ,
                       [[1, 1, 2, 3, 3],
                        [1, 1, 2, 3, 3],
                        [1, 1, 2, 3, 3]])

    m = np.array([[1, 2, 3]]).T     # Transpose
    result = num.expand_matrix_boundary(m)
    assert np.allclose(result ,
                       [[1, 1, 1],
                        [1, 1, 1],
                        [2, 2, 2],
                        [3, 3, 3],
                        [3, 3, 3]])

    m = np.array([[1, 2, 3],
                  [4, 5, 6]])
    result = num.expand_matrix_boundary(m)
    assert np.allclose(result ,
                       [[1, 1, 2, 3, 3],
                        [1, 1, 2, 3, 3],
                        [4, 4, 5, 6, 6],
                        [4, 4, 5, 6, 6]
                        ])