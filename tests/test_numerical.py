import numpy as np
from modules.numerical.numerical import Numerical as num


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