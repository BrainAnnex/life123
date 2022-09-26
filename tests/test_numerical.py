import numpy as np
from modules.numerical.numerical import Numerical as num



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
