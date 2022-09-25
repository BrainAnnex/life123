# General numerical methods

import numpy as np


class Numerical:

    @classmethod
    def foo(cls, x):
        return x*x


    @classmethod
    def gradient_order4_1d(cls, arr, dx: float, dtype='float'):
        """
        Compute the gradient, from the values in the given array,
        using the 5-point Central Difference, which produces an accuracy of order 4.

        At the boundary, or close to it, different 5-point stencils are used,
        still resulting in an accuracy of order 4.

        Returns the same shape as the input array.

        For multi-dimensional cases, use gradient_order4()

        :param arr:     One-dimensional array of numbers
        :param dx:      Delta_x (assumed constant)
        :param dtype:   Data type to use for the elements of the returned array
        :return:
        """
        arr_size = len(arr)

        result = np.zeros(arr_size, dtype=dtype)

        # For the leftmost boundary point, use the 5-point forward difference
        i = 0
        result[i] =(-25*arr[i] +48*arr[i+1] -36*arr[i+2] +16*arr[i+3] -3*arr[i+4])/12

        # For the 2nd leftmost boundary point, use skewed 5-point central differences
        i = 1
        result[i] = (-3*arr[i-1] -10*arr[i] +18*arr[i+1] -6*arr[i+2] +arr[i+3])/12


        # Now do the interior points

        # Coefficients for the 2nd order central finite differences for the first derivative
        C2 = -1/12
        C1 = 2/3

        for i in range(2, arr_size-2):
            result[i] = -C2 * arr[i-2] \
                        -C1 * arr[i-1] \
                        +C1 * arr[i+1] \
                        +C2 * arr[i+2]


        # For the 2nd rightmost boundary point, use skewed 5-point central differences
        i = arr_size-2
        result[i] = (-arr[i-3] +6*arr[i-2] -18*arr[i-1] +10*arr[i] +3*arr[i+1])/12

        # For the rightmost boundary point, use the 5-point backward difference
        i = arr_size-1
        result[i] =(3*arr[i-4] -16*arr[i-3] +36*arr[i-2] -48*arr[i-1] +25*arr[i])/12

        return result / dx



    @classmethod
    def gradient_order4(cls, f, *varargs):
        """
        Compute the gradient, from the values in the given (possibly multidimensional) array,
        using the 5-point Central Difference, which produces an accuracy of order 4.

        At the boundary, simple 2-point forward or backward differences are used (accuracy order 1.)
        TODO: clarify what's happening at the points immediately next to the boundary.

        Returns the same shape as the input array.

        For 1-dimensional cases, can also use gradient_order4_1d(), which produces accuracy order 4 at ALL points.

        (ADAPTED FROM https://gist.github.com/deeplycloudy/1b9fa46d5290314d9be02a5156b48741 , which is
        based on 2nd order version from http://projects.scipy.org/scipy/numpy/browser/trunk/numpy/lib/function_base.py

        :param f:       An N-dimensional array giving samples of a scalar function
        :param varargs: 0, 1, or N scalars giving the sample distances in each direction
        :return:        N arrays of the same shape as f giving the derivative of f with respect to each dimension
        """

        N = len(f.shape)            # Number of dimensions
        n_args = len(varargs)       # Number of arguments passed after the first one (f)
        if n_args == 0:
            dx = [1.0]*N            # List with the element 1.0 repeated N times
        elif n_args == 1:
            dx = [varargs[0]]*N
        elif n_args == N:
            dx = list(varargs)
        else:
            raise SyntaxError(f"invalid number of arguments (n_args)")


        # Use Central Differences on interior, and first differences on endpoints

        outvals = []

        # create slice objects --- initially all are [:, :, ..., :]
        slice0 = [slice(None)]*N
        slice1 = [slice(None)]*N
        slice2 = [slice(None)]*N
        slice3 = [slice(None)]*N
        slice4 = [slice(None)]*N

        out_type = f.dtype.char
        if out_type not in ['f', 'd', 'F', 'D']:   # float or double (real or complex)
            out_type = 'f'      # Use single-precision floats as a default

        for axis in range(N):
            # select out appropriate parts for this dimension
            out = np.zeros(f.shape, dtype=out_type)

            slice0[axis] = slice(2, -2)
            slice1[axis] = slice(None, -4)
            slice2[axis] = slice(1, -3)
            slice3[axis] = slice(3, -1)
            slice4[axis] = slice(4, None)
            # 1D equivalent -- out[2:-2] = (f[:4] - 8*f[1:-3] + 8*f[3:-1] - f[4:])/12.0
            out[slice0] = (f[slice1] - 8.0*f[slice2] + 8.0*f[slice3] - f[slice4])/12.0

            slice0[axis] = slice(None, 2)
            slice1[axis] = slice(1, 3)
            slice2[axis] = slice(None, 2)
            # 1D equivalent -- out[0:2] = (f[1:3] - f[0:2])
            out[slice0] = (f[slice1] - f[slice2])

            slice0[axis] = slice(-2, None)
            slice1[axis] = slice(-2, None)
            slice2[axis] = slice(-3, -1)
            ## 1D equivalent -- out[-2:] = (f[-2:] - f[-3:-1])
            out[slice0] = (f[slice1] - f[slice2])


            # divide by step size
            outvals.append(out / dx[axis])

            # reset the slice object in this dimension to ":"
            slice0[axis] = slice(None)
            slice1[axis] = slice(None)
            slice2[axis] = slice(None)
            slice3[axis] = slice(None)
            slice4[axis] = slice(None)

        if N == 1:
            return outvals[0]
        else:
            return outvals
