# General numerical methods

from typing import Union
import numpy as np
import pandas as pd


class Numerical:
    """
    Assorted, general numerical methods
    """


    @classmethod
    def curve_intersect_interpolate(cls, df: pd.DataFrame, row_index :int, x, var1, var2) -> (float, float):
        """
        Fine-tune the intersection point between 2 Pandas dataframe columns,
        using linear interpolation

        :param df:          A Pandas dataframe with at least 3 columns
        :param row_index:   The index of the Pandas dataframe row
                                with the smallest absolute value of difference in concentrations
        :param x:           The name of the dataframe column with the independent variable
        :param var1:        The name of the dataframe column with the data from the 1st curve
        :param var2:        The name of the dataframe column with the data from the 2nd curve
        :return:            The pair (time of intersection, common value)
        """
        assert len(df) >= 3, \
            f"the Pandas dataframe must have at least 3 rows (it has {len(df)})"

        col_names = list(df.columns)
        assert x in col_names, f"the given Pandas dataframe lacks a column named `{x}`"
        assert var1 in col_names, f"the given Pandas dataframe lacks a column named `{var1}`"
        assert var2 in col_names, f"the given Pandas dataframe lacks a column named `{var2}`"

        row = df.loc[row_index]
        t = row[x]
        val1 = row[var1]
        val2 = row[var2]

        avg_val = (val1+val2) / 2.

        #print("Coarse intersection coordinates: ",  (t, avg_val))


        next_row = df.loc[row_index+1]      # TODO: check whether it exists
        next_t = next_row[x]
        next_val1 = next_row[var1]
        next_val2 = next_row[var2]

        prev_row = df.loc[row_index-1]      # TODO: check whether it exists
        prev_t = prev_row[x]
        prev_val1 = prev_row[var1]
        prev_val2 = prev_row[var2]

        if np.allclose(val1, val2):     # The intersection occurs at the given middle point
            print("val1 is very close to, or equal to, val2")
            if np.allclose(next_val1, next_val2) or np.allclose(prev_val1, prev_val2):
                raise Exception("the intersection isn't well-defined because there's extended overlap")
            else:
                return (t, avg_val)


        later_inversion = False
        earlier_inversion = False

        if (    (val2 > val1) and (next_val2 < next_val1)
                or
                (val2 < val1) and (next_val2 > next_val1)
            ):   # If there's an inversion in height of points on the two curves, at times t vs. next_t
            later_inversion = True


        # No luck finding a curve inversion with the next point; try the previous point instead

        if (    (val2 > val1) and (prev_val2 < prev_val1)
                or
                (val2 < val1) and (prev_val2 > prev_val1)
            ):   # If there's an inversion in height of points on the two curves, at times t vs. prev_t
            earlier_inversion = True


        if later_inversion and earlier_inversion:
            raise Exception("2 intersections detected")

        if later_inversion:
            return cls.segment_intersect(  t=(t, next_t),
                                           y1=(val1, next_val1),
                                           y2=(val2, next_val2))

        if  earlier_inversion:
            return cls.segment_intersect(  t=(prev_t, t),
                                           y1=(prev_val1, val1),
                                           y2=(prev_val2, val2))


        print("Unable to locate a crossover in the curves; "
              "a coarser value is returned as intersection point")

        return (t, avg_val)



    @classmethod
    def segment_intersect(cls, t, y1, y2) -> (float, float):
        """
        Find the intersection of 2 segments in 2D.
        Their respective endpoints share the same x-coordinates: (t_start, t_end)

        Note: for an alternate method, see line_intersect()

        :param t:   Pair with the joint x-coordinates of the 2 start points,
                        and the joint x-coordinates of the 2 end points
        :param y1:  Pair with the y_coordinates of the endpoints of the 1st segment
        :param y2:  Pair with the y_coordinates of the endpoints of the 2nd segment
                        Note: either segment - but not both - may be horizontal

        :return:    A pair with the (x,y) coord of the intersection;
                        if the segments don't meet their specs (which might result in the lack of an intersection),
                        an Exception is raised
        """
        t_start, t_end = t
        y1_start, y1_end = y1
        y2_start, y2_end = y2

        assert t_end > t_start, "t_end must be strictly > t_start"

        if np.allclose(y1_start, y1_end):   # if line1 is horizontal
            assert y2_start != y2_end, "cannot intersect 2 horizontal segments"

        if np.allclose(y2_start, y2_end):   # if line2 is horizontal
            assert y1_start != y1_end, "cannot intersect 2 horizontal segments"

        delta_t = t_end - t_start
        delta_r = y1_end - y1_start
        delta_d = y2_end - y2_start

        '''
        Solving for t_intersect the equation:
        (t_intersect - t_start) * (delta_r/delta_t) + y1_start 
            = (t_intersect - t_start) * (delta_d/delta_t) + y2_start
        '''
        factor = (y1_start - y2_start) / (delta_r - delta_d)
        t_intersect = t_start -  factor * delta_t
        y_intersect = y1_start - factor * delta_r

        return t_intersect, y_intersect



    @classmethod
    def line_intersect(cls, p1, p2, q1, q2) -> Union[tuple, None]:
        """
        Returns the coordinates of the intersection
        between the line passing thru the points p1 and p2,
        and the line passing thru the points q1 and q2.

        Based on https://stackoverflow.com/a/42727584/5478830

        Note: for an alternate method, see segment_intersect()

        :param p1:  Pair of (x,y) coord of a point on the 1st line
        :param p2:  Pair of (x,y) coord of another point on the 1st line
        :param q1:  Pair of (x,y) coord of a point on the 2nd line
        :param q2:  Pair of (x,y) coord of another point on the 2nd line

        :return:    A pair with the (x,y) coord of the intersection, if it exists;
                        or None if it doesn't exist
        """
        s = np.vstack([p1, p2, q1, q2])     # s for "stacked"
                                            #   This will be a matrix :
                                            #   [[p1_x, p1_y]
                                            #    [p2_x, p2_y]
                                            #    [q1_x, q1_y]
                                            #    [q2_x, q2_y]]

        vertical_vector = np.ones((4, 1))   # Four 1's stacked up vertically
        mat = np.hstack((s, vertical_vector)) # This is the earlier "s" matrix, with an extra columns of 1's
        l1 = np.cross(mat[0], mat[1])       # This is a vector representation of the 1st line
        l2 = np.cross(mat[2], mat[3])       # This is a vector representation of the 2nd line
        x, y, z = np.cross(l1, l2)          # Data about the point of intersection

        if np.allclose(z, 0.):              # Lines are parallel, or almost so
            return None

        return (x/z, y/z)



    @classmethod
    def deep_flatten(cls, items) -> list:
        """
        Completely flatten any lists/tuples of lists or tuples
        into a single list
        EXAMPLES:
            [[1,2,3], [4,5,6], [7], [8,9]] becomes [1, 2, 3, 4, 5, 6, 7, 8, 9]
            (1,2)                       becomes [1,2]
            [(1, 2), (3,4)]             becomes [1, 2, 3, 4]
            5                           becomes [5]
            "hello"                     becomes ["hello"]

        :param items:   A list or tuple possibly containing lists or tuples
        :return:        A flat list
        """
        if items == [] or items == ():
            return []
        elif (type(items) is not list) and (type(items) is not tuple):
            return [items]
        else:
            return cls.deep_flatten(items[0]) + cls.deep_flatten(items[1:])



    @classmethod
    def compare_results(cls, res1, res2) -> float:
        """
        Use a Euclidean distance metric to compare 2 sets of results - typically from 2 runs
        of different accuracy.
        Each result value may be a number or a list or tuple, possibly containing lists or tuples,
        of numbers.
        However, after the structure is flattened,
        the 2 result data sets must have the same number of points

        :param res1:
        :param res2:
        :return:    The distance between the two sets of value, based on an L2 (Euclidean) metric
        """
        flat1 =  cls.deep_flatten(res1)
        flat2 =  cls.deep_flatten(res2)
        assert len(flat1) == len(flat2), \
            "compare_results(): the 2 sets of data (after flattening) must have the same length"

        return cls.compare_vectors(v1=np.array(flat1), v2=np.array(flat2))



    @classmethod
    def compare_vectors(cls, v1: np.array, v2: np.array, metric=None, trim_edges=0) -> float:
        """
        Return the Euclidean distance between the two given same-sized Numpy arrays
        (TODO: offer option to use alternate metrics)

        :param v1:
        :param v2:
        :param metric:      NOT YET USED: the default L2 norm is always used
        :param trim_edges:  (OPTIONAL) number of elements to ditch at each end,
                                prior to comparing the vectors (default: 0)
        :return:            The distance between the two vectors based on an L2 (Euclidean) metric
        """
        # TODO: give a helpful error message if v1 or v2 aren't 1-d,
        #       or if they have different dimensions
        if trim_edges != 0:
            v1 = v1[trim_edges : -trim_edges]
            v2 = v2[trim_edges : -trim_edges]

        diff = v1 - v2
        return np.linalg.norm(diff, ord=None)



    @classmethod
    def compare_states(cls, state1: np.array, state2: np.array, verbose=False) -> None:
        """
        Print out various assessments of how similar two same-sized Numpy arrays are.
        Typically, those arrays will be system state snapshots - but could be anything.

        :param state1:
        :param state2:
        :param verbose: If True, additional output is shown
        :return:        None
        """
        diff = state1 - state2
        abs_diff = abs(diff)
        print("Max of unsigned absolute differences: ", np.max(abs_diff))
        print("L2 norm of differences (vectors) / Frobenius norm (matrices): ", np.linalg.norm(diff, ord=None))

        rel_diff = diff / state1

        if verbose:
            print("Relative differences: ", rel_diff)

        print("Max of unsigned relative differences: ", np.max(abs(rel_diff)))

        print("Mean of relative differences: ", np.mean(rel_diff))
        print("Median of relative differences: ", np.median(rel_diff))
        print("Standard deviation of relative differences: ", np.std(rel_diff))

        print("np.allclose with lax tolerance?  (rtol=1e-01, atol=1e-01) : ",
              np.allclose(state1 , state2, rtol=1e-01, atol=1e-01))
        print("np.allclose with mid tolerance?  (rtol=1e-02, atol=1e-03) : ",
              np.allclose(state1 , state2, rtol=1e-02, atol=1e-03))
        print("np.allclose with tight tolerance?  (rtol=1e-03, atol=1e-05) : ",
              np.allclose(state1 , state2, rtol=1e-03, atol=1e-05))
        print("np.allclose with extra-tight tolerance?  (rtol=1e-05, atol=1e-08) : ",
              np.allclose(state1 , state2, rtol=1e-05, atol=1e-08))



    @classmethod
    def gradient_order4_1d(cls, arr, dx=1.0, dtype='float'):
        """
        Compute the gradient, from the values in the given array,
        using the 5-point Central Difference, which produces an accuracy of order 4.

        At the boundary, or close to it, different 5-point stencils are used,
        still resulting in an accuracy of order 4.
        For the coefficients, see: https://web.media.mit.edu/~crtaylor/calculator.html

        Returns the same size as the input array.

        For multi-dimensional cases, use gradient_order4()
        For a simpler, but less accurate (order 2) approach, simply use Numpy gradient()

        :param arr:     One-dimensional Numpy array (or list, or tuple) of numbers,
                        with at least 5 elements
        :param dx:      Delta_x (assumed constant)
        :param dtype:   Data type to use for the elements of the returned array

        :return:        A Numpy array with the same size as the one passed in arr
        """
        arr_size = len(arr)

        assert arr_size >= 5, "gradient_order4_1d(): input numpy array must have at least 5 elements"

        result = np.zeros(arr_size, dtype=dtype)

        # For the leftmost boundary point, use the 5-point forward difference
        i = 0
        result[i] = (-25*arr[i] +48*arr[i+1] -36*arr[i+2] +16*arr[i+3] -3*arr[i+4])/12

        # For the 2nd leftmost boundary point, use the skewed 5-point central differences
        i = 1
        result[i] = (-3*arr[i-1] -10*arr[i] +18*arr[i+1] -6*arr[i+2] +arr[i+3])/12


        # Now process the interior points

        # Coefficients for the 2nd order central finite differences for the first derivative
        C2 = -1/12
        C1 = 2/3

        for i in range(2, arr_size-2):
            result[i] = -C2 * arr[i-2] \
                        -C1 * arr[i-1] \
                        +C1 * arr[i+1] \
                        +C2 * arr[i+2]


        # For the 2nd rightmost boundary point, use the skewed 5-point central differences
        i = arr_size-2
        result[i] = (-arr[i-3] +6*arr[i-2] -18*arr[i-1] +10*arr[i] +3*arr[i+1])/12

        # For the rightmost boundary point, use the 5-point backward difference
        i = arr_size-1
        result[i] = (3*arr[i-4] -16*arr[i-3] +36*arr[i-2] -48*arr[i-1] +25*arr[i])/12

        return result / dx



    @classmethod
    def gradient_order4(cls, f: np.array, *varargs) -> Union[np.array, list]:
        """
        Compute the gradient, from the values in the given (possibly multidimensional) array,
        using the 5-point Central Difference, which produces an accuracy of order 4.

        At the boundary, and at the points immediately adjacent to the boundary,
        simple 2-point forward or backward differences are used (accuracy order 1.)

        It returns a list of ndarrays (or a single ndarray if there is only 1 dimension).
        Each list element has the same shape as the original array,
        and corresponds to the derivatives of the passed array with respect to each dimension.

        For 1-dimensional cases, can also use gradient_order4_1d(), which produces accuracy order 4 at ALL points.

        (ADAPTED FROM https://gist.github.com/deeplycloudy/1b9fa46d5290314d9be02a5156b48741 , which is
        based on 2nd order version from http://projects.scipy.org/scipy/numpy/browser/trunk/numpy/lib/function_base.py

        :param f:       An N-dimensional array giving samples of a scalar function
        :param varargs: 0, 1, or N scalars giving the sample distances in each direction

        :return:        A list of N numpy arrays,
                        each of the same shape as f giving the derivative of f with respect to each dimension
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
            raise SyntaxError(f"gradient_order4(): invalid number of arguments ({n_args})")


        # Use Central Differences on interior, and first differences on endpoints

        out_vals = []       # List of numpy arrays

        # Create slice objects.
        # Initially, they all are "select everything", ie  [:, :, ..., :] across the N dimensions
        # EXAMPLES:
        #           in 1d, the 5 variables below will each be [slice(None, None, None)]
        #           in 2d, each will be [slice(None, None, None), slice(None, None, None)]
        slice0 = [slice(None)]*N
        slice1 = [slice(None)]*N
        slice2 = [slice(None)]*N
        slice3 = [slice(None)]*N
        slice4 = [slice(None)]*N

        out_type = f.dtype.char
        if out_type not in ['f', 'd', 'F', 'D']:   # float or double (real or complex)
            out_type = 'd'      # Use double-precision floats as a default

        for axis in range(N):
            # select out appropriate parts for this dimension
            out = np.zeros(f.shape, dtype=out_type)

            # First, process the Inner Points along this axis (i.e. not the first or last 2)
            slice0[axis] = slice(2, -2)     # Ditch the first 2 and last 2
            slice1[axis] = slice(None, -4)  # Ditch the last 4
            slice2[axis] = slice(1, -3)     # Ditch the 1st and the last 3
            slice3[axis] = slice(3, -1)     # Ditch the first 3 and the last one
            slice4[axis] = slice(4, None)   # Ditch the first 4
            # All the above slice objects will produce a size 4 shorter on the current axis

            # 1D equivalent -- out[2:-2] = (f[:4] - 8*f[1:-3] + 8*f[3:-1] - f[4:])/12.0
            out[slice0] = (f[tuple(slice1)] - 8.0*f[tuple(slice2)] + 8.0*f[tuple(slice3)] - f[tuple(slice4)])/12.0

            # Now process the first 2 points
            slice0[axis] = slice(None, 2)
            slice1[axis] = slice(1, 3)
            slice2[axis] = slice(None, 2)
            # 1D equivalent -- out[0:2] = (f[1:3] - f[0:2])
            out[slice0] = (f[tuple(slice1)] - f[tuple(slice2)])

            # And, finally, the last 3 points
            slice0[axis] = slice(-2, None)
            slice1[axis] = slice(-2, None)
            slice2[axis] = slice(-3, -1)
            # 1D equivalent -- out[-2:] = (f[-2:] - f[-3:-1])
            out[slice0] = (f[tuple(slice1)] - f[tuple(slice2)])


            # Divide by the step size used on this axis
            out_vals.append(out / dx[axis])

            # Reset the slice object in this dimension to ":" (ie, "select everything")
            slice0[axis] = slice(None)
            slice1[axis] = slice(None)
            slice2[axis] = slice(None)
            slice3[axis] = slice(None)
            slice4[axis] = slice(None)
        # END for axis

        if N == 1:
            return out_vals[0]
        else:
            return out_vals



    @classmethod
    def expand_matrix_boundary(cls, m):
        """
        Add a row at the top and at the bottom, and also add a column to the left and to the right,
        repeating the edge values of the matrix
        EXAMPLE:  [[1, 2],
                   [3, 4]]
                 will turn to
                  [[1, 1, 2, 2],
                   [1, 1, 2, 2],
                   [3, 3, 4, 4],
                   [3, 3, 4, 4]]
                 Note how the original matrix is "embedded" in the center of the larger one

        :param m:   A Numpy matrix
        :return:
        """
        # Stack up the first row, the matrix, and the last row
        tall = np.concatenate( ([m[0, :]] , m , [m[-1, :]]) )
        #print(tall)

        # Align sideways the 1st column, the matrix and the last column
        wide = np.concatenate( (tall[:,[0]] , tall , tall[:,[-1]]) , axis=1)
        return wide
