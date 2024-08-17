# General numerical methods

from typing import Union, List, Tuple
import numpy as np
import pandas as pd
import math


class Numerical:
    """
    Assorted, general numerical methods
    """


    @classmethod
    def simple_least_square(cls, X :np.ndarray, Y: np.ndarray) -> (float, float):
        """
        Give two numeric 1-D arrays X and Y, of the same size,
        do a least-square fit: Y = a + b * X , for some numbers a and b

        :param X:   A numpy numeric 1-D array
        :param Y:   A numpy numeric 1-D array, with the same number of elements as X
        :return:    A pair of numbers (a, b) that provide the least-square fit for Y = a + b * X
        """
        assert np.ndim(X) == 1, "simple_least_square(): the argument `X` must be a 1-D array"
        assert np.ndim(Y) == 1, "simple_least_square(): the argument `Y` must be a 1-D array"
        assert len(X) == len(Y), \
            "simple_least_square(): the two arguments must be numeric 1-D arrays with the SAME dimension"
        #TODO: more validations

        M = np.vstack([np.ones(len(X)), X]).T
        # M will be an n x 2 matrix , where n is the number of data points in each of the arguments.
        # The 1st column is all 1's ; the 2nd column contains the values of X
        # EXAMPLE  of M, if X = array([40., 35., 19.]):
        '''
        array([[ 1.        , 40.],
               [ 1.        , 35.],
               [ 1.        , 19.]])
        '''
        a, b = np.linalg.lstsq(M, Y, rcond=None)[0]     # Carry out the least-square fit  as: Y = a + b X
        return a, b


    @classmethod
    def two_vector_least_square(cls, V :np.ndarray, W :np.ndarray, Y :np.ndarray) -> (float, float):
        """
        Give three numeric 1-D arrays V, W and Y, of the same size,
        do a least-square fit: Y = a * V + b * W, for some numbers a and b

        :param V:
        :param W:
        :param Y:
        :return:
        """
        M = np.vstack([V, W]).T
        a, b = np.linalg.lstsq(M, Y, rcond=None)[0]
        return a, b



    @classmethod
    def reach_threshold(cls, df, x :str, y :str, y_threshold) -> Union[float, None]:
        """
        Given a set of 2-D points, whose x- and y- coordinates are stored, respectively, in the
        columns whose names are given by the arguments x and y, locate the x-coordinate at which
        the curve first reaches the given y-value threshold.
        That is, its intersection with the horizontal line y = y_threshold
        Linear interpolation is used.

        The dataframe is *assumed* sorted by the x variable, increasing monotonically; an Exception is
        raised if that's not the case.
        
        Note that if there are multiple intersections, only the leftmost (smallest x) is returned.

        In case no intersection is found, None is returned.

        :param df:          A Pandas dataframe with at least 2 columns
        :param x:           The name of the column with the x-values
        :param y:           The name of the column with the y-values
        :param y_threshold: A number with a y-value being sought in the xy curve
        :return:            The x-value at which the linearly-interpolated curve first reaches
                                the y-value y_threshold, or None if not found
        """

        # Iterate through the DataFrame rows to find the (first) segment where y0 lies between
        for i in range(len(df) - 1):
            # For each row, except the last one

            x1 = df.loc[i, x]
            x2 = df.loc[i + 1, x]
            assert x2 > x1, \
                f"reach_threshold(): the x-values (column `{x}`) are not monotonically-increasing, as expected"

            y1 = df.loc[i, y]                   # y-value at the start of the sliding interval
            y2 = df.loc[i + 1, y]               # y-value at the end of the sliding interval
            if min(y1, y2) <= y_threshold <= max(y1, y2):    # This checks if y0 is between y1 and y2

                if np.allclose(y1, y2):
                    return x1                   # A convention we're following : "when is the target value first reached?"

                # Perform linear interpolation
                x_intersection = x1 + (y_threshold - y1) * (x2 - x1) / (y2 - y1)
                return x_intersection

        return None  # Return None if no intersection is found



    @classmethod
    def curve_intersect(cls, df: pd.DataFrame,
                        x :str, var1 :str, var2 :str, explain=False) -> (float, float):
        """
        Determine the intersection point between 2 Pandas dataframe columns,
        using linear interpolation.
        If more than one is present, only the first (smallest value of x-coord) one will be returned;
        so, the specified time interval should be narrow enough to bracket the intersection of interest

        :param df:      A Pandas dataframe with at least 3 columns
        :param x:       The name of the dataframe column with the independent variable
        :param var1:    The name of the dataframe column with the data from the 1st curve
        :param var2:    The name of the dataframe column with the data from the 2nd curve
        :param explain: [OPTIONAL] If True, print out some details of the computation
        :return:        The pair (x-coordinate of intersection, common value)

        :return:        The pair (x-coordinate of intersection, common value)
        """
        df_size = len(df)   # Number of rows in the dataframe

        # Validation
        assert df_size >= 2, \
            f"curve_intersect(): the Pandas dataframe must have at least 2 rows (it has {df_size}). Not enough data..."

        col_names = list(df.columns)
        assert x in col_names, f"curve_intersect(): the given Pandas dataframe lacks a column named `{x}`"
        assert var1 in col_names, f"curve_intersect(): the given Pandas dataframe lacks a column named `{var1}`"
        assert var2 in col_names, f"curve_intersect(): the given Pandas dataframe lacks a column named `{var2}`"

        # The next 2 integers are the indices of the boundaries of the passed dataframe
        first_index = df.index[0]
        last_index = df.index[-1]
        if explain:
            print(f"first_index: {first_index}, last_index: {last_index}")

        # Examine the y-values of the two curves at the first (zero-th) point
        prev_val1 = df.loc[first_index, var1]
        prev_val2 = df.loc[first_index, var2]

        if np.allclose(prev_val1, prev_val2):     # Trivial intersection immediately found
            if explain:
                print(f"curve_intersect(): Curve intersection found at data row {first_index}, when `{x}` = {df.loc[first_index, x]}")
                if np.allclose(df.loc[1, var1], df.loc[1, var2]):
                    print(f"curve_intersect(): NOTICE - the intersection isn't well-defined because there's extended overlap "
                          f"of the two curves  between `{x}` = {df.loc[first_index, x]} and `{x}` = {df.loc[1, x]}")

            return (df.loc[first_index, x], prev_val1)


        # Make a note of which curve has a bigger value at the initial point
        if prev_val1 > prev_val2:
            initially_larger = 1
        else:
            initially_larger = 2

        row_index_curve_crossover = None

        # Iterate through the DataFrame rows to find the first occurrence of a curve crossover
        # (inversion about which one is larger)
        for i in range(first_index+1, last_index+1):
            # For each row, except the first one

            # y-coords of the two curves at the i-th point
            prev_val1 = df.loc[i, var1]
            prev_val2 = df.loc[i, var2]

            #print(f"val1: {val1} , val2: {val2}")

            if np.allclose(prev_val1, prev_val2):  # Trivial intersection found at the i-th point
                if explain:
                    print(f"curve_intersect(): Curve intersection found at data row {i}, when `{x}` = {df.loc[i, x]}")
                    if (i+1 < df_size) and np.allclose(df.loc[i+1, var1], df.loc[i+1, var2]):
                        print(f"curve_intersect(): NOTICE - the intersection isn't well-defined because there's extended overlap "
                              f"of the two curves between `{x}` = {df.loc[i, x]} and `{x}` = {df.loc[i+1, x]}")

                return (df.loc[i, x], df.loc[i, var1])


            if (prev_val1 > prev_val2) and (initially_larger == 2):
                # Curve 1 is now the larger one
                row_index_curve_crossover = i       # Value inversion found
                break

            if (prev_val2 > prev_val1) and (initially_larger == 1):
                # Curve 2 is now the larger one
                row_index_curve_crossover = i       # Value inversion found
                break
        # END for


        if row_index_curve_crossover is None:      # No crossover found between the two curves
            if explain:
                print("curve_intersect(): No curve crossover detected; i.e., one curve is always above the other one")
            return None         # No intersection located

        if explain:
            print(f"curve_intersect(): Curve crossover detected at data row: {row_index_curve_crossover}")

        before_index = row_index_curve_crossover - 1    # Dataframe row index just before the crossover
        after_index = row_index_curve_crossover         # Dataframe row index just after the crossover


        row_before = df.loc[before_index]
        prev_t = row_before[x]
        prev_val1 = row_before[var1]
        prev_val2 = row_before[var2]

        row_after = df.loc[after_index]
        next_t = row_after[x]
        next_val1 = row_after[var1]
        next_val2 = row_after[var2]

        if np.allclose(prev_t, next_t):     # Curve discontinuity that requires special handling
            print(f"WARNING: curve_intersect(): Discontinuity detected at `{x}` = {prev_t}")
            if np.allclose(prev_val1, next_val1):
                return (prev_t, prev_val1)
            if np.allclose(prev_val2, next_val2):
                return (prev_t, prev_val2)
            return None


        return cls.segment_intersect(  t=(prev_t, next_t),
                                       y1=(prev_val1, next_val1),
                                       y2=(prev_val2, next_val2)
                                    )



    @classmethod
    def segment_intersect(cls, t, y1, y2) -> (float, float):
        """
        Find the intersection of two segments in 2D.
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

        assert t_end > t_start, \
            f"segment_intersect(): t_end must be strictly > t_start.  Here, t_start = {t_start} , t_end = {t_end}"

        if np.allclose(y1_start, y1_end):       # If line1 is horizontal
            assert y2_start != y2_end, "segment_intersect(): cannot intersect 2 horizontal segments"

        if np.allclose(y2_start, y2_end):       # If line2 is horizontal
            assert y1_start != y1_end, "segment_intersect(): cannot intersect 2 horizontal segments"

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
    def gradient_order4_1d(cls, arr :Union[np.array, List, Tuple], dx=1.0, dtype='float') -> np.array:
        """
        Compute the gradient, from UNIT-SPACED x-values and y-values in the given 1-dimensional array,
        using the 5-point Central Difference, which produces an accuracy of order 4.

        At the boundary, or close to it, different 5-point stencils are used,
        still resulting in an accuracy of order 4.
        For the coefficients, see: https://web.media.mit.edu/~crtaylor/calculator.html

        Returns the same size as the input array.

        Note:
            * For multi-dimensional cases, use gradient_order4()
            * For a simpler, but less accurate (order 2) approach, use Numpy gradient(),
              which is also usable in case of uneven grids of x values

        :param arr:     One-dimensional Numpy array (or list, or tuple) of numbers,
                            with at least 5 elements
        :param dx:      Delta_x (assumed constant)
        :param dtype:   Data type to use for the elements of the returned array

        :return:        A Numpy array with the same size as the one passed in arr
        """
        arr_size = len(arr)

        assert arr_size >= 5, "gradient_order4_1d(): argument `arr` must have at least 5 elements"

        result = np.zeros(arr_size, dtype=dtype)    # The result is the same size as the input array


        # First, start with the 2 leftmost boundary points

        # For the leftmost boundary point, use the 5-point forward difference
        i = 0
        result[i] = (-25*arr[i] +48*arr[i+1] -36*arr[i+2] +16*arr[i+3] -3*arr[i+4])/12

        # For the 2nd leftmost boundary point, use the skewed 5-point central differences
        i = 1
        result[i] = (-3*arr[i-1] -10*arr[i] +18*arr[i+1] -6*arr[i+2] +arr[i+3])/12


        # Now, process the interior points

        # Coefficients for the 2nd order central finite differences for the first derivative
        C2 = -1/12
        C1 = 2/3

        for i in range(2, arr_size-2):
            result[i] = -C2 * arr[i-2] \
                        -C1 * arr[i-1] \
                        +C1 * arr[i+1] \
                        +C2 * arr[i+2]


        # Finally, conclude with the 2 rightmost boundary points

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

        All points need to be equally-spaced along each dimension.

        At the boundary, and at the points immediately adjacent to the boundary,
        the gradient is computed by the very simple 2-point forward or backward differences (accuracy order 1.)

        It returns a list of ndarrays (or a single ndarray if there is only 1 dimension).
        Each list element has the same shape as the original array,
        and corresponds to the derivatives of the passed array with respect to each dimension.

        For 1-dimensional cases, can also use gradient_order4_1d(), which produces accuracy order 4 at ALL points.

        (ADAPTED FROM https://gist.github.com/deeplycloudy/1b9fa46d5290314d9be02a5156b48741 ,
        which is based on 2nd order version from
        http://projects.scipy.org/scipy/numpy/browser/trunk/numpy/lib/function_base.py

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
    def gradient_uneven_grid(cls, x_values :np.array, f :np.array, stencil :int) -> np.ndarray:
        """
        Compute and return the gradient of a 1-D function f(x) at the given grid points,
        which may be uneven.  An arbitrary stencil ("sliding window") size may be used.

        For gradients (first derivatives), the accuracy will be dependent on the size of
        the stencil - though overly large stencils might lead to oscillations.

        Based on B. Fornberg, "Calculation of Weights in Finite Difference Formulas", 1998
        (https://epubs.siam.org/doi/abs/10.1137/S0036144596322507)

        :param x_values:    A Numpy array of x values (independent variable) that
                                may be unevenly-spaced
        :param f:           A Numpy array of the values of the function at the above grid point
        :param stencil:     An integer between 2 and the number of grid points.
        :return:            A Numpy array
        """
        # TODO: allow the gradients to also be computed at x values not on the grid
        assert type(stencil) == int, \
            "gradient_uneven_grid(): value for stencil must be an integer"

        assert stencil >= 2, \
            "gradient_uneven_grid(): value for stencil must be an integer >= 2"

        assert type(x_values) == np.ndarray, \
            f"gradient_uneven_grid(): argument `x_values` must be a Numpy array; instead, it is {type(x_values)}"
        assert type(f) == np.ndarray, \
            f"gradient_uneven_grid(): argument `f` must be a Numpy array; instead, it is {type(f)}"

        assert np.ndim(x_values) == 1, \
            "gradient_uneven_grid(): x_values must be a 1-dimensional Numpy array"
        assert np.ndim(f) == 1, \
            "gradient_uneven_grid(): f must be a 1-dimensional Numpy array"

        size = len(x_values)

        assert size == len(f), \
            "gradient_uneven_grid(): x_values and f must have the same dimension"

        assert stencil <= size, \
            f"gradient_uneven_grid(): value for stencil ({stencil}) cannot exceed the number of grid points ({size})"


        gradient_values = np.zeros(size, dtype=np.float64)

        # Start with the leftmost sliding window, of size stencil, across the array x_values
        # Extract the gradients at point in the left half of that window (rounded up if odd-sized)
        window_x_arr = x_values[ : stencil]   # Array of size equal to stencil value
        assert len(window_x_arr) == stencil
        window_f_arr = f[ : stencil]
        #print("- WINDOW: ", window_x_arr)
        grid_index = 0
        for grid_index in range(math.ceil(stencil / 2)):    # Ceiling of division)
            z = x_values[grid_index]
            #print(f"i : {grid_index} | z : {z}")
            numeric_1st_deriv = cls._compute_derivative(z, window_x_arr, window_f_arr)
            gradient_values[grid_index] = numeric_1st_deriv


        # Advance (by steps of 1) the sliding window across the array x_values, as far as it can go
        for window_start_pos in range(1, size - stencil + 1):
            window_x_arr = x_values[window_start_pos : window_start_pos+stencil]  # Array of size equal to stencil value
            window_f_arr = f[window_start_pos : window_start_pos+stencil]
            assert len(window_x_arr) == stencil
            #print("- WINDOW: ", window_x_arr)
            # Also advance by 1 the grid point at which we compute the gradient
            grid_index += 1
            z = x_values[grid_index]
            #print(f"i : {grid_index} | z : {z}")
            numeric_1st_deriv = cls._compute_derivative(z, window_x_arr, window_f_arr)
            gradient_values[grid_index] = numeric_1st_deriv


        # Process the last remaining grid points, using the last stencil
        # (which is at the end of the x_values array)
        grid_index += 1
        while grid_index < size:    # This loop stops at index size-1
            z = x_values[grid_index]
            #print(f"i : {grid_index} | z : {z}")
            numeric_1st_deriv = cls._compute_derivative(z, window_x_arr, window_f_arr)
            gradient_values[grid_index] = numeric_1st_deriv
            grid_index += 1

        return gradient_values



    @classmethod
    def _compute_derivative(cls, z, window_x_arr, window_f_arr, order=1) -> np.float64:
        """
        Compute and return the specified derivative order (by default the gradient)
        of a 1-D function f(x) at the  point x=z.
        The entire grid of values is used as a stencil.

        :param z:               x value (in grid or not) at which to compute the gradient of f
        :param window_x_arr:    Grid point locations
        :param window_f_arr:    Values of the function at the above grid points
        :param order:           Order of the derivative (by default 1st derivative, i.e. gradient)
        :return:                A floating-point value
        """
        c = cls._finite_diff_weights(z=z, x=window_x_arr, m=order)
        col = c[:, 1]       # Extract 2nd column (k = 1, to be used for 1st derivatives)
        #print("col: ", col)
        #print("window_f_arr: ", window_f_arr)
        numeric_1st_deriv = np.dot(col, window_f_arr)
        return numeric_1st_deriv



    @classmethod
    def _finite_diff_weights(cls, z, x, m) -> np.ndarray:
        """
        Based on B. Fornberg, "Calculation of Weights in Finite Difference Formulas", 1998
        (https://epubs.siam.org/doi/abs/10.1137/S0036144596322507)

        :param z:   x value (in grid or not) where approximations are to be accurate
        :param x:   Grid point locations, x0 thru xn : (n+1) of them
        :param m:   Highest derivative for which weights are sought.
                        EXAMPLE : m=1 to compute up to 1st derivative (gradient)

        :return:    (n+1) x (m+1) Numpy array:
                        c(i, j) stores the weights at grid locations x(i) for derivatives of order j,
                        with 0 <= i <= n , and 0 <= j < m, and n being one less than total number of grid points
        """
        n = len(x) - 1      # One less than total number of grid points

        c = np.zeros((n+1, m+1), dtype=np.float64)

        c1 = 1.0
        c4 = x[0] - z

        c[0, 0] = 1.0

        for i in range(1, n+1):
            mn = min(i, m)
            c2 = 1.0
            c5 = c4
            c4 = x[i] - z

            for j in range(i):
                c3 = x[i] - x[j]
                c2 *= c3

                if j == i - 1:
                    for k in range(mn, 0, -1):
                        c[i, k] = c1 * (k * c[i-1, k-1] - c5 * c[i-1, k]) / c2

                    c[i, 0] = -c1 * c5 * c[i-1, 0] / c2

                for k in range(mn, 0, -1):
                    c[j, k] = (c4 * c[j, k] - k * c[j, k-1]) / c3

                c[j, 0] = c4 * c[j, 0] / c3
            # END for j in range(i)

            c1 = c2

        return c



    @classmethod
    def expand_matrix_boundary(cls, m) -> np.ndarray:
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
        :return:    A Numpy matrix, with 2 extra rows (at top and bottom) and 2 extra columns (at left and right)
        """
        # Stack up the first row, the matrix, and the last row
        tall = np.concatenate( ([m[0, :]] , m , [m[-1, :]]) )
        #print(tall)

        # Align sideways the 1st column, the matrix and the last column
        wide = np.concatenate( (tall[:,[0]] , tall , tall[:,[-1]]) , axis=1)
        return wide
