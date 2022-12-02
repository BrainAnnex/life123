import pytest
import numpy as np
from modules.movies.movies import Movie, MovieArray


###  For class MovieGeneral  ###

# TODO



###  For class MovieTabular  ###

# TODO



###  For class MovieArray  ###

def test_MovieArray():
    m = MovieArray()

    m.store(pars=10, data_snapshot=np.array([1., 2., 3.]), caption="first entry")
    assert m.parameters == [10]
    assert m.captions == ["first entry"]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.movie, [[1., 2., 3.]])
    assert len(m) == 1
    assert str(m) == "MovieArray object with 1 snapshot(s) parametrized by `SYSTEM TIME`"

    m.store(pars=20, data_snapshot=np.array([10., 11., 12.]), caption="second entry")
    assert m.parameters == [10, 20]
    assert m.captions == ["first entry", "second entry"]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.movie, [[1., 2., 3.],
                                 [10., 11., 12.]])
    assert len(m) == 2
    assert str(m) == "MovieArray object with 2 snapshot(s) parametrized by `SYSTEM TIME`"

    m.store(pars=30, data_snapshot=np.array([-10., -11., -12.]))
    assert m.parameters == [10, 20, 30]
    assert m.captions == ["first entry", "second entry", ""]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.movie, [[1., 2., 3.],
                                 [10., 11., 12.],
                                 [-10., -11., -12.]])
    assert len(m) == 3
    assert str(m) == "MovieArray object with 3 snapshot(s) parametrized by `SYSTEM TIME`"

    with pytest.raises(Exception):
        m.store(pars=666, data_snapshot=np.array([1., 2., 3., 4., 5.]),
                caption="doesn't conform to shape of earlier entries!")


    # Again, in higher dimensions (and using methods to fetch the attributes)
    m_2D = MovieArray(parameter_name="a,b values")

    m_2D.store(pars={"a": 4., "b": 12.3},
               data_snapshot=np.array([[1., 2., 3.],
                                       [10., 11., 12.]]))
    assert m_2D.get_parameters() == [{"a": 4., "b": 12.3}]
    assert m_2D.get_captions() == [""]
    assert m_2D.get_shape() == (2, 3)
    assert np.allclose(m_2D.movie, [[1., 2., 3.],
                                    [10., 11., 12.]])
    assert len(m_2D) == 1
    assert str(m_2D) == "MovieArray object with 1 snapshot(s) parametrized by `a,b values`"

    m_2D.store(pars={"a": 400., "b": 123},
               data_snapshot=np.array([[-1., -2., -3.],
                                       [-10., -11., -12.]]
                                      ),
               caption="2nd matrix")
    assert m_2D.get_parameters() == [{"a": 4., "b": 12.3}, {"a": 400., "b": 123}]
    assert m_2D.get_captions() == ["", "2nd matrix"]
    assert m_2D.get_shape() == (2, 3)
    expected = np.array([
                            [[1., 2., 3.],
                             [10., 11., 12.]]
                            ,
                            [[-1., -2., -3.],
                             [-10., -11., -12.]]
                        ])
    assert np.allclose(m_2D.movie, expected)
    assert len(m_2D) == 2
    assert str(m_2D) == "MovieArray object with 2 snapshot(s) parametrized by `a,b values`"

    with pytest.raises(Exception):
        m_2D.store(pars=666, data_snapshot=np.array([10., 11., 12.]), caption="wrong shape for the data!")