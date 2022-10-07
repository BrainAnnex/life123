import pytest
import numpy as np
from modules.movies.movies import Movie, MovieArray


###  For class MovieArray  ###

# TODO



###  For class MovieArray  ###

def test_store():
    m = MovieArray()

    m.store(pars=1, data_snapshot=np.array([1., 2., 3.]), caption="first entry")
    assert m.parameters == [1]
    assert m.captions == ["first entry"]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.arr, [[1., 2., 3.]])

    m.store(pars=2, data_snapshot=np.array([10., 11., 12.]), caption="second entry")

    assert m.parameters == [1, 2]
    assert m.captions == ["first entry", "second entry"]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.arr, [[1., 2., 3.],
                               [10., 11., 12.]])

    m.store(pars=3, data_snapshot=np.array([-10., -11., -12.]))
    assert m.parameters == [1, 2, 3]
    assert m.captions == ["first entry", "second entry", ""]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.arr, [[1., 2., 3.],
                               [10., 11., 12.],
                               [-10., -11., -12.]])

    with pytest.raises(Exception):
        m.store(pars=666, data_snapshot=np.array([1., 2., 3., 4., 5.]),
                caption="doesn't conform to shape of earlier entries!")


    # Again, in higher dimensions
    m_2D = MovieArray()

    m_2D.store(pars={"a": 4., "b": 12.3},
               data_snapshot=np.array([[1., 2., 3.],
                                       [10., 11., 12.]]))
    assert m_2D.parameters == [{"a": 4., "b": 12.3}]
    assert m_2D.captions == [""]
    assert m_2D.snapshot_shape == (2, 3)
    assert np.allclose(m_2D.arr, [[1., 2., 3.],
                                  [10., 11., 12.]])

    m_2D.store(pars={"a": 400., "b": 123},
               data_snapshot=np.array([[-1., -2., -3.],
                                       [-10., -11., -12.]]
                                      ),
               caption="2nd matrix")
    assert m_2D.parameters == [{"a": 4., "b": 12.3}, {"a": 400., "b": 123}]
    assert m_2D.captions == ["", "2nd matrix"]
    assert m_2D.snapshot_shape == (2, 3)
    expected = np.array([
                            [[1., 2., 3.],
                             [10., 11., 12.]]
                            ,
                            [[-1., -2., -3.],
                             [-10., -11., -12.]]
                        ])

    assert np.allclose(m_2D.arr, expected)

    with pytest.raises(Exception):
        m_2D.store(pars=666, data_snapshot=np.array([10., 11., 12.]), caption="wrong shape!")