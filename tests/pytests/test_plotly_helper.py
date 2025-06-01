import pytest
import numpy as np
from life123 import PlotlyHelper
import plotly



def test_plot_curves():
    # TODO: expand the testing
    fig = PlotlyHelper.plot_curves(x = np.array([1, 2]),
                                   y = np.array([88, 99])
                                  )

    assert type(fig) == plotly.graph_objs._figure.Figure
    
    #print(fig.data)
    assert type(fig.data) == tuple
    assert len(fig.data) == 1
    data0 = fig.data[0]
    assert type(data0) == plotly.graph_objs._scatter.Scatter
    assert np.allclose(data0["x"], np.array([1, 2]))
    assert np.allclose(data0["y"], np.array([88, 99]))

    #print(fig.layout)
    assert type(fig.layout) == plotly.graph_objs._layout.Layout



def test__optimal_subplot_grid():
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=1, max_n_cols=4) == (1, 1)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=2, max_n_cols=4) == (1, 2)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=3, max_n_cols=4) == (1, 3)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=4, max_n_cols=4) == (1, 4)

    assert PlotlyHelper._optimal_subplot_grid_size(ncells=5, max_n_cols=4) == (2, 3)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=6, max_n_cols=4) == (2, 3)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=7, max_n_cols=4) == (2, 4)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=8, max_n_cols=4) == (2, 4)

    assert PlotlyHelper._optimal_subplot_grid_size(ncells=9, max_n_cols=4) == (3, 3)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=10, max_n_cols=4) == (3, 4)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=11, max_n_cols=4) == (3, 4)
    assert PlotlyHelper._optimal_subplot_grid_size(ncells=12, max_n_cols=4) == (3, 4)

    with pytest.raises(Exception):
        PlotlyHelper._optimal_subplot_grid_size(ncells=0)

    with pytest.raises(Exception):
        PlotlyHelper._optimal_subplot_grid_size(ncells=10, max_n_cols=0)
