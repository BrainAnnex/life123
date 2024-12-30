import pytest
from life123 import PlotlyHelper



def test_get_default_colors():
    assert PlotlyHelper.get_default_colors(1) == ['darkturquoise']
    assert PlotlyHelper.get_default_colors(2) == ['darkturquoise', 'green']



def test_lighten_color():
    assert PlotlyHelper.lighten_color("yellow", factor=1) == "rgb(255,255,255)"
    assert PlotlyHelper.lighten_color("green", factor=1) == "rgb(255,255,255)"
    assert PlotlyHelper.lighten_color("black", factor=1) == "rgb(255,255,255)"

    assert PlotlyHelper.lighten_color("red", factor=0) == "rgb(255,0,0)"
    # Notice that (0,255,0) is the CSS color "lime", rather than "green"!
    assert PlotlyHelper.lighten_color("lime", factor=0) == "rgb(0,255,0)"
    assert PlotlyHelper.lighten_color("green", factor=0) == "rgb(0,128,0)"
    assert PlotlyHelper.lighten_color("blue", factor=0) == "rgb(0,0,255)"

    assert PlotlyHelper.lighten_color("blue", factor=.8) == "rgb(204,204,254)"



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
