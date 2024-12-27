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
