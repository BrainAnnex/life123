import numpy as np
from life123.visualization.colors import Colors


def test_colors_class():
    """
    print()
    print(Colors.BASE_COLORS)
    print(Colors.TABLEAU_COLORS)
    print(Colors.XKCD_COLORS)
    print(Colors.CSS4_COLORS)
    print(Colors.get_named_colors_mapping())
    """
    assert len(Colors.get_named_colors_mapping()) == 1155



def test_get_named_colors_mapping():
    assert Colors.get_named_colors_mapping()["yellow"] == "#FFFF00"



def test_to_rgb():
    assert Colors.to_rgb("red") == (1.0, 0.0, 0.0)
    assert Colors.to_rgb("lime") == (0.0, 1.0, 0.0)
    assert Colors.to_rgb("blue") == (0.0, 0.0, 1.0)
    assert Colors.to_rgb("yellow") == (1.0, 1.0, 0.0)
    assert np.allclose(Colors.to_rgb("turquoise"),
                       (0.25098039215686274, 0.8784313725490196, 0.8156862745098039))
    assert np.allclose(Colors.to_rgb("lavender"),
                       (0.9019607843137255, 0.9019607843137255, 0.9803921568627451))
    assert np.allclose(Colors.to_rgb("green"),
                       (0.0, 0.5019607843137255, 0.0))



def test_lighten_color():
    assert Colors.lighten_color("yellow", factor=1) == "rgb(255,255,255)"
    assert Colors.lighten_color("green", factor=1) == "rgb(255,255,255)"
    assert Colors.lighten_color("black", factor=1) == "rgb(255,255,255)"

    assert Colors.lighten_color("red", factor=0) == "rgb(255,0,0)"
    # Notice that (0,255,0) is the CSS color "lime", rather than "green"!
    assert Colors.lighten_color("lime", factor=0) == "rgb(0,255,0)"
    assert Colors.lighten_color("green", factor=0) == "rgb(0,128,0)"
    assert Colors.lighten_color("blue", factor=0) == "rgb(0,0,255)"

    assert Colors.lighten_color("blue", factor=.8) == "rgb(204,204,254)"



def test_assign_default_colors():
    assert Colors.assign_default_colors(1) == ['darkturquoise']
    assert Colors.assign_default_colors(2) == ['darkturquoise', 'green']
