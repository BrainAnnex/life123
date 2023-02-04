import pytest
from src.modules.heuristics.heuristics import Heuristics as heur



def test_relative_significance():
    # Assess the relative significance of various quantities
    #       relative to a baseline value of 10
    assert heur.relative_significance(1, 10) == "S"
    assert heur.relative_significance(4.9, 10) == "S"
    assert heur.relative_significance(5.1, 10) == "C"
    assert heur.relative_significance(10, 10) == "C"
    assert heur.relative_significance(19.9, 10) == "C"
    assert heur.relative_significance(20.1, 10) == "L"
    assert heur.relative_significance(137423, 10) == "L"
