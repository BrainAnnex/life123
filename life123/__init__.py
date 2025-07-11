"""
    life123
    ~~~~~~~

    Dynamical Modeling of Biological Systems

	https://Life123.science

    :copyright:     (c) 2022-2025 by Julian West and the Life123 project.
    :license:       MIT.  See LICENSE file for more details.
"""

__version__ = "1.0.0rc4"


from life123.bio_sim_1d import BioSim1D
from life123.bio_sim_2d import BioSim2D
from life123.bio_sim_3d import BioSim3D
from life123.chem_data import ChemData
from life123.html_log import HtmlLog
from life123.collections import (
    CollectionTabular,
    CollectionArray,
    Collection
)
from life123.history import (
    HistoryBinConcentration,
    HistoryUniformConcentration,
    HistoryReactionRate
)
from life123.numerical import Numerical
from life123.reactions import ReactionGeneric, ReactionEnz, Reactions
from life123.thermodynamics import ThermoDynamics
from life123.uniform_compartment import UniformCompartment
from life123.reaction_kinetics import (ReactionKinetics, VariableTimeSteps)

from life123.visualization.graphic_log import GraphicLog
from life123.visualization.plotly_helper import PlotlyHelper
from life123.visualization.py_graph_visual import PyGraphVisual
from life123.visualization.colors import Colors


__all__ = [
    'BioSim1D',
    'BioSim2D',
    'BioSim3D',
    'ChemData',
    'HtmlLog',
    'CollectionTabular',
    'CollectionArray',
    'Collection',
    'Colors',
    'HistoryBinConcentration',
    'HistoryUniformConcentration',
    'HistoryReactionRate',
    'Numerical',
    'ReactionGeneric',
    'Reactions',
    'ReactionEnz',
    'ReactionKinetics',
    'ThermoDynamics',
    'VariableTimeSteps',
    'UniformCompartment',
    'GraphicLog',
    'PlotlyHelper',
    'PyGraphVisual'
]



def version():
    return __version__


def check_version(expected :str, enforce=False) -> None:
    """
    Check the passed version number against the actual version number of this library

    :param expected:A string with the expected version number
    :param enforce: If True, a mismatch of versions will result in an Exception
    :return:        None
    """
    if __version__ == expected:
        print("OK")
    else:
        print(f"*** CAUTION: the installed version of the life123 library ({__version__}) "
              f"doesn't match the stated expected version ({expected}).\nIn case of errors, "
              f"change the installed library version, or modify your code to conform to the installed library."
              f"\nChangelog: https://life123.science/history\n\n")

        if enforce:
            raise Exception(f"Using version '{__version__}' instead of the expected version '{expected}'")


