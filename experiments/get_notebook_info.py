"""
IMPORTANT: meant *only* for JupyterLab 3 notebooks.

Return the name - WITHOUT extension - of the JupyterLab notebook invoking this function.
This version is based on a modified version of the library "ipynbname",
and appears to work on a variety of systems.
"""

from experiments.ipynbname3 import name



def get_notebook_basename() -> str:
    """
    Return the name - WITHOUT extension - of the JupyterLab notebook invoking this function.
    Example:  "reach_equilibrium_1"

    :return:
    """
    return name()
