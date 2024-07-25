"""
IMPORTANT: meant *only* for JupyterLab 3 notebooks.

Return the name - WITHOUT extension - of the JupyterLab notebook invoking this function.
This version is based on a modified version of the library "ipynbname",
and appears to work on a variety of systems.

NOTE - 
    MAY NO LONGER BE NECESSARY: ipynbname.name() now appears to work for JupyterLab 3 notebooks.
    To use the ipynbname package from the PyPI release:
        Pip install the following version: ipynbname==2024.1.0.0
        
        import ipynbname
        notebook_name = ipynbname.name()
        
        See: https://github.com/msm1089/ipynbname
"""

from experiments.ipynbname3 import name



def get_notebook_basename() -> str:
    """
    Return the name - WITHOUT extension - of the JupyterLab notebook invoking this function.
    Example:  "reach_equilibrium_1"

    :return:
    """
    return name()
