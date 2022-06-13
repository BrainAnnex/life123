"""
A copy of this file should be present in every folder that is not at the project's root level.

Only needed when JupyterLab isn't started with an option that adds the project root to sys.path
(as done, for example, in the "quick.bat" script at the top level.)
There's no harm in using this module under all circumstances.


TO USE:
Insert the following two lines at the top of your notebooks:

import set_path
set_path.add_ancestor_dir_to_syspath(3)     # IMPORTANT: change this number as needed
"""

import sys
import pathlib


def add_ancestor_dir_to_syspath(n: int, verbose=True) -> None:
    """
    Add the n-th ancestor of the current working directory to the system path

    :param n:
    :param verbose:
    :return:        None
    """
    directory = pathlib.Path.cwd()      # Start with the current working directory

    for _ in range(n):
        directory = directory.parent    # Navigate up one level in the file hierarchy

    #PROJECT_HOME_DIR = pathlib.Path.cwd().parent.parent.parent   # IMPORTANT: change this relative path as needed!
    sys.path.append(str(directory))

    if verbose:
        print(f"Added '{directory}' to sys.path")


#add_ancestor_dir_to_syspath(3)
