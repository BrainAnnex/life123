"""
A copy of this file should be present in every notebook folder
that isn't the Life123 project's root folder.

Only needed when JupyterLab isn't started with an option
that adds the project's root folder to sys.path
(as done, for example, in the "quick.bat" script in Life123 root folder.)

However, there's no harm in using this module under all circumstances,
even if the project's root folder is already present in sys.path.

Adding the project's root folder to sys.path permits
the use of absolute paths in import of modules into notebooks.

TO USE:
------
Insert the following two lines at the top of your notebooks:

import set_path
set_path.add_ancestor_dir_to_syspath(3)     # IMPORTANT: change this number as needed
"""

import sys
import pathlib


def add_ancestor_dir_to_syspath(n: int, verbose=True) -> None:
    """
    Add the n-th ancestor of the current working directory to the system path.
    Doing so will allow the use of absolute paths (starting at that ancestor folder)
    in import of modules into notebooks.

    :param n:       The number of levels to go up,
                        to reach the project's home from the folder containing this notebook
    :param verbose: Flag indicating whether to print out the directory newly added to syspath
    :return:        None
    """
    directory = pathlib.Path.cwd()      # Start with the current working directory

    for _ in range(n):
        directory = directory.parent    # Navigate up one level in the file hierarchy

    sys.path.append(str(directory))

    if verbose:
        print(f"Added '{directory}' to sys.path")
