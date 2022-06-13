"""
A copy of this file should be present in every folder that is not at the project's root level.

Only needed in notebooks that weren't started with a script adding the project root to sys.path
(but no harm if used in those as well.)

*** IMPORTANT: change this relative path as needed in different folder!  (In the line below
that sets PROJECT_HOME_DIR

To use, simply insert the following line at the top of your notebooks:
import set_path
"""

import sys
import pathlib


def add_home_dir_to_syspath():
    PROJECT_HOME_DIR = pathlib.Path.cwd().parent.parent.parent   # IMPORTANT: change this relative path as needed!
    sys.path.append(str(PROJECT_HOME_DIR))
    print(f"Added '{PROJECT_HOME_DIR}' to sys.path")
 

add_home_dir_to_syspath()
