"""
    life123
    ~~~~~~~

    Dynamical Modeling of Biological Systems

	https://Life123.science

    :copyright:     (c) 2022-2026 by Julian West and the Life123 project.
    :license:       MIT.  See LICENSE file for more details.
"""

__version__ = "1.0.0rc9"    # IN-PROGRESS




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
