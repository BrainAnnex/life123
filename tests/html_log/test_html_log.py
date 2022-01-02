import pytest
from modules.html_log.html_log import HtmlLog as log
import os


def test_next_available_file():
    # Get rid of test log files
    for i in range(3):  # 0 thru 2, inclusive
        filename = f"test{i+1}.htm"
        if os.path.exists(filename):
            os.remove(filename)

    name = log._next_available_filename("test.htm", 3)
    assert name == "test1.htm"

    with open("test1.htm", "w") as fh:
        fh.write("junk")

    name = log._next_available_filename("test.htm", 3)
    assert name == "test2.htm"  # Because "test1.htm" already exists

    with open("test2.htm", "w") as fh:
        fh.write("junk")

    name = log._next_available_filename("test.htm", 3)
    assert name == "test3.htm"  # Because "test1.htm" and "test2.htm" already exist

    with open("test3.htm", "w") as fh:
        fh.write("junk")

    name = log._next_available_filename("test.htm", 3)
    assert name == "test_overflow.htm"  # Because the index would have to be 4 (higher than allowed max)

    # Cleanup, to avoid leaving test files around
    for i in range(3):  # 0 thru 2, inclusive
        filename = f"test{i+1}.htm"
        if os.path.exists(filename):
            os.remove(filename)