import pytest
from src.modules.html_log.html_log import HtmlLog as log
import os



def read_file(filename):
    """
    Read in, and return, the contents of the requested file

    :param filename:
    :return:
    """
    with open(filename, "r") as fh:
        file_contents = fh.read()

    return file_contents



def test_next_available_file():
    # Get rid of test log files that might linger around from earlier tests
    for i in range(1, 5):  # 1 thru 4, inclusive
        filename = f"test{i}.htm"
        if os.path.exists(filename):
            os.remove(filename)


    name = log._next_available_filename("test.htm")
    assert name == "test1.htm"

    with open("test1.htm", "w") as fh:
        fh.write("junk")        # This will be created in the same folder

    name = log._next_available_filename("test.htm")
    assert name == "test2.htm"  # Because "test1.htm" already exists
    with open("test2.htm", "w") as fh:
        fh.write("junk2")

    name = log._next_available_filename("test.htm")
    assert name == "test3.htm"  # Because "test1.htm" and "test2.htm" already exist
    with open("test3.htm", "w") as fh:
        fh.write("junk3")

    # Eliminate one of the files just created
    os.remove("test2.htm")

    name = log._next_available_filename("test.htm")
    assert name == "test2.htm"  # Because "test1.htm" already exists, while "test2.htm" had been deleted
    with open("test2.htm", "w") as fh:
        fh.write("junk2 re-created")

    name = log._next_available_filename("test.htm")
    assert name == "test4.htm"  # Because "test1.htm", "test2.htm" and "test3.htm" already exist
    with open("test4.htm", "w") as fh:
        fh.write("junk4")

    # Cleanup, to avoid leaving test files around
    for i in range(1, 5):  # 1 thru 4, inclusive
        filename = f"test{i}.htm"
        if os.path.exists(filename):
            os.remove(filename)



def test_write_to_file():
    filename = "test_write_to_file.htm"
    if os.path.exists(filename):
        os.remove(filename)

    with open(filename, "w") as fh:
        log.file_handler = fh
        log._write_to_file("this is a test")

    assert read_file(filename) == "this is a test"


    with open(filename, "w") as fh:
        log.file_handler = fh
        log._write_to_file("this is a test", newline=True)

    assert read_file(filename) == "this is a test<br>\n"


    with open(filename, "w") as fh:
        log.file_handler = fh
        log._write_to_file("this is a test", blanks_before=2)

    assert read_file(filename) == "<br>\n<br>\nthis is a test"


    with open(filename, "w") as fh:
        log.file_handler = fh
        log._write_to_file("this is a test", blanks_after=3)

    assert read_file(filename) == "this is a test<br>\n<br>\n<br>\n"

    os.remove(filename)



def test__convert_data():
    d = { "outer_radius": 200 }
    result = log._convert_data(d)
    assert result == "outer_radius_json: 200"

    d = { "outer_radius": 200, "width": True }
    result = log._convert_data(d)
    assert result == "outer_radius_json: 200,\nwidth_json: true"

    d = { "outer_radius": 200, "width": True, "a": [1, 'yes', "no"] }
    result = log._convert_data(d)
    assert result == '''outer_radius_json: 200,\nwidth_json: true,\na_json: [1, "yes", "no"]'''

    d = { "outer_radius": 200, "width": True, "a": [1, 'yes', "no"] , "b": {"x": 8, "y": False} }
    result = log._convert_data(d)
    assert result == '''outer_radius_json: 200,\nwidth_json: true,\na_json: [1, "yes", "no"],\nb_json: {"x": 8, "y": false}'''

    d = { "outer_radius": 200, "width": True, "a": [1, 'yes', "no"] , "b": {"x": 8, "y": False}, "c": "it's \"almost\" true" }
    result = log._convert_data(d)
    assert result == '''outer_radius_json: 200,\nwidth_json: true,\na_json: [1, "yes", "no"],\nb_json: {"x": 8, "y": false},\nc_json: "it's \\"almost\\" true"'''
    print(result)
