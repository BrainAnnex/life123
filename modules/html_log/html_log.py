import os.path                  # To check if a file exists
import re                       # Used to strip off HTML
from datetime import datetime



class HtmlLog:
    """
    An HTML logger to file, plus optional plain-text printing to standard output
    """

    #####################
    #  Class variables  #
    #####################

    config_lock = False             # Lock to prevent multiple calls to the config() method
    log_filename = ""               # Including paths and automatically-generated numerical suffixes (eg, "log8.htm"), if applicable
    file_handler = None             # To write into the log file


    # The following class variables (in CAPS) contain the default configuration
    LOG_DIRECTORY = ""              # TODO: test
    LOG_FILENAME_BASE = "log.htm"   # Name for the log file (if using just 1 file), or basename (if using multiple log files)


    SEPARATE_LOGS = False           # Flag indicating whether each run's output should go into a separate file (consecutively numbered);
                                    #   if False, a single file is used

    OVERWRITE_SINGLE_LOG = False    # Only applicable if SEPARATE_LOGS is False; otherwise, it's ignored;
                                    #   in scenarios where a single log file is used, determines if multiple runs over-write or append

    MAX_LOG_FILE_NUMBER = 100       # The max number of log files to create (if SEPARATE_LOGS is True);
                                    #   if the number is exceed, the logs are appended to a file with an "_overflow" suffix in the name

    ALSO_PRINT = True               # If True, anything sent to the log will out be sent to stdout
    EXTRA_HANDLERS = []             # NOT IN CURRENT USE.  List of functions to invoke with after sending anything to the log



    #########################################
    #            CONFIGURATIONS             #
    #########################################

    @classmethod
    def config(cls, filename=None, multiple=None, overwrite=None, max_files=None) -> None:
        """
        Roughly, the counterpart of basicConfig in the default Python logger

        :param filename:
        :param multiple:
        :param overwrite:
        :param max_files:
        :return:                None
        """
        assert not cls.config_lock, "the config() method can only be invoked once per run"

        # Change whatever configuration variables are being passed, if any
        if filename:
            cls.LOG_FILENAME_BASE = filename

        if multiple:
            cls.SEPARATE_LOGS = multiple

        if overwrite:
            cls.OVERWRITE_SINGLE_LOG = overwrite

        if max_files:
            cls.MAX_LOG_FILE_NUMBER = max_files


        # Generate, if applicable, a new name for the log file (applicable in the case of multiple logs),
        #       and create a file handle to write into it
        if cls.SEPARATE_LOGS:   # Create new files, with an auto-increment in the filename.  Example: "log8.htm"
            # Retrieve the first available filename
            cls.log_filename = cls._next_available_filename(cls.LOG_DIRECTORY + cls.LOG_FILENAME_BASE,
                                                            cls.MAX_LOG_FILE_NUMBER)
            cls.file_handler = open(cls.log_filename, "a")  # Create a new file, to write to
                                                            # (the append part is only relevant if we reach the max # of files)

        else:                   # Using a single log file
            cls.log_filename = cls.LOG_DIRECTORY + cls.LOG_FILENAME_BASE
            if cls.OVERWRITE_SINGLE_LOG:
                cls.file_handler = open(cls.log_filename, "w")     # Create a new file, to write into: over-written if present
            else:
                cls.file_handler = open(cls.log_filename, "a")     # Create a new file, to append to


        print(f"Output will be logged into the file '{cls.log_filename}'")

        cls.config_lock = True       # This intercepts and prevents additional calls to this method in the same run



    #########################################
    #            PRIMARY METHOD             #
    #########################################
    @classmethod
    def write(cls, msg: str, plain="", also_print=False,
              blanks_before = 0, newline = True, blanks_after = 0,
              style = None):
        """

        :param msg:             A string (possibly with HTML markup - but NOT recommended),
                                or object that can be sent as input to the str() function,
                                to be sent out to the specified log(s) and/or to printed (newline automatically added by default)

        ALL THE FOLLOWING ARGUMENTS ARE OPTIONAL
        :param plain:           A plain-text version (possibly blank) of the `msg` argument.
                                If an empty string, then the value in the `msg` argument is used,
                                after stripping off any HTML that might be present
        :param also_print:      Only applicable if `plain` is blank.  Flag indicating whether to also print an HTML-stripped version of msg

        :param blanks_before:   Number of blank lines to precede the output with
        :param newline:         Flag indicating whether the output should be terminated by a newline
        :param blanks_after:    Number of blank lines to place after the output

        :param style:           Name of a function (or list/tuple of function names) to apply to the message string prior to sending it to HTML logs.
                                Example:  HtmlLog.bold   , or   [HtmlLog.italic, HtmlLog.red]
        :return:
        """
        # Check if the configuration was run
        if not cls.file_handler:
            cls.config()            # Set up the default configuration

        msg = str(msg)              # To convert numbers, or objects that can be sent as input to the str() function


        # Determine if something should also be sent to the standard output
        if plain or also_print or cls.ALSO_PRINT:
            if not plain:
                # If no plain-message version was provided, make it equal to the message string,
                # stripped on any HTML that might be present
                plain = cls.strip_html(msg)

            # Dispatch the request to the handler for the standard output
            cls._write_to_console(plain, blanks_before=blanks_before, newline=newline, blanks_after=blanks_after)

        # The remainder of the function is for processing the HTML logging

        # For the HTML version, append <br> to all newline characters
        msg = msg.replace("\n", "\n<br>")

        # Apply any requested HTML formatting to the HTML version of the message
        if style is not None:
            #if isinstance(style, list):
            if type(style) == list or type(style) == tuple:
                # If multiple style functions were passed (e.g., one for boldface and one for color), apply each in turn
                for fn in style:
                    msg = fn(msg)
            else:
                msg = style(msg)    # A single style function is applied

        # Dispatch the request to the handler for the HTML file output
        cls._write_to_file(msg, blanks_before=blanks_before, newline=newline, blanks_after=blanks_after)



    ##########################################
    #           HIGH-LEVEL LOGGING           #
    ##########################################

    @classmethod
    def new_run(cls) -> None:
        """
        If the HTML log is cumulative, annotate it with a prominent header line (with a timestamp),
        to indicate the start of a new run (to easily distinguish this run from earlier ones.)
        Typically done early on in a run.

        :return:    None
        """
        #if cls.SEPARATE_LOGS:
            #return              # Not needed in case of separate logs

        # Prepare timestamp
        now = datetime.now()
        dt_string = now.strftime("%m/%d/%Y %H:%M:%S")   # mm/dd/YY H:M:S

        # Insert into the HTML logs a prominent HTML header, including a timestamp.  TODO: switch from style to class ?
        msg = "<div style='color:white; background-color:brown; font-size:24px; padding:10px; margin-top:25px'>NEW RUN &nbsp; ("
        msg += dt_string
        msg += ")</div>\n"

        plain_message = f"NEW RUN ({dt_string})"

        cls.write(msg, plain=plain_message)



    @classmethod
    def separator(cls, caption = "", also_print=False) -> None:
        """
        Append to the appropriate logs a separator line and an optional caption

        :param caption:     Optional caption string after the separator line
        :param also_print:  Flag indicating whether to also print a text version
        :return:            None
        """
        msg = "<hr>\n" + cls.bold(caption)
        if also_print or cls.ALSO_PRINT:
            plain_message = "------------------------------------------------------------------------------------------------------------\n" \
                            + caption
            cls.write(msg, plain=plain_message)
        else:
            cls.write(msg)



    ##########################################   Styling : HTML FORMATTING FUNCTIONS  ##########################################

    @staticmethod
    def bold(s: str) -> str:
        # EXAMPLE:  to use this as an argument to write(),
        #           pass a parameter such as style=HtmlLog.bold   , or   style=[HtmlLog.bold, HtmlLog.red]
        return "<b>" + s + "</b>"

    @staticmethod
    def italic(s: str) -> str:
        return "<i>" + s + "</i>"

    @staticmethod
    def h1(s: str) -> str:
        return "<h1>" + s + "</h1>"

    @staticmethod
    def h2(s: str) -> str:
        return "<h2>" + s + "</h2>"

    @staticmethod
    def h3(s: str) -> str:
        return "<h3>" + s + "</h3>"

    @staticmethod
    def gray(s: str) -> str:
        return "<span style='color:gray'>" + s + "</span>"

    @staticmethod
    def red(s: str) -> str:
        return "<span style='color:red'>" + s + "</span>"

    @staticmethod
    def green(s: str) -> str:
        return "<span style='color:green'>" + s + "</span>"

    @staticmethod
    def blue(s: str) -> str:
        return "<span style='color:blue'>" + s + "</span>"

    @staticmethod
    def color(s: str, color_value) -> str:
        # Note that this function takes a 2nd argument
        return f"<span style='color:{color_value}'>" + s + "</span>"

    @staticmethod
    def reverse(s: str) -> str:
        return "<span style='color:white; background-color:black'>" + s + "</span>"



    @staticmethod
    def html_indent(indent: int) -> str:
        """
        Compose and return a SPAN html element to create a left indentation by the specified value (expressed in "ch" units)
        :param indent:  In "ch" units (the width of the zero character "0")
        :return:        A string with the HTML code for the formatted element
        """
        return f"<span style='padding-left:{indent}ch'>&nbsp;</span>"



    @staticmethod
    def link(name: str, url: str, new_window=True) -> str:
        """

        :param name:
        :param url:
        :param new_window:
        :return:
        """
        if new_window:
            return f"<a href ='{url}' target='_blank'>{name}</a>"
        else:
            return f"<a href ='{url}'>{name}</a>"


    @staticmethod
    def button_post(text: str, url: str, post_data: str) -> str:
        """

        :param text:
        :param url:
        :param post_data: TODO: it currently cannot have double quotes
        :return:
        """
        return f'''
        <form method="POST" action='{url}'>
        <input type="submit" value="{text}">
        <input type="hidden" name="post_data" value="{post_data}">
        </form>
        '''




    ##########################################   Low-level logging (handlers for the logging)   ##########################################

    @classmethod
    def _write_to_file(cls, msg: str, blanks_before = 0, newline = True, blanks_after = 0) -> None:
        """
        Write the given message into the file managed by the File Handler stored in the class property "file_handler"

        :param msg:             String to write to the designated log file
        :param blanks_before:   Number of blank lines to precede the output with
        :param newline:         Flag indicating whether the output should be terminated by a newline
        :param blanks_after:    Number of blank lines to place after the output
        :return:                None
        """

        for i in range(blanks_before):
            cls.file_handler.write("<br>\n")

        cls.file_handler.write(msg)

        if newline:
            blanks_after += 1

        for i in range(blanks_after):
            cls.file_handler.write("<br>\n")



    @classmethod
    def _write_to_console(cls, msg: str, blanks_before = 0, newline = True, blanks_after = 0) -> None:
        """
        Print out (to the standard output) the given message

        :param msg:             Plain-text string to print out to the standard output
        :param blanks_before:   Number of blank lines to precede the output with
        :param newline:         Flag indicating whether the output should be terminated by a newline
        :param blanks_after:    Number of blank lines to place after the output
        :return:                None
        """
        for i in range(blanks_before):
            print()         # Output a new line

        if newline:
            print(msg)              # It includes a newline at the end
        else:
            print(msg, end = '')    # Suppress the default newline at the end

        for i in range(blanks_after):
            print()         # Output a new line



##########################################   Utilities   ##########################################

    @classmethod
    def strip_html(cls, html_str) -> str:
        """
        A light-duty HTML stripper that preserves newlines.
        EXAMPLE: "An <b>important</b> list:<br>A<br>B" results in "An important list:\nA\nB"
        :param html_str:
        :return:
        """
        html_str = html_str.replace('<br>', '\n')
        html_str = html_str.replace('<br\>', '\n')
        plain_str = re.sub('<[^<]+?>', '', html_str)

        return plain_str



    @classmethod
    def _next_available_filename(cls, filename_stem: str, max_index: int) -> str:
        """
        Given a particular file name (referred to as "stem", because we'll be added a numeric index),
        locate the first index such that no file with that name is present in the same directory as the original file.
        If the index would exceed the specified max, use the string "_overflow" in lieu of a numeric index.
        EXAMPLES:
            given "test.htm", it might return "test1.htm" or "test12.htm" or "test_overflow.htm"
            given "/my_path/test.htm", it might return "/my_path/test1.htm" or "/my_path/test12.htm" or "/my_path/test_overflow.htm"
            given "myfile", it might return "myfile1" or "myfile3" or "myfile_overflow"

        :param filename_stem:   A filename with or without a path.  If no path is given, the current directory is understood
        :param max_index:       The maximum number of individual file we want
        :return:                The modified filename, with either a numeric index or "_overflow" added at the end
                                (file suffices such as ".htm" aren't modified if present)
        """
        # Split the given filename into a base and a suffix
        (basename, suffix) = os.path.splitext(filename_stem)    # EXAMPLES: "test.jpg" becomes ("test", ".jpg")
        #           "myfile" becomes ("myfile", "")
        #           "/my_path/example.txt" becomes ("/my_path/example", ".txt")

        index = 1    # Tentatively try the first index number

        filename = basename + str(index) + suffix

        while os.path.exists(filename):
            index += 1      # Keep increasing the index in the filename, until an unused filename is found
            if index > max_index:
                break

            filename = basename + str(index) + suffix

        if index > max_index:
            filename = basename + "_overflow" + suffix
            print(f"Exceeded the maximum number of allowable log files ({max_index}).  From now on, using the file '{filename}'")
        else:
            filename = basename + str(index) + suffix

        return filename
