################################################################################
#       NOT used in this project!
#       An older library, parts of which were integrated into
#       the class HtmlLog
################################################################################

import logging
import os.path              # To check if a file exists
import traceback
import sys
import html
from datetime import datetime


class old_log:      #####  OBSOLETE!  DON'T USE  #####
    """
    Multi-output channel logging, with optional over-riding directives to include or suppress logging.
    It provides class methods and class properties.
    The default outputs, file names, log directory, etc. can be set in the CONFIGURATION section, below

    OUTPUT MODES:
    1=print ; 2=log file ; 3=UI console

    EXAMPLES:
            log.write("Instantiating the Simulation")       # Writes to the default output modes (they can be set in the CONFIGURATION section, below)

            log.setOutput("debug1", [2])                                        # Anything tagged with the "debug1" context will be output in Mode 2
            log.write("Instantiating the Simulation", context = "debug1")       # Output will be only sent to output mode 2 (i.e., the log file)

            log.setOutput("debug2", None)                                       # (or simply MISSING.)  By default, this will always be assumed whenever contexts are provided
            log.write("Instantiating the Simulation", context = "debug2")       # Nothing will be output this time

            log.setOutput("debug3", "default")                                  # The default output modes will be used for anything tagged with the "debug3" context
            log.write("Instantiating the Simulation", context = "debug3")

    """



    #########################     START OF CONFIG    #########################

    logFilePath = "D:\\Docs\\NeuroComputing\\Valrese\\logs\\" # Default source for the HTML log.  Must include the final slash
    # USE DOUBLE BACK SLASHES FOR WINDOWS PATHS

    logFile = "log"                     # Name of file where to write the HTML log (will be created if not present)
    logFileExtension = ".htm"               # Filename extension for the HTML log file

    separateLogs = True                     # Flag indicating whether each run's output should go into a separate file (consecutively numbered)

    defaultLogs = (1, 2)                    # Log ID or tuple of log ID's to use by default:  1=print ; 2=log file ; 3=UI console

    HTMLindentValue = 10                    # A standardized unit of indentation for the HTML logs (expressed in "ch" units, i.e. the width of a "0")

    MAX_LOG_FILE_NUMBER = 100               # The max number of log files to create (if in separateLogs mode)

    #########################     END OF CONFIG    #########################



    # PRIVATE PROPERTIES
    verbose = {}                    # To affect the printing of any output tagged with a "context" string
    _incrementalTextLog = ""        # Class property storing the Latest round of text log - generally, a condensed logging meant to be shown on the UI (e.g. on a "Log" tab)
    _next_anchor_number = 1         # Used by startLogGroup/endLogGroup, to provide unique anchor points through the HTML log, used to implement "Skip to end" links
    _stack_of_anchor_numbers = []   # Stack implemented with lists.  A stack data structure is needed to deal with nested startLogGroup/endLogGroup pairs



    # Class initialization: set the file-handle class property (fh)         TODO: clarify when this code gets run
    if not separateLogs:    # do an append if the log file already exists, or create it otherwise
        logFileFullName = logFilePath + logFile + logFileExtension
        fh = open(logFileFullName, "a")
    else:                   # create new files, with an auto-increment in the filename.  Example: "log8.htm"
        # Retrieve the first available filename
        index = 1
        logFileFullName: str = logFilePath + logFile + str(index) + logFileExtension

        while os.path.exists(logFileFullName):
            index += 1      # Keep increasing the index in the filename, until an unused filename is found
            logFileFullName: str = logFilePath + logFile + str(index) + logFileExtension
            if index > MAX_LOG_FILE_NUMBER:
                raise Exception("Exceeded the maximum number of allowable log files (" + str(MAX_LOG_FILE_NUMBER) + ")")

        fh = open(logFileFullName, "w")     # Create a new file, to write into

    #print("During import of the 'log' class, the property 'logFileFullName' got set to: " + logFileFullName)


    @classmethod
    def setOutput(cls, context: str, option):
        """

        :param context: A string
        :param option:  Either a list or the string "DEFAULT" or "NONE" (case in-sensitive)
        :return:
        """
        if (option is None) or (option.lower() == "None"):
            cls.verbose[context] = []
        elif option.lower() == "default":
            cls.verbose[context] = cls.defaultLogs
        else:
            # if we get here, option should be either 1, 2 or 3
            # TODO: validate option; if an integer, turn into a list
            cls.verbose[context] = option



    @classmethod
    def error(cls, message, source ="") -> None:
        """
        Log an error message, and then raise an Exception.
        Notice the absence of an HTML vs. plain-text option, compared to the write() method

        :param message: Error message to display
        :param source:  String used to provide information about the source within the class that the call originates from - for example, the name of the method.
                            If present, the message gets prefixed a simple display of the source
        :return:        None
        """
        cls.write("*****  ERROR!  *****\n" + message, source = source, blanksBefore = 2)
        raise Exception(message)



    @classmethod
    def prettyExceptionHandler(cls, exctype, value, tb):
        """
        Pretty-pretty the data from a Python exception, and then TERMINATE the program with exit code 666

        :param exctype:
        :param value:
        :param tb:
        :return:        None (execution gets terminated with exit code 666)
        """

        print("\n==========  HANDLING UNCAUGHT EXCEPTION  ==========")
        #print("Type:", exctype)
        print(F"Type: {exctype.__name__}")
        print(F"Message: `{value}`")
        #print("Traceback:", tb)
        print("TRACEBACK (in order of execution - most recent last):")

        traceback.print_tb(tb, None, sys.stdout)

        print("Exiting....")
        print("===================================================")
        exit(666)



    ###########  PRIMARY METHOD  ############
    @classmethod
    def write(cls, message, plainMsg = None, className = "", source = "", context = None, source_context = None,
              out = defaultLogs, indent = 0, blanksBefore = 0, newline = True, blanksAfter = 0, form = None, showLocation = False) -> None:
        """
        Send the given message to a set of logging locations, specified by the `out` argument,
        unless the current logging configuration (based on className and context) indicates to use different logging locations;
        if argument is missing, the class attribute "defaultLogs" is used instead.

        1=print ; 2=log file ; 3=UI console

        :param message:     A string (possibly with HTML markup), or object that can be sent as input to the str() function,
                            to be sent out to the specified log(s) and/or to printed (newline automatically added by default)

        ALL THE FOLLOWING ARGUMENTS ARE OPTIONAL:
        :param plainMsg:    A plain-text version (possibly blank) of the `message` argument.  If None, then the value in the `message` argument is used.
        :param className:   A string with the name of the invoking class.        Often, self.className is used, provided that the calling class has a defined "className" property
        :param source:      String used to provide information about the source within the class that the call originates from - for example, the name of the method.
                                If present, the message gets prefixed a simple display of the source
        :param context:     String used to implement a granular approach about which logs get withheld and which ones get included
        :param source_context:  String that, if set, overrides the source and context arguments with this same value, plus an extra "()" added to the source
        :param out:         Integer or list of Integers indicating which output(s) to send the message to (1=print ; 2=log file ; 3=UI console)  [Logging mode]
                                If not specified, the class attribute "defaultLogs" is used instead
                                NOTE: this value might get over-ridden by the current logging configuration (based on className and context)
        :param indent:      Integer indicating the number of blank spaces (or HTML left indent) to prefix to the message being output
        :param blanksBefore: Number of blank lines to precede the output with
        :param newline:     Flag indicating whether the output should be terminated by a newline
        :param blanksAfter: Number of blank lines to place after the output
        :param form:        Name of a function (or list of function names) to apply to the message string prior to sending it to HTML logs.
                                Example:  form=log.bold   , or   form=[log.italic, log.red]
        :param showLocation:Flag indicating whether to prominently show the className and the source in a preceding line in the log

        :return:            None
        """

        message = str(message)

        if source_context is not None:
            source = source_context + "()"
            context = source_context

        if context is not None:
            # Check if the logging outputs are to be modified because of directives set for the given context
            revisedMode = cls.directives(context)
            #print("--------------------- log.write revisedMode: ", revisedMode)
            out = revisedMode

        #print("FINAL out mode: ", out)


        # If no plain-message version was provided, make it equal to the message string; since this is meant ONLY for plain-text messages, protect
        # the original message for anything that could be misconstrued as HTML (for example "<some data in triangular bracket>")
        if plainMsg is None:
            plainMsg = message
            message = html.escape(message)


        # For the HTML version, append <br> to all newline characters
        message = message.replace("\n", "\n<br>")


        # Apply any requested HTML formatting to the HTML version of the message
        if form is not None:
            if isinstance(form, list):
                for fn in form:
                    message = fn(message)
            else:
                message = form(message)


        # In the absence of specs calling for prominent displays, do a simple display of the source (if available)
        if source != "" and not showLocation:
            message = "<b>" + source + "</b>: " + message     # Optionally prefix a string describing the "source" of the message (such as the name of the method that is performing the log operation)

            plainMsg = source + ": " + plainMsg


        # Take care of indent, if applicable
        if indent > 0:
            message = log.htmlIndent(indent) + message
            indentString: str = " " * indent        # Repeat a blank character the specified number of times
            plainMsg =  indentString + plainMsg


        if showLocation:
            #  Prominently show the className and the source in a preceding line in the log
            message = "<br><span style='color:#DDDDDD; background-color:#555555; font-style:italic; padding-right:5px; font-size:20px'>" \
                      + "CLASS: " + className + " // METHOD: " + source + " </span><br>" \
                      + message

            plainMsg = "\n[ CLASS: " + className + " // METHOD: " + source + " ]\n" + plainMsg


        # If the `out` argument isn't a tuple, make a tuple of it
        if not isinstance(out, tuple):
            out = (out)


        # DISPATCH TO ALL REQUESTED LOG CHANNELS (1=print ; 2=log file ; 3=UI console)

        if 1 in out:    # Print to console (stdout)
            cls.writeToConsole(plainMsg, blanksBefore, newline, blanksAfter)

        if 2 in out:    # Log to file
            cls.writeToFile(message, blanksBefore, newline, blanksAfter)

        if 3 in out:    # Log to UI console
            cls.writeToUI(message, blanksBefore, newline, blanksAfter)



    @classmethod
    def blankLine(cls, out = defaultLogs, className = "", context = None):
        """
        Append a blank line to the appropriate logs (as specified in the `out` argument, unless over-ridden in logging directives)

        :param out:         Integer or list of Integers indicating which output(s) to send the message to (1=print ; 2=log file ; 3=UI console)  [Logging mode]
                                If not specified, the class attribute "defaultLogs" is used instead
                                NOTE: this value might get over-ridden by the current logging configuration (based on className and context)
        :param className:   A string with the name of the invoking class        Often, self.className is used, provide that the calling class has a defined "className" property
        :param context:     String used to implement a granular approach about which logs get withheld and which ones get included
        :return:            None
        """

        cls.write("", "", out=out, className=className, context=context)



    @classmethod
    def directives(cls, context = ""):
        """
        Consult the class property "verbose" and, based on the given context,
        return a list or tuple of logging modes to use (if a directive was found) or [] , meaning no logging at all (if no applicable directive is present)
        Examples:   [1, 2, 3] ,  [] , [2, 3],  [1], (2, 3)
        An empty list means no logging.

        :param context:     A string used to implement a granular approach about which logs get included
        :return:            A (possibly empty) list or tuple of integers, coding for logging modes, OR None of no directives are present
        """

        if cls.verbose is None or cls.verbose == {}:    # If there are no directives at all
            return []


        # if there are context-based directives
        if context in cls.verbose:
            return cls.verbose[context]

        # If we get thus far, then no applicable directives were found
        return []



    @classmethod
    def toInclude(cls, context = "") -> bool:
        """
        Consult the logging directives for the given context, to determine whether ANY logging of any type
        is to be performed.  Return True if at least 1 type of logging is to be done, or False otherwise.
        This is handy in situations where some data prepping is needed ahead of logging - and it can be skipped if all logging is turned off

        :param context:     A string used to implement a granular approach about which logs get included
        :return:            True if at least 1 type of logging is to be done, or False otherwise
        """
        if cls.directives(context) == []:
            return False
        else:
            return True




    ##########################################   Low-level logging   ##########################################

    @classmethod
    def clearUIbuffer(cls) -> None:
        """
        Clear the UI buffer string
        :return:    None
        """
        cls._incrementalTextLog = ""


    @classmethod
    def UIbuffer(cls) -> str:
        """
        Look up the UI buffer string
        :return:
        """
        return cls._incrementalTextLog


    @classmethod
    def printUIbuffer(cls) -> None:
        """
        Print out the contents of the UI buffer
        :return:    None
        """
        print("_________________________________________________")
        print("BUFFER for the UI:")
        print(cls.UIbuffer())
        print("_________________________________________________")



    @classmethod
    def writeToUI(cls, message: str, blanksBefore = 0, newline = True, blanksAfter = 0) -> None:
        """
        Write the given message into the file managed by the File Handler stored in the class property "fh"

        :param message:     String to send to the UI log
        :param blanksBefore: Number of blank lines to precede the output with
        :param newline:     Flag indicating whether the output should be terminated by a newline
        :param blanksAfter: Number of blank lines to place after the output
        :return:            None
        """

        for i in range(blanksBefore):
            cls._incrementalTextLog += "<br>\n"

        cls._incrementalTextLog += message

        if newline:
            cls._incrementalTextLog += "<br>\n"

        for i in range(blanksAfter):
            cls._incrementalTextLog += "<br>\n"



    @classmethod
    def writeToFile(cls, message: str, blanksBefore = 0, newline = True, blanksAfter = 0) -> None:
        """
        Write the given message into the file managed by the File Handler stored in the class property "fh"

        :param message:     String to write to the designated log file
        :param blanksBefore: Number of blank lines to precede the output with
        :param newline:     Flag indicating whether the output should be terminated by a newline
        :param blanksAfter: Number of blank lines to place after the output
        :return:            None
        """

        for i in range(blanksBefore):
            cls.fh.write("<br>\n")

        cls.fh.write(message)

        if newline:
            cls.fh.write("<br>\n")

        for i in range(blanksAfter):
            cls.fh.write("<br>\n")



    @classmethod
    def writeToConsole(cls, message: str, blanksBefore = 0, newline = True, blanksAfter = 0) -> None:
        """
        Print out (to the standard output) the given message

        :param message:     String to print out to the standard output
        :param blanksBefore: Number of blank lines to precede the output with
        :param newline:     Flag indicating whether the output should be terminated by a newline
        :param blanksAfter: Number of blank lines to place after the output
        :return:            None
        """

        for i in range(blanksBefore):
            print()       # Output a new line

        if newline:
            print(message)              # It includes a newline at the end
        else:
            print(message, end = '')    # Suppress the default newline at the end

        for i in range(blanksAfter):
            print("")   # Output a new line




    ##########################################   High-level logging   ##########################################

    @classmethod
    def newRun(cls) -> None:
        """
        If the HTML log is cumulative, annotate it with a prominent header line (with a timestamp), to indicate the start of a new simulation run
        (to easily distinguish this run from earlier ones.)
        No action taken in case of separate logs.
        Typically done early on in a run.
        Console logs are not affected.

        :return:    None
        """
        if cls.separateLogs:
            return              # Not needed in case of separate logs

        # Prepare timestamp
        now = datetime.now()
        dt_string = now.strftime("%m/%d/%Y %H:%M:%S")   # mm/dd/YY H:M:S

        # Insert into the HTML logs a prominent HTML header, including a timestamp
        msg = "<div style='color:white; background-color:brown; font-size:24px; padding:10px; margin-top:25px'>NEW RUN &nbsp; ("
        msg += dt_string
        msg += ")</div>\n"

        plainMessage = ""   # No message needed for the console log

        cls.write(msg, plainMessage)



    @classmethod
    def separator(cls, caption = "") -> None:
        """
        Append to the appropriate logs a separator line and an optional caption

        :param caption:     Optional caption string after the separator line
        :return:            None
        """
        message = "<hr>\n" + caption
        plainMessage = "------------------------------------------------------------------------------------------------------------\n" + caption
        cls.write(message, plainMessage, form=cls.bold)



    @classmethod
    def logEnter(cls, currentClass: str, currentMethod: str, msg = "") -> None:
        """
        This logging function indicates the start of an important method,
        to be labeled by its class and method, plus an optional message
        Prior to exiting that method, a call should be made to logExit().

        :param currentClass:    The name of the class from which this method was invoked
        :param currentMethod:   The name of the invoking method
        :param msg:             Optional message to add to the log
        :return:                None
        """

        cls.write("<b>ENTERING</b> &nbsp; Class: " + currentClass + " // Method: " + currentMethod,
                  "[ENTERING]    Class: " + currentClass + " // Method: " + currentMethod, blanksBefore=1, className=currentClass)

        # Log the optional message, if present
        if msg != "":
            cls.write(msg, form=cls.italic, className=currentClass)



    @classmethod
    def logExit(cls, currentClass: str, currentMethod: str, msg = "") -> None:
        """
        This call ought to be made just prior to exiting a method that invoked logEnter()

        :param currentClass:    The name of the class from which this method was invoked
        :param currentMethod:   The name of the invoking method
        :param msg:             Optional message to add to the log
        :return:                None
        """

        cls.write("<b>EXITING</b> &nbsp; Class: " + currentClass + " // Method: " + currentMethod,
                  "[EXITING]    Class: " + currentClass + " // Method: " + currentMethod, className=currentClass)

        # Log the optional message, if present
        if msg != "":
            cls.write(msg, form=cls.italic, className=currentClass)

        cls.blankLine(className=currentClass)



    @classmethod
    def startLogGroup(cls, msg = "", include_link = True) -> None:
        """
        Register the start of a related group of logs (such as a listing),
        with an optional message placed at the top as a group header.
        A dashed line (in text mode) or gray DIV box (in HTML mode) is used to contain the log entries inside the group.
        A "Skip to end" link will appear just above the box in HTML mode.
        WARNING: the links will work correctly only if there is a single "run" in the file
        IMPORTANT: At the end, make sure to call endLogGroup()

        :param msg:             Optional message, as a header for the group
        :param include_link:
        """

        # Provide a "Skip to end" link
        # TODO: make it work in cases where there are multiple runs in the same file (maybe incorporate the file's update date in the link... or maybe use a large random number,
        #       shared between startLogGroup and endLogGroup, to try to avoid collisions
        if include_link:
            cls.write(F"<a href='#A_{cls._next_anchor_number}' style='margin-left:" + str(cls.HTMLindentValue) + "ch'>Skip to end</a>", plainMsg = "",
                      newline=False, blanksBefore=1)

        # Note: the anchor numbering is processed normally, regardless of whether the link is actually shown,
        #           so that no special action needs to be taken at endLogGroup()
        cls._stack_of_anchor_numbers.append(cls._next_anchor_number)    # Push number to stack
        cls._next_anchor_number += 1                                    # Autoincrement the number

        # Start a DIV element in HTML mode, and insert a dashed line in text mode
        cls.write("<br><div style='border:1px solid; margin-top:2; margin-left:" + str(cls.HTMLindentValue) + "ch; padding:5px; color:#555555; background-color:#EEEEEE; font-size:11px'>",
                  "------------------------------------------------------------------\n",
                  newline=False)        # Note: newline=False prevents the <br> at the start of the DIV element in the HTML version

        if msg:
            cls.write("<span style='font-size:14px; font-weight:bold'>" + msg + "</span>", msg)     # Optional message, as a header for the group



    @classmethod
    def endLogGroup(cls, msg = "") -> None:
        """
        Visually indicate the end of a related group of logs that were initiated by a call to startLogGroup(),
        with an optional message just before the closure

        :param msg:     Optional message, as a footer just before the closure of the group
        """

        if msg:
            cls.write(msg, form=log.italic)     # Optional message, as a footer just before the closure of the group

        # Provide a unique anchor point for the "Skip to end" link inserted by startLogGroup()
        anchor_number = cls._stack_of_anchor_numbers.pop()              # Pop the stack
        cls.write(F"<a name='A_{anchor_number}'></a>", plainMsg = "")

        # Terminate the DIV element in HTML mode and insert a dashed line in text mode
        cls.write("</div>",
                  "------------------------------------------------------------------")



    ##########################################    HANDY HTML FORMATTING STATIC FUNCTIONS   ##########################################

    @staticmethod
    def htmlIndent(indent: int) -> str:
        """
        Compose and return a SPAN html element to create a left indentation by the specified value (expressed in "ch" units)
        :param indent:  In "ch" units (the width of the "0", zero)
        :return:        A string with the HTML code for the formatted element
        """
        return "<span style='padding-left:" + str(indent) + "ch'>&nbsp;</span>"


    @staticmethod
    def bold(s: str):
        # Example:  form=log.bold   , or   form=[log.italic, log.red]
        return "<b>" + s + "</b>"

    @staticmethod
    def italic(s: str):
        return "<i>" + s + "</i>"

    @staticmethod
    def h1(s: str):
        return "<h1>" + s + "</h1>"

    @staticmethod
    def h2(s: str):
        return "<h2>" + s + "</h2>"

    @staticmethod
    def h3(s: str):
        return "<h3>" + s + "</h3>"

    @staticmethod
    def gray(s: str):
        return "<span style='color:gray'>" + s + "</span>"

    @staticmethod
    def red(s: str):
        return "<span style='color:red'>" + s + "</span>"

    @staticmethod
    def green(s: str):
        return "<span style='color:green'>" + s + "</span>"

    @staticmethod
    def blue(s: str):
        return "<span style='color:blue'>" + s + "</span>"

    @staticmethod
    def reverse(s: str):
        return "<span style='color:white; background-color:black'>" + s + "</span>"

# END class "old_log"



###################################################################################################################

class ExperimentalLog:
    """
    Unused test to replace the earlier class "old_log" (above) but utilizing Python's default logger
    """
    log_file = "old_log.htm"    # Default name
    logger = None


    @classmethod
    def initialize(cls, filename="old_log.htm"):
        cls.log_file = filename

        # Create a custom logger
        cls.logger = logging.getLogger("ExperimentalLog")
        cls.logger.setLevel(logging.INFO)

        # Create handlers
        #c_handler = logging.StreamHandler()
        f_handler = logging.FileHandler(filename)
        #c_handler.setLevel(logging.WARNING)
        #f_handler.setLevel(logging.INFO)

        # Create formatters and add it to handlers
        #c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        #c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)

        # Add handlers to the logger
        #logger.addHandler(c_handler)
        cls.logger.addHandler(f_handler)



    @classmethod
    def write(cls, msg, also_print=False):
        # For the HTML version, append <br> to all newline characters
        msg_html = msg.replace("\n", "\n<br>")

        # Append a final newline
        msg_html += "<br>"

        caller_info = cls.logger.findCaller(stack_info=True)

        print(caller_info)
        cls.logger.info(msg_html)
        cls.logger.warning(msg_html)

        if also_print:
            print(msg)