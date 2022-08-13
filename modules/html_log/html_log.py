import os.path                  # To check if a file exists
import re                       # Used to strip off HTML
from datetime import datetime
import json



class HtmlLog:
    """
    An HTML logger to file, plus optional plain-text printing to standard output
    TODO:   maybe rename "HtmlLogger" or "ReportMaker"
            Add feature to create a Table of Contents
    """

    #####################
    #  Class variables  #
    #####################

    config_lock = False             # Lock to prevent multiple calls to the config() method
    log_fullname = ""               # Including optional path and, if applicable, automatically-generated numerical suffixes (eg, "D:/Docs/log8.htm")
    file_handler = None             # To write into the log file


    # The following class variables (in CAPS) contain the default configuration
    LOG_DIRECTORY = ""              # TODO: test
    LOG_FILENAME_BASE = "default_log.htm"   # Name for the log file (if using just 1 file), or basename (if using multiple log files)


    SEPARATE_LOGS = False           # Flag indicating whether each run's output should go into a separate file (consecutively numbered);
                                    #   if False, a single file is used

    OVERWRITE_SINGLE_LOG = False    # Only applicable if SEPARATE_LOGS is False; otherwise, it's ignored;
                                    #   in scenarios where a single log file is used, determines if multiple runs over-write or append

    ALSO_PRINT = True               # Default about whether to also send plain-text version of the log messages to stdout
    EXTRA_HANDLERS = []             # NOT IN CURRENT USE.  List of functions to invoke after sending anything to the log

    use_D3 = False
    use_Vue= False
    VUE_COUNT = 1               # Auto-increment unique number for possible multiple Vue apps on the same page;
                                #   each <DIV> containing the Vue component will have an ID of the form "vue-root-VUE_ID"



    #########################################
    #            CONFIGURATIONS             #
    #########################################

    @classmethod
    def config(cls, filename=None, namestem=None, mode=None,
                    css=None, js=None, use_D3=False, Vue_lib=False) -> None:
        """
        It can only be called once in a run.
        If desired, it can be done explicitly - and should be done so if the Class defaults need to be changed;
        if not done explicitly, it'll be done automatically if calling any method in this library.

        Roughly, the counterpart of basicConfig in the default Python logger

        # ARGUMENTS OPTIONALLY USED TO CHANGE THE CLASS DEFAULTS
        :param filename:    Name for the log file (if using just 1 file), or basename (if using multiple log files);
                                EXAMPLES: "log.htm" or "reports/log1.htm"

        :param namestem:    EXPERIMENTAL, not yet in use.  A filename, or full/relative path, without a suffix;
                                ".log.htm" will automatically get added

        :param mode:        Must be one of:
                                "append"    - if the log file exists, append to it; otherwise, create it (AVOID if using Vue)
                                "overwrite" - if the log file exists, first clear it; otherwise, create it
                                "multiple"  - each run's output should go into a separate file (consecutively numbered)

        # ARGUMENTS OPTIONALLY USED TO PASS STYLING/JAVASCRIPT/GRAPHING PARAMETER
        :param css:         String, or list of strings, with name(s) of CSS files to include
        :param js:          String, or list of strings, with name(s) of extra JavaScript files to include
        :param use_D3:      Flag indicating whether D3 will be used.  If True,
                                https://d3js.org/d3.v7.min.js will be included
        :param Vue_lib:     Full name of js file that contains the desired version of the Vue library.
                                NOTE: at present, use of Vue isn't compatible with the "append" mode

        :return:            None
        """
        assert not cls.config_lock, "The config() method can only be invoked ONCE per run.  Did you already make a call to config() or write() ?"

        # Change the default values of any configuration variables that are explicitly being passed
        if filename:
            cls.LOG_FILENAME_BASE = filename

        if mode == "multiple":
            cls.SEPARATE_LOGS = True
        elif mode == "overwrite":
            cls.SEPARATE_LOGS = False
            cls.OVERWRITE_SINGLE_LOG = True
        elif mode == "append":
            cls.SEPARATE_LOGS = False
            cls.OVERWRITE_SINGLE_LOG = False
        elif mode is not None:
            raise Exception(f"The `mode` argument, if used, must be one of : 'multiple' , 'overwrite' , 'append'. The value passed was `{mode}`")


        # Generate, if applicable, a new name for the log file (applicable in the case of multiple logs),
        #       and create a file handle to write into it
        if cls.SEPARATE_LOGS:   # Create new files, with an auto-increment in the filename.  Example: "log8.htm"
            # Retrieve the first available filename
            cls.log_fullname = cls._next_available_filename(cls.LOG_DIRECTORY + cls.LOG_FILENAME_BASE)
            cls.file_handler = open(cls.log_fullname, "w")          # Create a new file, to write to: over-write if present

        else:                   # Using a single log file
            cls.log_fullname = cls.LOG_DIRECTORY + cls.LOG_FILENAME_BASE
            if cls.OVERWRITE_SINGLE_LOG:
                cls.file_handler = open(cls.log_fullname, "w")      # Create a new file, to write to; over-write if present
                cls._html_comment("Auto-generated log file", newline_before = False)
            else:
                cls.file_handler = open(cls.log_fullname, "a")      # Create a new file, to append to

                
        # Give feedback about where the output is being logged
        print(f"-> Output will be LOGGED into the file '{cls.log_fullname}'")


        # PREPARE THE HEADER

        # Handle the required CSS file(s)
        if not css:
            css_lines = ""
        elif type(css) == str:
            css_lines = f'    <link type="text/css" rel="stylesheet" href="{css}">'
        elif type(css) == list:
            css_list = [f'    <link type="text/css" rel="stylesheet" href="{css_file}">'
                            for css_file in css]
            css_lines = "\n".join(css_list)
        else:
            raise Exception(f"Argument css, if passed, must be a string or a list of strings; instead, it was a {type(css)}")


        # Handle the required JavaScript file(s)
        js_lines = ''

        if Vue_lib:
            if not cls.SEPARATE_LOGS and not cls.OVERWRITE_SINGLE_LOG:
                raise Exception("At present, Vue cannot be used in the 'append' mode: use mode='overwrite' or mode='multiple'")

            js_lines += f'\n    <script src="{Vue_lib}"></script>'
            use_Vue = True
        else:
            use_Vue = False


        if use_D3:
            js_lines += '\n    <script src="https://d3js.org/d3.v7.min.js"></script>'


        if js:
            if type(js) == str:
                js_lines += f'\n    <script src="{js}"></script>'
            elif type(js) == list:
                for js_file in js:
                    js_lines += f'\n    <script src="{js_file}"></script>'
        else:
            raise Exception(f"Argument js, if passed, must be a string or a list of strings; instead, it was a {type(js)}")


        cls.use_D3 = use_D3
        cls.use_Vue = use_Vue


        # Assemble the initial part of the HTML file, thru the opening of the <BODY> tag
        initial_html = f'''<!DOCTYPE html>
<html lang="en">
<head>  
    <meta charset="UTF-8">
    <title>Log</title>

{css_lines}
{js_lines}
</head>

<body>\n'''

        cls._write_to_file(initial_html)

        cls.config_lock = True       # This intercepts and prevents additional calls to this method in the same run




    ###############################################
    #               PRIMARY METHOD                #
    ###############################################

    @classmethod
    def write(cls, msg: str, plain="", also_print=ALSO_PRINT,
              indent = 0,
              blanks_before = 0, newline = True, blanks_after = 0,
              style = None, style_par = None):
        """

        :param msg:             A string (possibly with HTML markup - but NOT recommended),
                                or object that can be sent as input to the str() function,
                                to be sent out to the log, and possibly also printed
                                (newline automatically added by default)

        ALL THE FOLLOWING ARGUMENTS ARE OPTIONAL
        :param plain:           A plain-text version (possibly blank) of the `msg` argument.
                                If an empty string, then the value in the `msg` argument is used,
                                after stripping off any HTML that might be present
        :param also_print:      Only applicable if `plain` is blank.
                                Flag indicating whether to also print an HTML-stripped version of msg

        :param indent:          Integer indicating the number of blank spaces (or HTML left indent) to prefix to the message being output

        :param blanks_before:   Number of blank lines to precede the output with
        :param newline:         Flag indicating whether the output should be terminated by a newline (a <br> in the HTML version)
        :param blanks_after:    Number of blank lines to place after the output

        :param style:           Name of a function (or list/tuple of function names) to apply to the message string prior to sending it to HTML logs.
                                Example:  HtmlLog.bold   , or   [HtmlLog.italic, HtmlLog.red]
        :param style_par:

        :return:
        """

        # Check if the configuration was run
        if not cls.file_handler:
            print("The logger module is being AUTOMATICALLY set to its default configuration (to set it manually, make a call to config() first)")
            cls.config()            # Set up the default configuration

        msg = str(msg)              # To convert numbers, or objects that can be sent as input to the str() function


        # Determine if something should also be sent to the standard output
        if plain or also_print:
            if not plain:
                # If no plain-message version was provided, make it equal to the message string,
                # stripped on any HTML that might be present
                plain = cls.strip_html(msg)

            # Take care of indent, if applicable
            if indent > 0:
                indent_str = " " * indent        # Repeat a blank character the specified number of times
                plain =  indent_str + plain

            # Dispatch the request to the handler for the standard output
            cls._write_to_console(plain, blanks_before=blanks_before, newline=newline, blanks_after=blanks_after)


        # The remainder of the function is for processing the HTML logging

        # For the HTML version, append <br> to all newline characters
        msg = msg.replace("\n", "\n<br>")

        # Apply any requested HTML formatting to the HTML version of the message
        if style is not None:
            if type(style) == list or type(style) == tuple:
                # If multiple style functions were passed (e.g., one for boldface and one for color), apply each in turn
                for fn in style:
                    msg = fn(msg)

            else:   # A single style function is applied
                if style_par:
                    msg = style(msg, style_par)
                else:
                    msg = style(msg)


        # Take care of indent, if applicable
        if indent > 0:
            msg = cls.html_indent(indent) + msg

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
            #return              # Not needed in case of separate logs      TODO: check this

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
    def separator(cls, caption = "", also_print=ALSO_PRINT) -> None:
        """
        Append to the appropriate logs a separator line and an optional caption

        :param caption:     Optional caption string after the separator line
        :param also_print:  Flag indicating whether to also print a text version
        :return:            None
        """
        msg = "<hr>\n" + cls.bold(caption)
        if also_print:
            plain_message = "------------------------------------------------------------------------------------------------------------\n" \
                            + caption
            cls.write(msg, plain=plain_message)
        else:
            cls.write(msg)



    @classmethod
    def blank_line(cls, n_blanks = 1) -> None:
        """
        Append the specified number of blank lines to the logs

        :param n_blanks:    Desired number of blank lines
        :return:            None
        """

        cls.write("", "", blanks_after = n_blanks)



    @classmethod
    def export_plot_D3(cls, data, svg_id, js_file, js_func, D3_tag="svg") -> None:
        """
        Append to the log the HTML needed to produce a D3-based plot

        :param data:
        :param svg_id:
        :param js_file:
        :param js_func:
        :param D3_tag:  The tag of the container that D3 reaches into.
                        Typically, "svg" or "div"
        :return:        None
        """
        if not cls.use_D3:
            cls.write("ERROR: In order to utilize D3 plots, the  use_D3=True  option must be used in the call to config()", style=cls.red)
            return

        cls._html_comment(f"Start of D3 plot (id `{svg_id}`)")

        cls._external_js_file(js_file, defer=False)

        html = f"<{D3_tag} id='{svg_id}'></{D3_tag}>\n"
        cls._write_to_file(html)

        # Start of JS script
        cls._write_to_file("<script defer>\n")

        js_data = "var json_data = '"       # Note the opening single quote...
        cls._write_to_file(js_data)

        json.dump(data, cls.file_handler)       # Convert the date to JS, and write it to the log file
        #json_str = json.dumps(data, indent=4)  # Without an indent value, the whole string is in one line

        js_data = f"';\nvar DATA_1 = JSON.parse(json_data);\n{js_func}(\"{svg_id}\", DATA_1);\n</script>"
        cls._write_to_file(js_data)
        # End of JS script

        cls._html_comment("End of D3 plot")



    @classmethod
    def export_plot_Vue(cls, data: dict, component_name: str, component_file: str) -> None:
        """
        Append to the log the HTML needed to produce a Vue-based plot

        :param data:            A python dictionary of data to pass to the Vue component
        :param component_name:  A string with the name of the existing Vue.js component to use.
                                        EXAMPLE: "vue_curves_4" (assuming that a js file with such a component exists)
        :param component_file:  A string with the name of the .js file containing the needed Vue component above

        :return:        None
        """
        if not cls.use_Vue:
            cls.write("ERROR: In order to utilize Vue, the  use_Vue=True  option must be used in the call to config()", style=cls.red)
            return

        # Validate args
        assert type(data) == dict,\
            f"export_plot_Vue(): the `data` argument must be a dictionary; the value passed was {type(data)}"
        #TODO: more validation

        data["component_id"] = cls.VUE_COUNT

        vue_id = f"vue-root-{cls.VUE_COUNT}"
        cls._write_to_file(f'\n\n<div id="{vue_id}">   <!-- DIV container for the VUE COMPONENTS below : Vue ROOT element -->\n')

        component_call = f'''
    <{component_name} {cls._pass_props(data)}
    >        
    </{component_name}>
    '''

        cls._write_to_file(component_call)

        cls._write_to_file('\n</div>	<!--  ~~~~~~~~~~~~~~~~~~~~~  END of Vue root element  ~~~~~~~~~~~~~~~~~  -->\n')

        html = f'''

<!--
    Vue components (and other JS).  This must appear AFTER the Vue-containing elements
  -->
<script src="{component_file}"></script>
'''
        cls._write_to_file(html)

        instantiate_vue = f'''
<script>
// Instantiation of the Vue ROOT component
new Vue({{
    el: '#{vue_id}',

    data: {{
        {cls._convert_data(data)}
    }}

}});
</script>

'''
        cls.VUE_COUNT += 1     # Auto-increment the ID to use for possible multiple Vue apps on the same page

        cls._write_to_file(instantiate_vue, newline=False)



    @staticmethod
    def _pass_props(data: dict) -> str:
        """
        Prepare a string that goes inside the Vue component call

        :param data:    A python dictionary of data for the component
        :return:        A string that goes inside the Vue component call
        """
        if data is None:
            return ""

        s = ""
        for key in data:
            s +=  f'\n        v-bind:{key}="{key}_json"'
            #s +=  f'\n        v-bind:{key}="JSON.parse({key}_json)"'   # In some tests, this extra conversion could be avoided

        return s


    @staticmethod
    def _convert_data(data: dict) -> str:
        """
        From a Python dictionary, prepare a string that goes inside the "data:" section of the Vue element instantiation,
        and therefore needs to be valid JavaScript.

        EXAMPLE:  { "outer_radius": 200, "width": True, "a": [1, 'yes', "no"] , "b": {"x": 8, "y": False}, "c": "it's \"almost\" true" }
            returns the string:
            '''outer_radius_json: 200,\nwidth_json: true,\na_json: [1, "yes", "no"],\nb_json: {"x": 8, "y": false},\nc_json: "it's \\"almost\\" true"'''
            which prints out as [left indent added for readability]:

                outer_radius_json: 200,
                width_json: true,
                a_json: [1, "yes", "no"],
                b_json: {"x": 8, "y": false},
                c_json: "it's \"almost\" true"

        Notice the "_json" suffix added to the keys  (TODO: maybe this suffix could be ditched...)

        :param data:    A python dictionary
        :return:        A string that goes inside the "data:" section of the Vue element instantiation
        """
        if data is None:
            return ""

        """
        # In numerous tests, this extra conversion could be avoided
        l = [f"{key}_json: '{json.dumps(value)}'"
                        for (key, value) in data.items()]   # Notice the "_json" suffix added to keys
                                                            # and the values converted to JSON and placed in single quotes
        """
        l = []
        for (key, value) in data.items():
            l.append(f"{key}_json: {json.dumps(value)}")    # Notice the "_json" suffix added to keys

        #print(l)
        res = ",\n".join(l)
        return res





    ##########################################   Styling : HTML FORMATTING FUNCTIONS  ##########################################

    @staticmethod
    def bold(s: str) -> str:
        # EXAMPLE:  to use this as an argument to write(),
        #           pass a parameter such as  style=HtmlLog.bold   , or   style=[HtmlLog.bold, HtmlLog.red]
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
        # Note that this function takes a 2nd argument, which requires use of style_par argument.
        # EXAMPLES of usage:    style=HtmlLog.color, style_par='gray'
        #                       style=HtmlLog.color, style_par='#DDD'
        return f"<span style='color:{color_value}'>" + s + "</span>"

    @staticmethod
    def reverse(s: str) -> str:
        return "<span style='color:white; background-color:black'>" + s + "</span>"



    @staticmethod
    def html_indent(indent: int) -> str:
        """
        Compose and return a SPAN html element to create a left indentation
        by the specified value (expressed in "ch" units)

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



    ##########################################
    #           JavaScript support           #
    ##########################################

    @classmethod
    def _external_js_file(cls, filename: str, defer=False):
        """
        Write into the log file the necessary HTML code to include
        the given JavaScript file

        :param filename:
        :param defer:
        :return:
        """
        if defer:
            code = f"<script src='{filename}' defer></script>"
        else:
            code = f"<script src='{filename}'></script>"

        cls._write_to_file(code + "\n", newline=False)




    ###########################   Low-level logging (handlers for the logging)   ###########################

    @classmethod
    def _write_to_file(cls, msg: str, blanks_before = 0, newline = False, blanks_after = 0) -> None:
        """
        Write the given message (containing text and/or HTML code) into the file managed with
        the File Handler stored in the class property "file_handler"

        :param msg:             String to write to the designated log file
        :param blanks_before:   Number of blank lines ("<br>\n") to precede the output with
        :param newline:         Flag indicating whether the output should be terminated by an HTML newline ("<br>\n")
        :param blanks_after:    Number of blank lines ("<br>\n") to place after the output
        :return:                None
        """

        for i in range(blanks_before):
            cls.file_handler.write("<br>\n")

        cls.file_handler.write(msg)

        if newline:
            blanks_after += 1

        for i in range(blanks_after):
            cls.file_handler.write("<br>\n")

        cls.file_handler.flush()    # To avoid weird buffering issues seen in JupyterLab
        


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
    def close(cls):
        """
        Useful in case of ongoing file lock the prevents its deletion
        (for example, if the file was opened in JupyterLab)
        """
        cls.file_handler.close()
        
        
    
    @classmethod
    def _html_comment(cls, comment: str, newline_before = True) -> None:
        """
        Write to the log, an HTML comment (useful for readability), on a separate line

        :param comment:         String with the body of the comment.  EXAMPLE: "Start of data table"
        :param newline_before:
        :return:                None
        """
        s = f"<!-- {comment} -->\n"
        if newline_before:
            s = "\n" + s

        cls._write_to_file(s)



    @classmethod
    def strip_html(cls, html_str) -> str:
        """
        A light-duty HTML stripper that preserves newlines.
        EXAMPLE:    "An <b>important</b> list:<br>A<br>B"
                    results in
                    "An important list:\nA\nB"

        :param html_str:
        :return:
        """
        html_str = html_str.replace('<br>', '\n')
        html_str = html_str.replace('<br\>', '\n')
        plain_str = re.sub('<[^<]+?>', '', html_str)

        return plain_str



    @classmethod
    def _next_available_filename(cls, filename_stem: str) -> str:
        """
        Given a particular file name (referred to as "stem", because we'll be suffixing a numeric index),
        locate the first index such that no file with that name is present in the same directory as the original file.

        EXAMPLES:
            given "test.htm", it might return "test1.htm" or "test23.htm"
            given "/my_path/test.htm", it might return "/my_path/test1.htm" or "/my_path/test34.htm"
            given "myfile", it might return "myfile1" or "myfile8"

        :param filename_stem:   A filename with or without a path.  If no path is given, the current directory is understood

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
            print(f"`{filename}` already exists")
            index += 1      # Keep increasing the index in the filename, until an unused filename is found
            filename = basename + str(index) + suffix


        filename = basename + str(index) + suffix

        return filename
