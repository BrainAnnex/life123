# Two classes, "DisplayNetwork" and "GraphicLog" (deprecated!)

from typing import Union, List
from life123.html_log import HtmlLog as log
import json




class DisplayNetwork:
    """
    Simplified version of the old classes "GraphicLog" and "HtmlLog".

    Used to create an HTML file that graphically displays an interactive version
    of graph-network data.
    This HTML file is a scaffold for one or move Vue components that provide the desired functionality
    """

    @classmethod
    def _write_to_file(cls, file_handler, text :str) -> None:
        """
        Write the given text (containing simple text and/or HTML code) into the file managed with
        the given File Handler

        :param file_handler:A file handler, for example as returned by calls to open()
        :param text:        String (simple text or HTML code) to write
                                to the file managed by the above file handler
        :return:            None
        """
        file_handler.write(text)
        file_handler.flush()    # To avoid weird buffering issues seen in JupyterLab



    @classmethod
    def _html_header(cls, component_file_css :str) -> str:
        """
        Generate and return the text for the initial part of an HTML file, including the HEAD section.
        Various general JavaScript library files (Vue v.2, D3 and Cytoscape)
        are for now hardwired in the code.

        :param component_file_css:  The full URL of the .CSS file used by our Vue component
        :return:                    A string with HTML code
        """
        return f'''<!-- Auto-generated log file -->
<!DOCTYPE html>
<html lang="en">
<head>  
    <meta charset="UTF-8">
    <title>Cytoscape graphics</title>

    <link type="text/css" rel="stylesheet" href="{component_file_css}">

    <script src="https://life123.science/libraries/vue_2.6.12.js"></script>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js"></script>
</head>

'''


    @classmethod
    def _html_body_start(cls, caption=None) -> str:
        """
        Prepare and return the text for the early part (but post-headers) of an HTML file,
        with an optional text caption to be shown at the top

        :param caption: [OPTIONAL] A string with plain text or HTML code.
                            EXAMPLES:  "Figure 1" ,  "<b>My Header</b>"
        :return:        A string with the HTML code to incorporate into the HTML file being formed
        """
        return f'''
<body>

{caption}
        
'''


    @classmethod
    def _html_vue_container(cls, vue_component :str, component_file :str, graph_data :dict, vue_count=1) -> str:
        """
        Generate and return the HTML for a Vue container (a DIV element with the Vue ROOT component),
        plus a script to instantiate the above Vue root component

        :param vue_component:   A string with the name of the existing Vue.js component to use.
                                    EXAMPLE: "vue_cytoscape_3" (assuming that a .js file with such a Vue component exists)
        :param component_file:  The full name (including path) of the .js file that contains the above component
        :param graph_data:      A python dictionary of data to pass to the Vue component.
                                    It must contain 3 keys: 'structure', 'color_mapping', 'caption_mapping'.
                                    For more details, see export_plot()
        :param vue_count:       An integer used to differentiate between multiple Vue components in the same HTML file.
                                    By default, 1
        :return:                A string with HTML code,
                                    including a SCRIPT element that instantiates the Vue root element
        """

        vue_id = f"vue-root-{vue_count}"    # EXAMPLE: "vue-root-1"

        return f'''

<div id="{vue_id}">   <!-- DIV container for the VUE COMPONENTS below : Vue ROOT element -->

    <{vue_component} v-bind:graph_data="graph_data_json"
                     v-bind:component_id="{vue_count}">        
    </{vue_component}>

</div>    <!--  ~~~~~~~~~~~~~~~~~~~~~  END of Vue root element  ~~~~~~~~~~~~~~~~~  -->



<!--
    Vue components (and other JS).  These must appear AFTER the Vue-containing elements
  -->
<script src="{component_file}"></script>



<script>
    // Instantiation of the Vue ROOT element
    new Vue({{
        el: '#{vue_id}',

        data: {{
            graph_data_json: 
{json.dumps(graph_data, indent=4)}
        }}  // END of data

    }});
</script>
        '''



    @classmethod
    def export_plot(cls, graph_data :dict, graphic_component :str,
                    filename :str, caption="<h1>Interactive network plot</h1>",
                    vue_comps_dir="https://life123.science/libraries/vue_components/") -> None:
        """
        Send to the given file the HTML data to create a Vue-based display of a network.

        This is meant to work alongside a Vue.js component that expects 2 arguments ("props"):
            1) graph_data
            2) component_id

        :param graph_data:          A python dictionary of data to pass to the Vue component.
                                        It must contain 4 keys:
                                            'nodes'         A list of dicts with node data
                                                EXAMPLE (two nodes):
                                                    [{'id': 1, '_node_labels': ['PERSON'], 'name': 'Julian'},
                                                     {'id': 2, '_node_labels': ['CAR'], 'color': 'white'}]

                                            'edges'         A list of dicts with edge data
                                                EXAMPLE (an edge between the two nodes of the previous example):
                                                    [{'id': 'edge-1', 'source': 1, 'target': 2, 'name': 'OWNS'}]

                                                Note: 'id' values can be integers or strings

                                            'color_mapping'     A dict
                                                EXAMPLE: {"label_1": "#8DCC92", "label_2": "#D9C8AD"}

                                            'caption_mapping'   A dict
                                                EXAMPLE: {"label_1": "name", "label_2": "name"}

        :param graphic_component:   The basename of a existing JavaScript and CSS files a
                                        that provides the interactive visualization functionality.
                                        The JS file is expected to implement a Vue.js component by the same name
                                        (but with hyphens in lieu of any underscore in the name, if applicable.)
                                        EXAMPLE: "vue_cytoscape_5" (assuming that a "vue_cytoscape_5.js" file
                                                 and a "vue_cytoscape_5.css" file
                                                 exist in the directory specified by the argument `vue_comps_dir`,
                                                 and that the .js file implements a Vue component named "vue-cytoscape-5")

        :param filename:            The name of the file into which to place the HTML code
                                        to create the interactive network plot.
                                        The suffix ".htm" will be added if it doesn't end with ".htm" or ".html"
                                        If the file already exists, it will get overwritten.
                                        (Note: this file will automatically include an internal reference to the JavaScript
                                        file specified in `graphic_component`)

        :param caption:             [OPTIONAL] String (plain text or HTML) displayed at the top of the HTML file.
                                        By default, a prominent "Interactive network plot"

        :param vue_comps_dir:       [OPTIONAL] The full name of the directory where the Vue components reside.
                                        A final "/" in the name is optional.
                                        Note: when this function is used in a Jupyterlab notebook, it's best to use
                                              a URL.  EXAMPLE: "https://life123.science/libraries/vue_components/"
        :return:                    None
        """
        # TODO: expand pytest

        # Perform data validation
        assert type(graph_data) == dict, "export_plot(): argument `graph_data` must be a python dictionary"

        assert len(graph_data) == 4, \
                f"export_plot(): argument `graph_data` must contains exactly 4 keys, " \
                f"named 'nodes', 'edges', 'color_mapping', 'caption_mapping'.  " \
                f"Instead, it contains the following keys: {list(graph_data)}"

        assert ('nodes' in graph_data) and type(graph_data['nodes']) == list, \
                f"export_plot(): the argument `graph_data` must contain a key named 'nodes' whose value is a list.  " \
                f"Passed value: {graph_data.get('nodes')}"

        assert ('edges' in graph_data) and type(graph_data['edges']) == list, \
                f"export_plot(): the argument `graph_data` must contain a key named 'edges' whose value is a list.  " \
                f"Passed value: {graph_data.get('edges')}"

        assert ('color_mapping' in graph_data) and type(graph_data['color_mapping']) == dict, \
                f"export_plot(): the argument `graph_data` must contain a key named 'color_mapping' whose value is a python dictionary.  " \
                f"Passed value: {graph_data.get('color_mapping')}"

        assert ('caption_mapping' in graph_data) and type(graph_data['caption_mapping']) == dict, \
                f"export_plot(): the argument `graph_data` must contain a key named 'caption_mapping' whose value is a python dictionary.  " \
                f"Passed value: {graph_data.get('caption_mapping')}"

        assert type(filename) == str, "export_plot(): argument `filename` must be a string"

        # Add the suffix ".htm" to `filename`, unless it already end with ".htm" or ".html"
        if not filename.lower().endswith((".htm", ".html")):
            filename += ".htm"


        # Prepare writing into the file (OVERWRITE)
        file_handler = open(filename, "w")      # Create a new file, to write to; over-write if present


        # Export into the HTML file the various Vue-related parts
        if not vue_comps_dir.endswith("/"):
            vue_comps_dir += "/"        # Add the final slash, unless already present

        assert not graphic_component.endswith(".js") and not graphic_component.endswith(".css"), \
            "export_plot(): the argument `graphic_component` should be the BASE NAME of existing .js and .css files; " \
            "do not include any suffix"

        component_file_js = f"{vue_comps_dir}{graphic_component}.js"
        component_file_css = f"{vue_comps_dir}{graphic_component}.css"

        # By convention, Vue components use hyphens instead of underscores in their name
        vue_component_name = graphic_component.replace("_", "-")    # Replace each instance of "_" (if any) with "-"

        html = cls._html_header(component_file_css) + cls._html_body_start(caption=caption) + \
               cls._html_vue_container(vue_component=vue_component_name, vue_count=1,
                                       component_file=component_file_js, graph_data=graph_data)

        cls._write_to_file(file_handler, text = html)

        print(f"[GRAPHIC ELEMENT SENT TO FILE `{filename}`]")







#################################################################################################

class GraphicLog:
    """
    DEPRECATED!

    To simplify use of HtmlLog and Vue components from within this project.
    """

    ###  START OF CONFIGURATION  ###
    #  NOTE: Using local files causes issues in JupyterLab, and it's NOT recommended!
    LOCAL_FILES_ROOT = None     # The project's root folder.  This needs to be set by config().   EXAMPLE: "../../../"
                                # It's the relative path from location of *THE LOG FILE* to Life123's home
    LOCAL_VUE_LIBRARY_FILE = "vue2_lib/vue2.js"                 # Relative to LOCAL_FILES_ROOT
    LOCAL_SVG_HELPER_LIBRARY_FILE = "SVG_helper/svg_helper.js"  # Relative to LOCAL_FILES_ROOT
    LOCAL_VUE_COMPS_DIR = "vue_components/"                     # Relative to LOCAL_FILES_ROOT

    # NOTE: Using remote files is recommended, especially when using JupyterLab
    REMOTE_FILES_ROOT = "https://life123.science/libraries/"
    REMOTE_VUE_LIBRARY_FILE = "https://life123.science/libraries/vue_2.6.12.js"
    REMOTE_SVG_HELPER_LIBRARY_FILE = "https://life123.science/libraries/svg_helper_1.2.js"
    REMOTE_VUE_COMPS_DIR = "https://life123.science/libraries/vue_components/"

    FILES_ROOT = None
    VUE_LIBRARY_FILE = None
    SVG_HELPER_LIBRARY_FILE = None
    VUE_COMPS_DIR = None
    ###  End of configuration  ###


    @classmethod
    def is_initialized(cls) -> bool:
        """
        Return True if this module has been initialized, or False otherwise

        :return:
        """
        if cls.FILES_ROOT is None:    
            return False

        return True



    @classmethod
    def config(cls, filename, components: Union[str, List[str]], mode='overwrite', extra_js=None,
                    local_files=False, home_rel_path=None) -> None:
        """
        Initialize this library
        
        :param filename:        Name, with FULL path, of the desired log file
        :param components:      Either a string or list of strings,
                                with the name(s) of ALL the Vue.js graphic components that will be used in the log file
        :param mode:            A string with the desired logging mode.  Must be one of:
                                    "overwrite" - if the log file exists, first clear it; otherwise, create it
                                    "multiple"  - each run's output should go into a separate file (consecutively numbered)
                                    (note the "append" mode of the class HtmlLog cannot be used)
        :param extra_js:        Optional string with the name of an "extra" JavaScript file to include
                                    ("extra" refers to the fact that several JS files are already automatically included)
                                
        :param local_files:     If True, local files will be used for the Vue components (see arg "home_rel_path")
        :param home_rel_path:   Relative path from location of *THE LOG FILE* to the Life123's project home
                                    NEW: ONLY APPLICABLE IF local_files=True ; otherwise, it's simply ignored
                                    EXAMPLE: "../../.."

        :return:                None
        """
        if local_files:
            assert home_rel_path, \
                "config(): if the local_files=True flag is set, then an argument must be passed for `home_rel_path`"

            cls.FILES_ROOT = home_rel_path
            cls.VUE_LIBRARY_FILE = home_rel_path + "/" + cls.LOCAL_VUE_LIBRARY_FILE
            cls.SVG_HELPER_LIBRARY_FILE = home_rel_path + "/" + cls.LOCAL_SVG_HELPER_LIBRARY_FILE
            cls.VUE_COMPS_DIR = home_rel_path + "/" + cls.LOCAL_VUE_COMPS_DIR
        else:
            assert not home_rel_path, \
                "config(): if the local_files=False flag is set, then the `home_rel_path` must be None"

            cls.FILES_ROOT = cls.REMOTE_FILES_ROOT
            cls.VUE_LIBRARY_FILE = cls.REMOTE_VUE_LIBRARY_FILE
            cls.SVG_HELPER_LIBRARY_FILE = cls.REMOTE_SVG_HELPER_LIBRARY_FILE
            cls.VUE_COMPS_DIR = cls.REMOTE_VUE_COMPS_DIR


        if type(components) == str:
            css_files = f"{cls.VUE_COMPS_DIR}{components}.css"
        elif type(components) == list:
            css_files = [f"{cls.VUE_COMPS_DIR}{comp}.css"
                                        for comp in components]
        else:
            raise Exception("GraphicLog.config(): argument `components` must be either a string or a list of strings")


        assert mode == "overwrite" or mode == "append", \
                "GraphicLog.config(): argument `mode` must be either 'overwrite' or 'append'"

        # Assemble a list of all the needed JavaScript files
        js = cls.SVG_HELPER_LIBRARY_FILE
        if extra_js:
            js = [js, extra_js]

        # Note: paths are from the location of *THE LOG FILE*
        log.config(filename=filename, mode=mode,
                   use_D3=True,
                   Vue_lib = cls.VUE_LIBRARY_FILE,
                   js  = js,
                   css = css_files)



    @classmethod
    def export_plot(cls, graph_data: dict, graphic_component: str, print_notification=True, unpack=False) -> None:
        """
        Send to the log file the data to create a Vue-based plot

        :param graph_data:          A python dictionary of data to pass to the Vue component
        :param graphic_component:   A string with the name of the existing Vue.js component to use.
                                        EXAMPLES (assuming that a js file with such a Vue component exists)
                                            "vue-cytoscape-5"
                                            "vue_curves_4" (OLD naming convention)
        :param print_notification:  If True, something is printed to inform of what's happening with the log file
        :param unpack:              Use True for Vue components that require their data unpacked into individual arguments;
                                        False for that accept a single data argument, named "graph_data"
        :return:                    None
        """
        assert cls.FILES_ROOT is not None, \
            "export_plot(): must first initialize library with call to GraphicLog.config()"

        # Perform data validation
        assert type(graph_data) == dict, "export_plot(): argument `graph_data` must be a python dictionary"

        assert len(graph_data) == 3 or len(graph_data) == 4, \
                    "export_plot(): argument `graph_data` must contains exactly 3 or 4 keys"
                    # TODO: Change the above error message.  Change to 4 keys named 'nodes', 'edges', 'color_mapping', 'caption_mapping'

        #assert ('structure' in graph_data) and type(graph_data['structure']) == list, \
                #f"export_plot(): the argument `graph_data` must contain a key named 'structure' whose value is a list.  Passed value: {graph_data.get('structure')}"

        assert ('color_mapping' in graph_data) and type(graph_data['color_mapping']) == dict, \
                    f"export_plot(): the argument `graph_data` must contain a key named 'color_mapping' whose value is a python dictionary.  " \
                    f"Passed value: {graph_data.get('color_mapping')}"

        assert ('caption_mapping' in graph_data) and type(graph_data['caption_mapping']) == dict, \
                    f"export_plot(): the argument `graph_data` must contain a key named 'caption_mapping' whose value is a python dictionary.  " \
                    f"Passed value: {graph_data.get('caption_mapping')}"


        if unpack:
            log.export_plot_Vue_unpack_args(data=graph_data,
                                            component_name = graphic_component,
                                            component_file = f"{cls.VUE_COMPS_DIR}{graphic_component}.js")
        else:
            log.export_plot_Vue(graph_data=graph_data,
                                component_name = graphic_component,
                                component_file = f"{cls.VUE_COMPS_DIR}{graphic_component}.js")

        if print_notification:
            print(f"[GRAPHIC ELEMENT SENT TO LOG FILE `{log.log_fullname}`]")
