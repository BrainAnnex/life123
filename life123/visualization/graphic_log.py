# Two classes, "DisplayNetwork" and "GraphicLog" (deprecated!)

from typing import Union, List
from life123.html_log import HtmlLog as log
import json




class DisplayNetwork:
    """
    Simplified version of the old classes "GraphicLog" and "HtmlLog".

    Used to create a single HTML file that graphically displays an interactive version
    of graph-network data
    """

    @classmethod
    def _write_to_file(cls, file_handler, text :str) -> None:
        """
        Write the given text (containing simple text and/or HTML code) into the file managed with
        the given File Handler

        :param file_handler:
        :param text:        String to write to the file managed by the above file handler
        :return:            None
        """
        file_handler.write(text)
        file_handler.flush()    # To avoid weird buffering issues seen in JupyterLab



    @classmethod
    def _html_header(cls, component_file_css :str) -> str:
        """
        Generate and return the text for the initial part of an HTML file, including the HEAD section.
        Various general JavaScript library files (Vue 2, D3 and Cytoscape) are for now hardwired in the code.

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

        :param caption:
        :return:        A string with HTML code
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

        :param vue_component:
        :param component_file:
        :param graph_data:
        :param vue_count:
        :return:                A string with HTML code
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
    // Instantiation of the Vue ROOT component
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
                    filename :str, caption="Interactive network plot") -> None:
        """
        Send to the given file the data to create a Vue-based display of a network.

        This is meant to work alongside a Vue.js component that expects 2 arguments ("props"):
            1) graph_data
            2) component_id

        :param graph_data:          A python dictionary of data to pass to the Vue component.
                                        It must contain 3 keys:
                                            'structure'         A list of dicts
                                                EXAMPLE, representing two nodes and an edge between them:
                                                    [     {'id': 'n-0',
                                                           'labels': ['label_1'],
                                                           'name': 'Company A',
                                                           'my_field': 1234
                                                          },
                                                          {'id': 'n-1',
                                                           'labels': ['label_2'],
                                                           'name': 'Customer B',
                                                           'my_field': 88
                                                          },
                                                          {'name': 'SELLS_TO',
                                                           'id': 'edge-1',
                                                           'source': 'n-0',
                                                           'target': 'n-1',
                                                           'edge_prop': 'some value'
                                                          }
                                                    ]

                                            'color_mapping'     A dict
                                                EXAMPLE: {"label_1": "#8DCC92", "label_2": "#D9C8AD"}

                                            'caption_mapping'   A dict
                                                EXAMPLE: {"label_1": "name", "label_2": "name"}

        :param graphic_component:   A string with the name of the existing Vue.js component to use.
                                        EXAMPLE: "vue_cytoscape_2" (assuming that a .js file with such a Vue component exists)

        :param filename:            The name of the file into which to create the interactive network plot.
                                        The suffix ".htm" will be added if it doesn't end with ".htm" or ".html"

        :param caption:             String displayed at the top of the HTML file

        :return:                    None
        """
        # Perform data validation
        assert type(graph_data) == dict, "export_plot(): argument `graph_data` must be a python dictionary"

        assert len(graph_data) == 3, \
                "export_plot(): argument `graph_data` must contains exactly 3 keys, named 'structure', 'color_mapping', 'caption_mapping'"

        assert ('structure' in graph_data) and type(graph_data['structure']) == list, \
                f"export_plot(): the argument `graph_data` must contain a key named 'structure' whose value is a list.  Passed value: {graph_data.get('structure')}"

        assert ('color_mapping' in graph_data) and type(graph_data['color_mapping']) == dict, \
                f"export_plot(): the argument `graph_data` must contain a key named 'color_mapping' whose value is a python dictionary.  Passed value: {graph_data.get('color_mapping')}"

        assert ('caption_mapping' in graph_data) and type(graph_data['caption_mapping']) == dict, \
                f"export_plot(): the argument `graph_data` must contain a key named 'caption_mapping' whose value is a python dictionary.  Passed value: {graph_data.get('caption_mapping')}"

        assert type(filename) == str, "export_plot(): argument `filename` must be a string"

        # Add the suffix ".htm" to `filename`, unless it already end with ".htm" or ".html"
        if not filename.lower().endswith((".htm", ".html")):
            filename += ".htm"


        # Prepare writing into the file (OVERWRITE)

        file_handler = open(filename, "w")      # Create a new file, to write to; over-write if present


        # Export into the HTML file the various Vue-related parts

        VUE_COMPS_DIR = "https://life123.science/libraries/vue_components/"
        component_file_js = f"{VUE_COMPS_DIR}{graphic_component}.js"
        component_file_css = f"{VUE_COMPS_DIR}{graphic_component}.css"

        html = cls._html_header(component_file_css) + cls._html_body_start(caption=f"<h1>{caption}</h1>") + \
               cls._html_vue_container(vue_component=graphic_component, vue_count=1, component_file=component_file_js, graph_data=graph_data)

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
                "config(): if the local_files=Flase flag is set, then the `home_rel_path` must be None"

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
                                        EXAMPLE: "vue_curves_4" (assuming that a js file with such a Vue component exists)
        :param print_notification:  If True, something is printed to inform of what's happening with the log file
        :param unpack:              Use True for Vue components that require their data unpacked into individual arguments;
                                        False for that accept a single data argument, named "graph_data"
        :return:                    None
        """
        assert cls.FILES_ROOT is not None, \
            "export_plot(): must first initialize library with call to GraphicLog.config()"

        # Perform data validation
        assert type(graph_data) == dict, "export_plot(): argument `graph_data` must be a python dictionary"

        assert len(graph_data) == 3, \
                "export_plot(): argument `graph_data` must contains exactly 3 keys, named 'structure', 'color_mapping', 'caption_mapping'"

        assert ('structure' in graph_data) and type(graph_data['structure']) == list, \
                f"export_plot(): the argument `graph_data` must contain a key named 'structure' whose value is a list.  Passed value: {graph_data.get('structure')}"

        assert ('color_mapping' in graph_data) and type(graph_data['color_mapping']) == dict, \
                f"export_plot(): the argument `graph_data` must contain a key named 'color_mapping' whose value is a python dictionary.  Passed value: {graph_data.get('color_mapping')}"

        assert ('caption_mapping' in graph_data) and type(graph_data['caption_mapping']) == dict, \
                f"export_plot(): the argument `graph_data` must contain a key named 'caption_mapping' whose value is a python dictionary.  Passed value: {graph_data.get('caption_mapping')}"


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

