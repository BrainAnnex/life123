from typing import Union, List
from src.modules.html_log.html_log import HtmlLog as log


class GraphicLog:
    """
    To simplify use of HtmlLog and Vue components from within this project
    """

    ###  START OF CONFIGURATION  ###
    LOCAL_FILES_ROOT = None     # This needs to be set by config().   EXAMPLE: "../../../"
                                #   For local files (NOT recommended), it's the relative path from location of *THE LOG FILE* to Life123's home
                                #   Using local files causes issues in JupyterLab
    LOCAL_VUE_LIBRARY_FILE = "modules/Vue2_lib/vue2.js"                 # Relative to LOCAL_FILES_ROOT
    LOCAL_SVG_HELPER_LIBRARY_FILE = "modules/SVG_helper/svg_helper.js"  # Relative to LOCAL_FILES_ROOT
    LOCAL_VUE_COMPS_DIR = "modules/visualization/vue_components/"       # Relative to LOCAL_FILES_ROOT

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
    def config(cls, filename, components: Union[str, List[str]], mode='overwrite', extra_js=None, local_files=False, home_rel_path=None) -> None:
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
            assert home_rel_path, "config(): if the local_files=True flag is set, then an argument must be passed for `home_rel_path`"

            cls.FILES_ROOT = home_rel_path
            cls.VUE_LIBRARY_FILE = home_rel_path + "/" + cls.LOCAL_VUE_LIBRARY_FILE
            cls.SVG_HELPER_LIBRARY_FILE = home_rel_path + "/" + cls.LOCAL_SVG_HELPER_LIBRARY_FILE
            cls.VUE_COMPS_DIR = home_rel_path + "/" + cls.LOCAL_VUE_COMPS_DIR
        else:
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
    def export_plot(cls, data: dict, graphic_component: str, print_notification=True) -> None:
        """
        Send a Vue-based plot to the log file

        :param data:                A python dictionary of data to pass to the Vue component
        :param graphic_component:   A string with the name of the existing Vue.js component to use.
                                        EXAMPLE: "vue_curves_4" (assuming that a js file with such a component exists)
        :param print_notification:  If True, something is printed to inform of what's happening with the log file
        :return:                    None
        """
        assert cls.FILES_ROOT is not None, "export_plot(): must first initialize library with call to GraphicLog.config()"


        log.export_plot_Vue(data=data,
                            component_name = graphic_component,
                            component_file = f"{cls.VUE_COMPS_DIR}{graphic_component}.js")

        if print_notification:
            print(f"[GRAPHIC ELEMENT SENT TO LOG FILE `{log.log_fullname}`]")
