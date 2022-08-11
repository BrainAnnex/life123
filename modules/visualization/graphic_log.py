# To simplify use of HtmlLog and Vue components from within this project

from typing import Union, List
from modules.html_log.html_log import HtmlLog as log


class GraphicLog:

    #COMPONENT_NAME = None
    HOME_REL_PATH = None        # Relative path from location of *THE LOG FILE* to Life123's home
                                # EXAMPLE: "../../.."


    @classmethod
    def is_initialized(cls) -> bool:
        """
        Return True if this module has been initialized, or False otherwise

        :return:
        """
        #if cls.COMPONENT_NAME is None or cls.HOME_REL_PATH is None:
        if cls.HOME_REL_PATH is None:
            return False
        return True



    @classmethod
    def config(cls, filename, components: Union[str, List[str]], home_rel_path, mode='overwrite', extra_js=None):
        """

        :param filename:        Name, with full path, of the desired log file
        :param components:      Either a string or list of strings,
                                with the name(s) of ALL the graphic components that will be used in the log file
        :param home_rel_path:   Relative path from location of *THE LOG FILE* to Life123's home
                                EXAMPLE: "../../.."
        :param mode:            A string with the desired logging mode.  Must be one of:
                                    "overwrite" - if the log file exists, first clear it; otherwise, create it
                                    "multiple"  - each run's output should go into a separate file (consecutively numbered)
                                (note the "append" mode of the class HtmlLog cannot be used)
        :param extra_js:        Optional string with the name of an "extra" JavaScript file to include
                                ("extra" refers to the fact that several JS files are already automatically included)

        :return:
        """
        if type(components) == str:
            css_files = f"{home_rel_path}/modules/visualization/vue_components/{components}.css"
        elif type(components) == list:
            css_files = [f"{home_rel_path}/modules/visualization/vue_components/{comp}.css"
                                        for comp in components]
        else:
            raise Exception("GraphicLog.config(): argument `components` must be either a string or a list of strings")


        assert mode == "overwrite" or mode == "append", \
                "GraphicLog.config(): argument `mode` must be either 'overwrite' or 'append'"

        js = f"{home_rel_path}/modules/SVG_helper/svg_helper.js"
        if extra_js:
            js = [js, extra_js]

        # Note: paths are from the location of *THE LOG FILE*
        log.config(filename=filename, mode=mode,
                   use_D3=True,
                   Vue_lib = f"{home_rel_path}/modules/Vue2_lib/vue2.js",
                   js  = js,
                   css = css_files)

        #cls.COMPONENT_NAME = components
        cls.HOME_REL_PATH = home_rel_path



    @classmethod
    def export_plot(cls, data, graphic_component, print_notification=True) -> None:
        """

        :param data:
        :param graphic_component:
        :param print_notification:
        :return:                    None
        """
        assert cls.HOME_REL_PATH is not None, "Must first call GraphicLog.config(), and pass a `home_rel_path` argument"


        log.export_plot_Vue(data=data,
                        component_name = graphic_component,
                        component_file = f"{cls.HOME_REL_PATH}/modules/visualization/vue_components/{graphic_component}.js")

        if print_notification:
            print(f"[GRAPHIC ELEMENT SENT TO LOG FILE `{log.log_fullname}`]")
