# To simplify use of HtmlLog and Vue components from within this project

from typing import Union, List
from modules.html_log.html_log import HtmlLog as log


class GraphicLog:

    #COMPONENT_NAME = None
    HOME_REL_PATH = None


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
    def config(cls, filename, components: Union[str, List[str]], home_rel_path, mode='overwrite'):
        """

        :param filename:
        :param components:      Either a string or list of strings,
                                with the name(s) of ALL the graphic components that will be used in the log file
        :param home_rel_path:   Relative path from location of the LOG file to Life123's home
                                EXAMPLE: "../../.."
        :param mode:
        :return:
        """
        if type(components) == str:
            css_files = f"{home_rel_path}/modules/visualization/vue_components/{components}.css"
        elif type(components) == list:
            css_files = [f"{home_rel_path}/modules/visualization/vue_components/{comp}.css"
                                        for comp in components]
        else:
            raise Exception("GraphicLog.config(): argument `components` must be either a string or a list of strings")


        # Note: paths are from the location of THE LOG FILE
        log.config(filename=filename, mode=mode,
                   use_D3=True,
                   Vue_lib = f"{home_rel_path}/modules/Vue2_lib/vue2.js",
                   js  = f"{home_rel_path}/modules/SVG_helper/svg_helper.js",
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
        #assert cls.COMPONENT_NAME is not None, "Must first call GraphicLog.config(), and pass a `component_name` argument"
        assert cls.HOME_REL_PATH is not None, "Must first call GraphicLog.config(), and pass a `home_rel_path` argument"


        log.export_plot_Vue(data=data,
                        component_name = graphic_component,
                        component_file = f"{cls.HOME_REL_PATH}/modules/visualization/vue_components/{graphic_component}.js")

        if print_notification:
            print("[Graphic element sent to log file]")
