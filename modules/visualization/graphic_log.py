# EXPERIMENTAL - to simplify use of HtmlLog and Vue components from within this project

from modules.html_log.html_log import HtmlLog as log


class GraphicLog:

    COMPONENT_NAME = None
    HOME_REL_PATH = None



    @classmethod
    def config(cls, filename, component_name, home_rel_path, mode='overwrite'):
        """

        :param filename:
        :param component_name:
        :param home_rel_path:   Relative path from location of the LOG file to Life123's home
                                EXAMPLE: "../../.."
        :param mode:
        :return:
        """
        # Note: paths are from the location of THE LOG FILE
        log.config(filename=filename, mode=mode,
                   use_D3=True,
                   Vue_lib = f"{home_rel_path}/modules/Vue2_lib/vue2.js",
                   js  = f"{home_rel_path}/modules/SVG_helper/svg_helper.js",
                   css = f"{home_rel_path}/modules/visualization/vue_components/{component_name}.css")

        cls.COMPONENT_NAME = component_name
        cls.HOME_REL_PATH = home_rel_path



    @classmethod
    def export_plot(cls, data):
        """

        :param data:
        :return:
        """
        assert cls.COMPONENT_NAME is not None, "Must first call GraphicLog.config(), and pass a `component_name` argument"
        assert cls.HOME_REL_PATH is not None, "Must first call GraphicLog.config(), and pass a `home_rel_path` argument"


        log.export_plot_Vue(data=data,
                        component_name = cls.COMPONENT_NAME,
                        component_file = f"{cls.HOME_REL_PATH}/modules/visualization/vue_components/{cls.COMPONENT_NAME}.js")
