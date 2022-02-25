# Examples of creating a log with D3 graphics in a Vue environment


from modules.html_log.html_log import HtmlLog as log

log.config(filename="logs/test_with_D3_plus_Vue.htm", overwrite=True,
           use_D3=True,
           Vue_lib = "../../../Vue2_lib/vue2.js",
           css=["../D3_interactive_star.css", "../D3_heatmap.css"],
           js="../../../SVG_helper/svg_helper.js")
# Note: paths are from the location of THE LOG FILE

log.new_run()


########################   D3 Plot with Vue  ###############################

log.write("Time to do a D3 plot with Vue:", style=log.h1, newline=False)

vue_id = "vue-root-1"     # Unique ID to use for the <DIV> containing the Vue component

log.export_plot_Vue(data=None, vue_id=vue_id,
                    component_name="vue-interactive-star", component_file="../Vue_components/interactive_star.js")


###################   A repeat of the same D3 Plot with Vue  ###########################

log.write("A repeat of the same plot:", style=log.h1, newline=False)

vue_id = "vue-root-2"     # Unique ID to use for the <DIV> containing the Vue component

my_data = {"outer_radius": 23}

log.export_plot_Vue(data=my_data, vue_id=vue_id,
                    component_name="vue-interactive-star", component_file="../Vue_components/interactive_star.js")


########################   Yet another D3 Plot with Vue  ###############################

log.write("Another D3 plot with Vue:", style=log.h1, newline=False)

vue_id = "vue-root-3"     # Unique ID to use for the <DIV> containing the Vue component

all_data = {
                "my_groups": ["A", "B", "C"],
                "my_vars": ["v1", "v2"],
                "my_data": [
                                { "group": "A", "variable": "v1", "value": "30" },
                                { "group": "A", "variable": "v2", "value": "95" },
                                { "group": "B", "variable": "v1", "value": "37" },
                                { "group": "B", "variable": "v2", "value": "50" },
                                { "group": "C", "variable": "v1", "value": "96" },
                                { "group": "C", "variable": "v2", "value": "13" }
                           ],
                "outer_width": 850,
                "outer_height": 450,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 30}
            }

log.export_plot_Vue(data=all_data, vue_id=vue_id,
                    component_name="vue-heatmap-8", component_file="../Vue_components/heatmap8.js")
