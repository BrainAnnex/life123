# Examples of creating a log with D3 graphics

from life123 import HtmlLog as log

log.config(filename="logs/test_with_D3.htm", mode='overwrite', use_D3=True, css="../D3_barchart.css")

log.new_run()


###################   1st plot   ##########################

log.write("Time to do a D3.js plot:", style=log.h1)

my_data = [
    {"id": 'd1', "index": 0, "value": 10, "region":'USA'},
    {"id": 'd2', "index": 1, "value": 11, "region":'Italy'},
    {"id": 'd3', "index": 2, "value": 12, "region":'Malta'},
    {"id": 'd4', "index": 3, "value": 6, "region":'Germany'}
]

handler_file = "../D3_barchart.js"
handler_func = "draw_barchart"

svg_id = "svg1"     # ID to use for the <SVG element>

log.export_plot_D3(data=my_data, svg_id=svg_id, js_file=handler_file, js_func=handler_func)


###################   2nd plot   ##########################

log.write("A repeat of the same plot:", style=log.h2, blanks_before=2)

svg_id = "svg2"     # ID to use for the <SVG element>

log.export_plot_D3(data=my_data, svg_id=svg_id, js_file=handler_file, js_func=handler_func)


###################   3rd plot   ##########################

my_heatmap_data = [
    { "group": "A", "variable": "v1", "value": "30" },
    { "group": "A", "variable": "v2", "value": "95" },
    { "group": "B", "variable": "v1", "value": "37" },
    { "group": "B", "variable": "v2", "value": "50" },
    { "group": "C", "variable": "v1", "value": "96" },
    { "group": "C", "variable": "v2", "value": "13" }
]

log.write("And now, a heatmap:", style=log.h2, blanks_before=2)

handler_file = "../D3_heatmap.js"
handler_func = "draw_heatmap"
svg_id = "svg3"                      # ID to use for the <DIV element>

log.export_plot_D3(data=my_heatmap_data, svg_id=svg_id, js_file=handler_file, js_func=handler_func, D3_tag="div")
