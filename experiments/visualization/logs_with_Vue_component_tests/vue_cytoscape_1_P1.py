# Python test "P1" of using the graphic module (Vue component) "vue_cytoscape_1" in an HTML log

from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

COMPONENT_NAME = "vue_cytoscape_1"       # CHANGE AS NEEDED
LOG_FILENAME= "vue_cytoscape_1_P1.htm"   # CHANGE AS NEEDED


# Initialize the HTML logging.
# Note: the relative path is from the location of *THE LOG FILE* to the project's home
GraphicLog.config(filename=LOG_FILENAME,
                  components=COMPONENT_NAME,
                  home_rel_path="../../..")


log.write(f'Example of Network Diagram, using the module "{COMPONENT_NAME}":',
          style=log.h2, blanks_after=1)


# EXAMPLE OF DATA VISUALIZATION (with arbitrary data, hardwired below)

all_data = {
#Data to define the nodes and edges of the network
    'graph':
    [
        {'id': 1, 'label': 'Reactant', 'name': 'A', 'diff_rate': 0.4},
        {'id': 2, 'label': 'Reactant', 'name': 'B', 'diff_rate': 1.8},
        {'id': 3, 'label': 'Reactant', 'name': 'C', 'diff_rate': 0.57},

        {'id': 4, 'label': 'Product' , 'name': 'X', 'diff_rate': 3.5},
        {'id': 5, 'label': 'Product' , 'name': 'Y', 'diff_rate': 0.1},

        {'id': 6, 'label': 'Reaction' , 'name': 'RXN', 'Rf': 24.2, 'Rb': 4.1},

        {'id': 10, 'source': 1, 'target': 6, 'name': 'reacts'},
        {'id': 11, 'source': 2, 'target': 6, 'name': 'reacts'},
        {'id': 12, 'source': 3, 'target': 6, 'name': 'reacts'},

        {'id': 13, 'source': 6, 'target': 4, 'name': 'produces'},
        {'id': 14, 'source': 6, 'target': 5, 'name': 'produces'}
    ],

    # Mapping the node label to its interior color
    'color_mapping':  {
        'Reactant': 'neo4j_green',
        'Product': 'neo4j_red',
        'Reaction': 'neo4j_lightbrown'
    }
}

# Send the plot to the HTML log file
GraphicLog.export_plot(all_data, COMPONENT_NAME)
