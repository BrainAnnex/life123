# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

import ipynbname

from life123 import check_version, ChemData, UniformCompartment, PlotlyHelper, GraphicLog, DisplayNetwork

# %%
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Specify the chemicals and the reactions; this data structure will get re-used in 
# the various simulations below
chem_data = ChemData()

uc = UniformCompartment(chem_data=chem_data, preset="slower")

# Reaction A <-> B (fast)
uc.add_reaction(reactants="A", products="B",
                       kF=64., kR=8.) 

# Reaction B <-> C (slow)
uc.add_reaction(reactants="B", products="C",
                       kF=12., kR=2.) 

print("Number of reactions: ", uc.number_of_reactions())

# %%
uc.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
rnxs = uc.get_reactions()

rnxs.plot_reaction_network("vue_cytoscape_2")

# %%

# %%
rxn_graph_data = rnxs.prepare_graph_network()
rxn_graph_data

# %%

# %%
structure = [{'id': 'n1', 'labels': ['label_1'], 'name': 'A', 'field_1': 123},
             {'id': 'n2', 'labels': ['label_2'], 'name': 'B', 'field_1': 8},

             {'id': 'edge-1', 'name': 'produces', 'source': 'n1', 'target': 'n2', 'link_prop': "hello"}
           ]

# %%
color_mapping = {'label_1': '#8DCC92', 'label_2': '#D9C8AD'}
caption_mapping = {'label_1': 'name', 'label_2': 'name'}

# %%
fig = {"structure": structure, "color_mapping": color_mapping, "caption_mapping": caption_mapping}
fig

# %%

# %%
# Send a plot of the network of reactions to the HTML log file
GraphicLog.export_plot(graph_data=fig, graphic_component="vue_cytoscape_2", unpack=False)

# %%

# %%
simplified = {'structure': 
              [
                  {'id': 'n-0',
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
                   'source': 'n-0',
                   'target': 'n-1',
                   'id': 'edge-1',
                   'edge_prop': "some value for the link"
                  }
              ],
              
             'color_mapping': {'label_1': '#8DCC92', 'label_2': '#D9C8AD'},
             'caption_mapping': {'label_1': 'name', 'label_2': 'name'}
             }

# %%

# %%
GraphicLog.export_plot(simplified, "vue_cytoscape_2", unpack=False)


# %%

# %%

# %%

# %%
class SimpleGraphic:
    
    @classmethod
    def _write_to_file(cls, file_handler, text :str) -> None:
        file_handler.write(text)
        file_handler.flush()    # To avoid weird buffering issues seen in JupyterLab
        
      
    
    @classmethod  
    def _html_header(cls) -> str:
        return '''<!-- Auto-generated log file -->
<!DOCTYPE html>
<html lang="en">
<head>  
    <meta charset="UTF-8">
    <title>Cytoscape graphics</title>

    <link type="text/css" rel="stylesheet" href="https://life123.science/libraries/vue_components/vue_cytoscape_2.css">

    <script src="https://life123.science/libraries/vue_2.6.12.js"></script>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <script src="https://life123.science/libraries/svg_helper_1.2.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js"></script>
</head>

'''

    
    @classmethod  
    def _html_body_start(cls, caption=None) -> str:
        return f'''
<body>

{caption}
        
'''

    
    @classmethod 
    def _html_vue_container(cls, vue_component :str, component_file :str, graph_data :dict, vue_count=1) -> str:
        """
        Generate and return the HTML for a Vue container (a DIV element with the Vue ROOT component),
        plus a script to instantiate the above Vue root component
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
            graph_data_json: {json.dumps(graph_data, indent=4)}
        }}  // END of data

    }});
</script>
        '''
    
    
    @classmethod
    def export_plot(cls, graph_data :dict, graphic_component :str, filename :str) -> None:
        """
        Send to the log file the data to create a Vue-based plot.
        This is meant to work alongside Vue components that expects 2 arguments ("props"):
            1) graph_data
            2) component_id

        :param graph_data:          A python dictionary of data to pass to the Vue component
        :param graphic_component:   A string with the name of the existing Vue.js component to use.
                                        EXAMPLE: "vue_curves_4" (assuming that a js file with such a Vue component exists)
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
       
        # TODO: make sure that `filename` ends with ".htm"

        
        # Prepare writing into the file (OVERWRITE)
        
        file_handler = open(filename, "w")      # Create a new file, to write to; over-write if present
                
        
        # Export into the HTML log file the various Vue-related parts
        
        VUE_COMPS_DIR = "https://life123.science/libraries/vue_components/"
        component_file = f"{VUE_COMPS_DIR}{graphic_component}.js"
        
        html = cls._html_header() + cls._html_body_start(caption="<h1>Interactive graph plot</h1>") + \
               cls._html_vue_container(vue_component=graphic_component, vue_count=1, component_file=component_file, graph_data=graph_data)
        
        cls._write_to_file(file_handler, text = html)
             
                           
        print(f"[GRAPHIC ELEMENT SENT TO FILE `{filename}`]")

# %%

# %%
import json

# %%
SimpleGraphic.export_plot(graph_data=simplified, graphic_component="vue_cytoscape_2", filename="cascade_test.SIMPLIFIED.htm")

# %%

# %%
DisplayNetwork.export_plot(graph_data=simplified, graphic_component="vue_cytoscape_2", filename="cascade_my_test")

# %%
