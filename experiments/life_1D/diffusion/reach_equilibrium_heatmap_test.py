"""
Exploring reaching equilibrium, first on a shorter timescale and then a longer one (but
both with identical time steps.)

The system starts out with a pulse in bin 2 (the 3rd bin from the left)

Notice the diffusing pulse "bouncing" off the left wall after total time 30
"""

from life_1D.bio_sim_1d import BioSim1D as bio
from modules.html_log.html_log import HtmlLog as log



bio.initialize_universe(n_bins=10, n_species=1)

bio.set_uniform_concentration(species_index=0, conc=0.)
bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)

bio.set_diffusion_rates([0.1])

bio.describe_state(show_diffusion_rates=True)



log.config(filename="../logs/reach_equilibrium.htm", overwrite=True,
           use_D3=True,
           Vue_lib = "../../../modules/Vue2_lib/vue2.js",
           js = "../../../modules/SVG_helper/svg_helper.js",
           css="../../../modules/visualization/D3_heatmap.css")
# Note: paths are from the location of THE LOG FILE


log.write("Dtime=10, with time steps of 0.1 ...", blanks_before=2)


total_time = 0.
for i in range(2):
    delta_time = 10.
    status = bio.diffuse(time_duration=delta_time, time_step=0.1)
    total_time += delta_time
    log.write(f"After Delta time {delta_time}.  TOTAL TIME {total_time}  ({status['steps']} steps taken):",
              blanks_before=1)
    bio.describe_state(concise=True)


log.blank_line()


# DATA VISUALIZATION WITH A HEATMAP

vue_id = "vue-root-1"     # Unique ID to use for the <DIV> containing the Vue component

all_data = {
    "my_groups": ["A", "B", "C"],
    "my_vars": ["v1", "v2"],
    "my_data": [
        { "group": "A", "variable": "v1", "value": "0" },
        { "group": "A", "variable": "v2", "value": "82" },
        { "group": "B", "variable": "v1", "value": "37" },
        { "group": "B", "variable": "v2", "value": "50" },
        { "group": "C", "variable": "v1", "value": "100" },
        { "group": "C", "variable": "v2", "value": "13" }
    ],
    "range_min": 0,
    "range_max": 100,
    "outer_width": 850,
    "outer_height": 450,
    "margins": {"top": 30, "right": 30, "bottom": 30, "left": 30}
}

log.export_plot_Vue(data=all_data, vue_id=vue_id,
                    component_name="vue-heatmap-9",
                    component_file="../../../modules/visualization/vue_components/heatmap9.js")



###################   A repeat of the same type of plot, with different data  ###########################

log.write("A repeat of the plot:", style=log.h1, newline=False)

vue_id = "vue-root-2"     # Unique ID to use for the <DIV> containing the Vue component


my_groups = [str(i) for i in range(bio.n_bins)]
print()
print(my_groups)

my_data = [{"group": str(i), "variable": "Mol 0", "value": str(bio.univ[0, i])}
           for i in range(bio.n_bins)]
print(my_data)

all_data = {
    "my_groups": my_groups,
    "my_vars": ["Mol 0"],
    "my_data": my_data,
    "range_min": 0,
    "range_max": 2.5,
    "outer_width": 850,
    "outer_height": 450,
    "margins": {"top": 30, "right": 30, "bottom": 30, "left": 30}
}

log.export_plot_Vue(data=all_data, vue_id=vue_id,
                    component_name="vue-heatmap-9",
                    component_file="../../../modules/visualization/vue_components/heatmap9.js")
