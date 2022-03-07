"""
Exploring reaching equilibrium, including a heatmap depiction of some of the steps.

The system starts out with a pulse in bin 2 (the 3rd bin from the left)
"""

from modules.chemicals.chemicals import Chemicals as chem
from life_1D.bio_sim_1d import BioSim1D as bio
from modules.html_log.html_log import HtmlLog as log


chem_data = chem(diffusion_rates=[0.1])
bio.initialize_universe(n_bins=10, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=0.)
bio.inject_conc_to_cell(species_index=0, bin=2, delta_conc=10.)

bio.describe_state(show_diffusion_rates=True)



log.config(filename="../logs/reach_equilibrium.htm", mode='overwrite',
           use_D3=True,
           Vue_lib = "../../../modules/Vue2_lib/vue2.js",
           js = "../../../modules/SVG_helper/svg_helper.js",
           css="../../../modules/visualization/D3_heatmap.css")
# Note: paths are from the location of THE LOG FILE


log.write("1-D diffusion to equilibrium of a single species, with Diffusion rate 0.1. Time steps of 0.1",
          style=log.bold, blanks_before=2)
log.write("Heatmap with log scale in domain [0.5-10]", style=log.color, style_par='#BBB')



def visualize_state(time: float) -> None:
    """
    NOTE: a related function is now available as a method of BioSim1D

    :param time:
    :return:
    """
    log.write(f"Time : {time}", style=log.h1, newline=False)

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
        "range_max": 10,
        "outer_width": 850,
        "outer_height": 50,
        "margins": {"top": 10, "right": 30, "bottom": 18, "left": 30}
    }

    log.export_plot_Vue(data=all_data,
                        component_name="vue-heatmap-9",
                        component_file="../../../modules/visualization/vue_components/heatmap9.js")



#############################################


delta_time = 10.

total_time = 0.

visualize_state(total_time)

for i in range(50):
    status = bio.diffuse(time_duration=delta_time, time_step=0.1)
    total_time += delta_time

    print(f"\nAfter Delta time {delta_time}.  TOTAL TIME {total_time}  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)

    if i<2 or i==6 or i>=49:
        visualize_state(total_time)
