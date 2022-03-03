"""
Exploring reaching equilibrium, including a heatmap depiction of some of the steps.

The system starts out with a pulse in bins near the left and the right endpoints
"""

from life_1D.bio_sim_1d import BioSim1D as bio
from modules.html_log.html_log import HtmlLog as log



# Initialize the system
bio.initialize_universe(n_species=1, n_bins=9)

bio.set_uniform_concentration(species_index=0, conc=0.)

bio.inject_conc_to_cell(species_index=0, bin=2, delta_conc=10.)
bio.inject_conc_to_cell(species_index=0, bin=6, delta_conc=10.)

bio.set_diffusion_rates([0.1])

bio.describe_state(show_diffusion_rates=True)



log.config(filename="../logs/reach_equilibrium.htm", mode='overwrite',
           use_D3=True,
           Vue_lib = "../../../modules/Vue2_lib/vue2.js",
           js = "../../../modules/SVG_helper/svg_helper.js",
           css="../../../modules/visualization/D3_heatmap.css")
# Note: paths are from the location of THE LOG FILE




delta_time = 3.

total_time = 0.

heatmap_pars = {"range": [1.5, 3.5],
                "outer_width": 850, "outer_height": 50,
                "margins": {"top": 10, "right": 30, "bottom": 18, "left": 30}
                }

log.write("1-D diffusion to equilibrium of a single species, with Diffusion rate 0.1. Time steps of 0.1",
          style=log.bold, blanks_before=2)
log.write(f"Heatmap with log scale in domain {heatmap_pars['range']}", style=log.color, style_par='#BBB')


bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {total_time}")

for i in range(15):
    status = bio.diffuse(time_duration=delta_time, time_step=0.1)
    total_time += delta_time

    print(f"\nAfter Delta time {delta_time}.  TOTAL TIME {total_time}  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)

    #if i<2 or i==6 or i>=14:
    #visualize_state(total_time)
    bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {total_time}")
