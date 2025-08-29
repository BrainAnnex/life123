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

# %% [markdown]
# ### Demonstration of file storage of system concentration history, and stop-restart of simulation, for the reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium.
#
# Same as experiment `react_1_a`, but with CSV file storage of the concentration history, and reaction stop-restart.

# %% [markdown]
# ### TAGS :   "basic", "uniform compartment"

# %%
LAST_REVISED = "Aug. 28, 2025"
LIFE123_VERSION = "1.0.0rc5"         # Library version this experiment is based on

# %%
#import set_path            # Using MyBinder?  Uncomment this before running the next cell!
                            # Importing this module will add the project's home directory to sys.path

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

import ipynbname
import life123

# %%
life123.check_version(LIFE123_VERSION)    # To check compatibility

# %%
# Initialize logging (for the system state)
csv_log_file = ipynbname.name() + "_system_log.csv"   # Use the notebook base filename 
                                                      # IN CASE OF PROBLEMS, set manually to any desired name

# %%

# %%

# %% [markdown]
# ## Initialize the Uniform-Compartment Simulation

# %%
# Instantiate the simulator and specify the chemicals
uc = life123.UniformCompartment()  

# %%
# We're now requesting that all System Concentration Data get logged in our previously-specified CSV file
uc.start_csv_log(csv_log_file)

# %%
# Reaction A <-> B , with 1st-order kinetics in both directions
uc.add_reaction(reactants="A", products="B", kF=3., kR=2.)

uc.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
uc.set_conc({"A": 80., "B": 10.})

# %%

# %% [markdown]
# #### This time (contrasted to experiment `react_1_a`) we'll be running the simulation in two parts

# %% [markdown]
# ## Part 1 (early run)

# %%
uc.single_compartment_react(initial_step=0.1, target_end_time=0.2)   # The first part of our run

# %%
uc.get_history()      # The system's concentrations history, saved during the run of single_compartment_react()

# %%
uc.plot_history(show_intervals=True)

# %%
# We're nowhere near equilibrium yet!
uc.is_in_equilibrium()

# %%

# %%

# %% [markdown]
# ## Part 2 (late run)

# %%
initial_step = 0.013759414272    # We're choosing this value simply FOR DEMONSTRATION PURPOSES,
                                 # to remain in exact lockstep with the time course of experiment `react_1_a`

'''
If you run experiment `react_1_a`, you can determine what the next time step would have been, had we not stopped early this time.

    In experiment `react_1_a`, after running the simulation, you can issue:

    list(uc.get_history(t_start=0.21, t_end=0.23, columns="SYSTEM TIME"))

    and you will get:  [0.211700785152, 0.225460199424]

    Their difference is the initial_step we're using here.
''';

# %%
uc.single_compartment_react(initial_step=0.013759414272, target_end_time=1.0)   # The 2nd part of our run, to the final target end time

# %%
uc.get_history()   # The cumulative system's history

# %%
uc.plot_history(show_intervals=True)

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %%

# %% [markdown]
# #### As we requested, a log of the concentration data, in CSV format, has ben saved in the following file:

# %%
csv_log_file

# %%
# Here's dump of the contents of that log file
with open(csv_log_file, 'r', encoding='utf8') as fh:
    file_contents = fh.read()
    print(file_contents)

# %%
