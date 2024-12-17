import numpy as np
import pandas as pd
from typing import Union
from life123.collections import CollectionTabular



class Diagnostics:
    """
    For the management of reaction diagnostic data
    """

    def __init__(self, reactions):

        assert reactions is not None, \
            "Diagnostics class cannot be instantiated with a missing value for the argument `reactions`"

        self.reactions = reactions

        self.chem_data = reactions.get_chem_data()


        # TODO: maybe drop the "diagnostic_" from the names, or rename it to "historic_"
        self.diagnostic_conc_data = CollectionTabular(parameter_name="TIME")
                                        # An expanded version of the normal System History.
                                        #   Columns of the dataframes:
                                        #       'TIME' 	'A' 'B' ...  'caption'
                                        #
                                        #   Note: if an interval run is aborted, NO entry is created here
                                        #         (this approach DIFFERS from that of other diagnostic data)


        self.diagnostic_rxn_data = {}   # "Diagnostic reaction data", PER REACTION: a dict with as many entries as reactions.
                                        #   The keys are the reaction indices; the values are objects of type "MovieTabular",
                                        #   which contain Pandas dataframes with the following columns
                                        #   (referring to one specific reaction):
                                        #           'START TIME' , 'time_step' , 'aborted', 'Delta A' , 'Delta B' , ... , 'rate', 'caption'
                                        #
                                        #   Note:   - entries are always added, even if a time-step run is aborted
                                        #           - the various 'Delta concentrations' are for the chemicals involved in the reaction
                                        #             over the time interval that *STARTS* at the value in the "START TIME" column


        self.diagnostic_decisions_data = CollectionTabular(parameter_name="START_TIME")
                                        #   Columns of the dataframes:
                                        #       'START_TIME' 	'Delta A' 'Delta B' ...
                                        #               [plus, if applicable, other fields such as
                                        #               'action', 'norm_A', 'norm_B', 'step_factors']
                                        #
                                        #   Note: entries are always added, even if an interval run is aborted




    #####  1. diagnostic_rxn_data  #####

    def save_rxn_data(self, rxn_index :int, system_time, time_step,
                      increment_vector_single_rxn: Union[np.array, None],
                      aborted=False,
                      rate=None,
                      caption="") -> None:
        """
        Save up diagnostic data for 1 reaction, for a simulation step
        (by convention, regardless of whether the step is completed or aborted)

        :param rxn_index:                   The integer index (0-based) to identify the reaction of interest
        :param system_time:                 The START time of the reaction step
        :param time_step:                   The duration of the current simulation step
        :param increment_vector_single_rxn: A Numpy array of size equal to the total number of chemical species,
                                                containing the "delta concentrations" caused by this reaction
                                                for ALL the chemicals (whether involved in the reaction or not);
                                                it may be None if the current reaction step was aborted
                                                Note: diagnostic data is saved only for the chemicals involved in this reaction
        :param aborted:                     True is the current reaction step was aborted (i.e. will get repeated);
                                                False (default) for normal situations
        :param rate:                        [OPTIONAL] The value of the reaction rate (aka reaction velocity) at this step
        :param caption:                     [OPTIONAL] string to describe the snapshot
        :return:                            None
        """
        # Validate the reaction index
        self.reactions.assert_valid_rxn_index(rxn_index)

        # Sorted list of the indexes of all the chemicals participating in this reaction
        indexes = self.reactions.get_chemicals_indexes_in_reaction(rxn_index)

        # Validate increment_vector_single_rxn
        if increment_vector_single_rxn is not None:     # If the values are available (not the case in aborts)
            assert len(increment_vector_single_rxn) == self.chem_data.number_of_chemicals(), \
                f"save_diagnostic_rxn_data(): the length of the Numpy array increment_vector_single_rxn " \
                f"({len(increment_vector_single_rxn)}) doesn't match the overall number of chemical ({self.chem_data.number_of_chemicals()})"


        # Initialize a "MovieTabular" object for this reaction, if needed
        if rxn_index not in self.diagnostic_rxn_data:
            self.diagnostic_rxn_data[rxn_index] = CollectionTabular(parameter_name="START_TIME")


        data_snapshot = {"time_step": time_step, "aborted": aborted}      # Dict being prepared to add a new row to a Pandas dataframe

        if increment_vector_single_rxn is None:     # If the values aren't available (as is the case in aborts)
            for index in indexes:       # For all the chemicals participating in this reaction
                data_snapshot["Delta " + self.chem_data.get_label(index)] = np.nan
        else:
            for index in indexes:       # For all the chemicals participating in this reaction
                data_snapshot["Delta " + self.chem_data.get_label(index)] = increment_vector_single_rxn[index]

        if rate is not None:
            data_snapshot["rate"] = rate

        self.diagnostic_rxn_data[rxn_index].store(par=system_time,
                                                  data_snapshot=data_snapshot, caption=caption)



    def save_diagnostic_aborted_rxns(self, system_time, time_step,
                                    caption="") -> None:
        """
        Save up diagnostic data, for ALL reactions, in case of aborted simulation step

        :param system_time: The START time of the reaction step
        :param time_step:   The duration of the current simulation step
        :param caption:     [OPTIONAL] string to describe the snapshot
        :return:            None
        """
        for rxn_index in range(self.reactions.number_of_reactions()):
            self.save_rxn_data(system_time=system_time, time_step=time_step,
                               increment_vector_single_rxn=None,
                               aborted=True,
                               rxn_index=rxn_index, caption=caption)



    def annotate_abort_rxn_data(self, msg: str) -> None:
        """
        Store the given value in the caption field of the last record
        for EACH of the reaction-specific dataframes.
        Also, set the "aborted" fields to "True"

        :param msg: Value to set the comment field to
        :return:    None
        """
        for rxn_index in range(self.chem_data.number_of_reactions()):
            self.diagnostic_rxn_data[rxn_index].set_caption_last_snapshot(msg)
            self.diagnostic_rxn_data[rxn_index].set_field_last_snapshot(field_name="aborted", field_value="True")
            # TODO: investigate why using "True" rather then True



    def get_system_history_with_rxn_rates(self, rxn_index :int) -> pd.DataFrame:
        """
        Return a Pandas dataframe containing the saved System History,
        plus an extra column named "rate", for the reaction rates of the specified reaction.
        Notice: the last row of the System History does NOT get included,
                because no rate information is known about the next simulation step not taken

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :return:            A Pandas data frame with the following columns:
                                "TIME"
                                A series of columns for each of the registered chemicals
                                "caption"
                                "rate"
        """
        # TODO: probably ditch, because no longer needed.  (Also, the dataframe merge might not always work...)
        # TODO: the times had better match up!  Currently not validated

        rates = self.get_rxn_rates(rxn_index=rxn_index)     # A Pandas dataframe with 2 columns: "START_TIME" and "rate"
        system_history = self.get_diagnostic_conc_data()    # Note that is a copy of the dataframe; so, no harm in changing it, below
        # Dropping the last row, because no rate information is known about the next simulation step not taken!
        system_history = system_history[:-1]

        # They'd better match up!
        assert len(system_history) == len(rates), \
            f"get_system_history_with_rxn_rates() : unable to reconcile " \
            f"the system history data with the reaction data for the requested reaction (index {rxn_index}). " \
            f"The diagnostic concentration data (excluding the last entry, which can't have a rate value) contains {len(system_history)} records, " \
            f"while the diagnostic rates data contains {len(rates)} values"

        system_history["rate"] =  rates["rate"]

        return system_history



    def get_rxn_rates(self, rxn_index :int) -> Union[pd.DataFrame, None]:
        """
        Return a Pandas dataframe with 2 columns: "START_TIME" and "rate",
        with the rates of the specified reaction at those times.

        Note that the dataframe will have 1 LESS ROW than a corresponding dataframe of concentration history,
        because the rate at the final time won't be known (since that time would be the start of a simulation step
        that has not been done)

        Also, aborted time steps are eliminated (only the rate of the final "re-do" is included)

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :return:            A Pandas dataframe with times and their corresponding reaction rates
        """
        rxn_df = self.get_rxn_data(rxn_index=rxn_index)
        if rxn_df is None:
            return None

        # Ditch any row where the "aborted" field has a True value; only take 2 columns of the resulting dataframe
        df_filtered = rxn_df[rxn_df['aborted'] == False][["START_TIME", "rate"]]

        return df_filtered.reset_index(drop=True)   # Re-number rows (drop=True avoids adding the old index as a new column)



    def get_rxn_data(self, rxn_index :int, head=None, tail=None,
                     t=None, print_reaction=True) -> Union[pd.DataFrame, None]:
        """
        Return a Pandas dataframe with the diagnostic run data of the requested SINGLE reaction,
        from the time that the diagnostics were enabled by instantiating this class.

        In particular, the dataframe contains the reaction rate at the start time,
        and the "Delta" values for each of the chemicals
        involved in the reaction - i.e. the change in their concentrations
        over the time interval that *STARTS* at the value in the "START_TIME" column.
        (So, there'll be no row with the final current System Time)

        Note: entries are always added, even if an interval run is aborted and automatically re-done;
              therefore, some entries will be duplicates in the "START_TIME" column.
              Those entries with have a value of True in the "aborted" column.

        Optionally, print out a brief description of the reaction.

        Optionally, limit the dataframe to a specified numbers of rows at the end,
        or just return one entry corresponding to a specific time
        (the row with the CLOSEST time to the requested one, which will appear in an extra column
        called "search_value")

        :param rxn_index:       The integer index (0-based) to identify the reaction of interest

        :param head:            (OPTIONAL) Number of records to return,
                                    from the start of the diagnostic dataframe.
        :param tail:            (OPTIONAL) Number of records to return,
                                    from the end of the diagnostic dataframe.
                                    If either the "head" arguments is passed, this argument will get ignored

        :param t:               (OPTIONAL) Individual time to pluck out from the dataframe;
                                    the row with closest time will be returned.
                                    If this parameter is specified, an extra column - called "search_value" -
                                    is inserted at the beginning of the dataframe
                                    If either the "head" or the "tail" arguments are passed, this argument will get ignored

        :param print_reaction:  (OPTIONAL) If True (default), concisely print out the requested reaction

        :return:                If present, return a Pandas data frame with (all or some of)
                                    the diagnostic data of the specified reaction.
                                    Columns of the dataframes:
                                    'START_TIME','time_step','aborted','Delta A','Delta B'...,'rate','caption'
                                If not present, return None
        """
        # Validate the reaction index
        self.reactions.assert_valid_rxn_index(rxn_index)

        if print_reaction:
            print("Reaction: ", self.reactions.single_reaction_describe(rxn_index=rxn_index, concise=True))

        movie_obj = self.diagnostic_rxn_data.get(rxn_index)    # Object of type "MovieTabular"

        if movie_obj is None:
            print(f"get_diagnostic_rxn_data(): no diagnostics data exists for reaction index {rxn_index} ; "
                  f"did you run the reaction simulation prior to calling this function?" )
            return None


        return movie_obj.get_dataframe(head=head, tail=tail, search_col="START_TIME", search_val=t, return_copy=True)




    #####  2. diagnostic_conc_data  #####

    def save_diagnostic_conc_data(self, system_data, system_time, caption="") -> None:
        """
        Save the diagnostic concentration data during the run, indexed by the current System Time.
        Note: if an interval run is aborted, by convention NO entry is created here

        :return: None
        """
        self.diagnostic_conc_data.store(par=system_time,
                                        data_snapshot=system_data, caption=caption)



    def get_diagnostic_conc_data(self) -> pd.DataFrame:
        """
        Return the diagnostic concentration data saved during the run.
        This will be a complete set of simulation steps,
        even if the user elected to only save part of the history during the run.

        Note: if the run of an interval step is aborted,
              by convention NO entry is created here

        :return: A Pandas dataframe, with the columns:
                    'TIME' 	'A' 'B' ... 'caption'
                    where 'A', 'B', ... are the labels of all the chemicals
        """
        return self.diagnostic_conc_data.get_dataframe(return_copy=True)




    #####  3. diagnostic_decisions_data  #####

    def save_diagnostic_decisions_data(self, system_time, data :dict,
                                       delta_conc_arr :Union[np.ndarray, None], caption="") -> None:
        """
        Used to save the diagnostic concentration data during the run, indexed by the given System Time.
        Note: if an interval run is aborted, by convention an entry is STILL created here

        :param system_time:
        :param data:
        :param delta_conc_arr:  A Numpy array of "delta concentrations".  EXAMPLE: array[1.23, 52.2]
        :param caption:         [OPTIONAL] String with a caption for this record
        :return:                None
        """
        delta_conc_dict = {}
        if delta_conc_arr is not None:
            delta_conc_dict = self._delta_conc_dict(delta_conc_arr)
            # EXAMPLE:  {"Delta A": 1.23, "Delta X": 52.2}


        delta_conc_dict.update(data)        # Merge the data dict into the delta_conc_dict

        self.diagnostic_decisions_data.store(par=system_time,
                                             data_snapshot=delta_conc_dict, caption=caption)



    def get_decisions_data(self) -> pd.DataFrame:
        """
        Determine and return the diagnostic data about concentration changes at every step - by convention,
        EVEN aborted ones
    
        :return:    A Pandas dataframe with a "TIME" column, and columns for all the "Delta concentration" values
        """
        return self.diagnostic_decisions_data.get_dataframe(return_copy=True)




    #############  EXPLAIN THINGS  #############


    def explain_reactions(self) -> bool:
        """
        Provide a detailed explanation of all the steps of the reactions,
        from the saved diagnostic data

        WARNING: Currently designed only for exactly 2 reactions!  TODO: generalize to any number of reactions

        TODO: test and validate usefulness, now that substeps got eliminated
        TODO: allow arguments to specify the min and max reaction times during which to display the explanatory data

        :return:    True if the diagnostic data is consistent for all the steps of all the reactions,
                    or False otherwise
        """
        number_reactions = self.chem_data.number_of_reactions()

        assert number_reactions == 2, \
            "explain_reactions() currently ONLY works when exactly 2 reactions are present. " \
            "Future versions will lift this restriction"

        row_baseline = 0
        row_list = [0, 0]       # TODO: generalize
        active_list = [0, 1]    # ALL the reactions.  TODO: generalize

        self._explain_reactions_helper(active_list=active_list,
                                       row_baseline=row_baseline, row_list=row_list)


        while row_baseline < len(self.get_diagnostic_conc_data()) - 2:
            row_baseline += 1

            print("Advance a single-step across all tables")

            for i in active_list:
                row_list[i] += 1

            status = self._explain_reactions_helper(active_list=active_list, row_baseline=row_baseline, row_list=row_list) # TODO: generalize


            if not status:
                return False    # Error termination

        # END while

        return True             # Successful termination



    def _explain_reactions_helper(self, active_list, row_baseline, row_list) -> bool:
        """
        Helper function for explain_reactions()

        :param active_list:
        :param row_baseline:
        :param row_list:
        :return:            True is the diagnostic data is consistent for this step, or False otherwise
        """
        print("-----------")
        print("ROW of baseline data: ", row_baseline)

        current_time = self.get_diagnostic_conc_data().loc[row_baseline]['TIME']
        print(f"TIME = {current_time:.5g}")

        print("row_list: ", row_list)
        print("active_list: ", active_list)

        chemical_list = self.chem_data.get_all_labels()
        chemical_delta_list = self._delta_names()

        conc_arr_before = self.get_diagnostic_conc_data().loc[row_baseline][chemical_list].to_numpy(dtype='float16')
        print("baseline concentrations: ", conc_arr_before)

        delta_cumulative = np.zeros(self.chem_data.number_of_chemicals(),
                                    dtype=float)  # One element per chemical species

        # For each reaction
        for rxn_index in range(self.chem_data.number_of_reactions()):
            if (rxn_index in active_list):  # If reaction is tagged as "fast"
                row = row_list[rxn_index]   # Row in the data frame for the diagnostic data on this reaction
                delta_rxn = self.get_rxn_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy(dtype='float16')
                print(f"From fast rxn {rxn_index}: delta_rxn = {delta_rxn}")
            else:                           # If reaction is tagged as "slow"
                row = row_list[rxn_index]   # Row in the data frame for the diagnostic data on this reaction
                delta_rxn = self.get_rxn_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy(dtype='float16')
                #delta_rxn = np.zeros(self.reaction_data.number_of_chemicals(),
                #                     dtype=float)   # This is the way it was in release beta 19
                print(f"From slow rxn {rxn_index}: delta_rxn = {delta_rxn}")


            delta_cumulative += delta_rxn
        # END for

        print("delta_cumulative: ", delta_cumulative)

        conc_after = conc_arr_before + delta_cumulative
        print("updated concentrations: ", conc_after)

        next_system_state = self.get_diagnostic_conc_data().loc[row_baseline + 1][chemical_list].to_numpy(dtype='float16')
        print(f"concentrations from the next row ({row_baseline + 1}) of the system state: ", next_system_state)

        status = np.allclose(conc_after.astype(float), next_system_state.astype(float))
        if status:
            print("Match OK")
        else:
            print("****************************************   MISMATCH!!!  ****************************************")

        print("-----------")

        return status



    def explain_time_advance(self, return_times=False, silent=False, sys_history=None) -> Union[None, tuple]:
        """
        Use the saved-up diagnostic data, to print out details of the timescales of the reaction run

        If diagnostics weren't enabled ahead of calling this function, an Exception is raised

        EXAMPLE of output:
            From time 0 to 0.0304, in 17 FULL steps of 0.0008
            (for a grand total of 38 FULL steps)

        :param return_times:    If True, all the critical times (times where the interval steps change)
                                    are saved and returned as a list
        :param silent:          If True, nothing gets printed out
        :param sys_history:     [OPTIONAL] The system history.
                                    If passed, use it in lieu of the diagnostic data;
                                    a consideration is the fact that the user might only have asked
                                    for PART of the history to be saved
        :return:                Depending on the argument "return_times", it returns either None, or a pair with 2 lists:
                                        1 - list of time values
                                        2 - list of step sizes  (will have one less element than the first list)
        """
        if sys_history:
            df = sys_history
            assert not df.empty , \
                "Diagnostics.explain_time_advance(): no history data found.  " \
                "Did you run the reaction simulation prior to calling this function?"
        else:
            df = self.get_diagnostic_conc_data()
            assert not df.empty , \
                "Diagnostics.explain_time_advance(): no diagnostic data found.  " \
                "Did you run the reaction simulation prior to calling this function?"

        n_entries = len(df)
        if sys_history:
            t = list(df["SYSTEM TIME"])     # List of times (simulation step points)
        else:
            t = list(df["TIME"])            # List of times (simulation step points)

        if n_entries == 1:
            print(f"Diagnostics.explain_time_advance(): no time advance found. "
                  f"Diagnostics only contain an initial System Time of {t[0]:.3g} ;  "
                  f"did you run the reaction simulation prior to calling this function?")
            return


        # If we get thus far, we have at least 2 entries in the list "t" of times
        grad = np.diff(t)
        grad_shifted = np.insert(grad, 0, 0)  # Insert a zero at the top (shifting down the array)
        #print(grad)
        #print(grad_shifted)
        start_interval = t[0]

        critical_times = [start_interval]  # List of times to report on (start time, plus times where the steps change)
        step_sizes = []

        total_number_full_steps = 0
        for i in range(1, n_entries-1):
            #print(i)
            if not np.allclose(grad[i], grad_shifted[i]):
                #print (f"   Detection at element {i} : time={t[i]}")
                #print (f"   From time {start_interval} to {t[i]}, in steps of {grad_shifted[i]}")
                n_full_steps_taken = self._explain_time_advance_helper(t_start=start_interval, t_end=t[i],
                                                                       delta_baseline=grad_shifted[i], silent=silent)

                start_interval = t[i]
                if n_full_steps_taken > 0:
                    total_number_full_steps += n_full_steps_taken
                    critical_times.append(t[i])
                    step_sizes.append(grad_shifted[i])


        # Final wrap-up of the interval's endpoint
        n_full_steps_taken = self._explain_time_advance_helper(t_start=start_interval, t_end=t[-1],
                                                               delta_baseline=grad_shifted[-1], silent=silent)

        if n_full_steps_taken > 0:
            total_number_full_steps += n_full_steps_taken
            critical_times.append(t[-1])
            step_sizes.append(grad_shifted[-1])

        if not silent:
            print(f"({total_number_full_steps} steps total)")

        if return_times:
            assert len(critical_times) == len(step_sizes) + 1, \
                "explain_time_advance(): validation error in the values to return"
            return (critical_times, step_sizes)



    def _explain_time_advance_helper(self, t_start, t_end, delta_baseline, silent: bool) -> Union[int, float]:
        """
        Using the provided data, about a group of same-size steps, create and print a description of it for the user

        :param t_start:
        :param t_end:
        :param delta_baseline:
        :param silent:          If True, nothing gets printed; otherwise, a line is printed out
        :return:                The corresponding number of FULL steps taken
        """
        if np.allclose(t_start, t_end):
            #print(f"   [Ignoring interval starting and ending at same time {t_start:.3g}]")
            return 0

        n_steps = round((t_end - t_start) / delta_baseline)
        step_s = "step" if n_steps == 1 else "steps"  # singular vs. plural
        if not silent:
            print(f"From time {t_start:.4g} to {t_end:.4g}, in {n_steps} {step_s} of {delta_baseline:.3g}")
        return n_steps



    def stoichiometry_checker(self, rxn_index :int, conc_arr_before: np.array, conc_arr_after: np.array) -> bool:
        """
        For the indicated reaction, investigate the change in the concentration of the involved chemicals,
        to ascertain whether the change is consistent with the reaction's stoichiometry.
        See https://life123.science/reactions

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        NOTE: the concentration changes in chemicals not involved in the specified reaction are ignored

        :param rxn_index:       The integer index (0-based) to identify the reaction of interest
        :param conc_arr_before: Numpy array with the concentrations of ALL the chemicals (whether involved
                                    in the reaction or not), in their index order, BEFORE the reaction step
        :param conc_arr_after:  Same as above, but after the reaction
                                TODO: maybe also accept a Panda's dataframe row

        :return:                True if the change in reactant/product concentrations is consistent with the
                                    reaction's stoichiometry, or False otherwise
        """
        assert self.reactions.number_of_reactions() == 1, \
            f"stoichiometry_checker() currently only works for 1-reaction simulations, but " \
            f"{self.chem_data.number_of_reactions()} reactions are present."

        assert len(conc_arr_before) == self.chem_data.number_of_chemicals(), \
            f"stoichiometry_checker() : the argument `conc_arr_before` must be a Numpy array " \
            f"with as many elements as the number of registered chemicals ({self.chem_data.number_of_chemicals()}); " \
            f"instead, it has {len(conc_arr_before)}"

        assert len(conc_arr_after) == self.chem_data.number_of_chemicals(), \
            f"stoichiometry_checker() : the argument `conc_arr_after` must be a Numpy array " \
            f"with as many elements as the number of registered chemicals ({self.chem_data.number_of_chemicals()}); " \
            f"instead, it has {len(conc_arr_after)}"

        self.reactions.assert_valid_rxn_index(rxn_index)

        return self._stoichiometry_checker_from_deltas(rxn_index = rxn_index,
                                                       delta_arr = conc_arr_after-conc_arr_before)



    def _stoichiometry_checker_from_deltas(self, rxn_index :int, delta_arr: np.array) -> bool:
        """
        Helper function.

        For the indicated reaction, investigate the change in the concentration of the involved chemicals,
        to ascertain whether the change is consistent with the reaction's stoichiometry.
        See https://life123.science/reactions

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        NOTE: the concentration changes in chemicals not involved in the specified reaction are ignored

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :param delta_arr:   Numpy array of numbers, with the concentrations changes of ALL the chemicals (whether involved
                                in the reaction or not), in their index order,
                                as a result of JUST the reaction of interest
        :return:            True if the change in reactant/product concentrations is consistent with the
                                reaction's stoichiometry, or False otherwise
                                Note: if any of the elements of the passed Numpy array is NaN, then True is returned
                                      (because NaN values are indicative of aborted steps; can't invalidate the stoichiometry
                                      check because of that)
        """
        #print("delta_arr: ", delta_arr)

        if np.isnan(delta_arr).any():
            return True         # The presence of a NaN, anywhere in delta_arr, is indicative of an aborted step

        rxn = self.reactions.get_reaction(rxn_index)
        reactants = self.reactions.get_reactants(rxn_index)
        products = self.reactions.get_products(rxn_index)

        # Pick (arbitrarily) the first reactant,
        # to establish a baseline change in concentration relative to its stoichiometric coefficient
        baseline_term = reactants[0]
        baseline_species_name = rxn.extract_species_name(baseline_term)
        baseline_species_index = self.chem_data.get_index(baseline_species_name)
        baseline_stoichiometry = rxn.extract_stoichiometry(baseline_term)
        baseline_ratio =  (delta_arr[baseline_species_index]) / baseline_stoichiometry
        #print("\nbaseline_ratio: ", baseline_ratio)

        for i, term in enumerate(reactants):
            if i != 0:
                species_name = rxn.extract_species_name(term)
                species_index = self.chem_data.get_index(species_name)
                stoichiometry = rxn.extract_stoichiometry(term)
                ratio =  (delta_arr[species_index]) / stoichiometry
                #print(f"ratio for `{self.chem_data.get_name(species)}`: {ratio}")
                if not np.allclose(ratio, baseline_ratio):
                    return False

        for term in products:
            species_name = rxn.extract_species_name(term)
            species_index = self.chem_data.get_index(species_name)
            stoichiometry = rxn.extract_stoichiometry(term)
            ratio =  - (delta_arr[species_index]) / stoichiometry     # The minus in front is b/c we're on the other side of the eqn
            #print(f"ratio for `{self.chem_data.get_name(species)}`: {ratio}")
            if not np.allclose(ratio, baseline_ratio):
                return False

        return True



    def stoichiometry_checker_entire_run(self) -> bool:
        """
        Verify that the stoichiometry is satisfied in all the reaction steps,
        using the diagnostic data from an earlier run

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        :return:    True if everything checks out, or False otherwise
        """
        # TODO: maybe make it also accept as inputs a reaction object and a Pandas dataframe
        if self.diagnostic_rxn_data == {}:
            print("WARNING *** In order to run stoichiometry_checker_entire_run(), "
                  "the diagnostics must be turned enabled prior to the simulation run!")
            return False

        number_rxns = self.chem_data.number_of_reactions()
        assert number_rxns == 1, \
            f"stoichiometry_checker_entire_run(): this function is currently designed for just 1 reaction " \
            f"(whereas {number_rxns} are present)"

        for rxn_index in range(number_rxns):
            diagnostic_data = self.get_rxn_data(rxn_index=rxn_index, print_reaction=False)
            for row_index in range(len(diagnostic_data)):
                df_row = diagnostic_data.loc[row_index]     # Row in the Panda's data frame of diagnostic data
                chemical_delta_list = self._delta_names()           # EXAMPLE: ["Delta A", "Delta B", "Delta C"]
                delta = df_row[chemical_delta_list].to_numpy(dtype='float32')      # Extract select columns from the data frame row, and turn into Numpy array
                status = self._stoichiometry_checker_from_deltas(rxn_index=rxn_index, delta_arr=delta)
                if not status:
                    print(f"Stoichiometry NOT satisfied by reaction # {rxn_index}: "
                          f"see row # {row_index} in the diagnostic data for that reaction, from get_diagnostic_rxn_data()")
                    print(df_row )
                    return False

        return True



    def _delta_names(self) -> [str]:
        """
        Return a list of strings, with the names of ALL the registered chemicals,
        in their index order, each prefixed by the string "Delta "
        EXAMPLE: ["Delta A", "Delta B", "Delta X"]

        :return:    A list of strings
        """
        chemical_list = self.chem_data.get_all_labels()      # EXAMPLE: ["A", "B", "X"]
        chemical_delta_list = ["Delta " + name
                               for name in chemical_list]
        return chemical_delta_list



    def _delta_conc_dict(self, delta_conc_arr :np.ndarray) -> dict:
        """
        Convert a Numpy array into a dict, based on all the registered chemicals.
        The keys are the chemical names, prefixed by "Delta "

        :param delta_conc_arr:  A Numpy array of "delta concentrations".  EXAMPLE: array[1.23, 52.2]
        :return:                A dictionary such as {"Delta A": 1.23, "Delta X": 52.2}
        """
        chemical_delta_list = self._delta_names()    # EXAMPLE: ["Delta A", "Delta X"]
        assert len(chemical_delta_list) == len(delta_conc_arr), \
            f"_delta_conc_dict(): mismatch in number of chemicals ({len(chemical_delta_list)} " \
            f"vs. {len(delta_conc_arr)})"

        delta_conc_dict = {}
        for i, delta_name in enumerate(chemical_delta_list):
            delta_conc_dict[delta_name] = delta_conc_arr[i]

        return delta_conc_dict
                            