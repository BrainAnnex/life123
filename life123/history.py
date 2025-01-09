import os
import csv
import pandas as pd
from life123.collections import CollectionTabular, Collection



class History:

    def __init__(self, active=False, chem_labels=None, frequency=1):      # obj
        """
        For the management of historical data

        :param chem_labels:
        :param active:
        :param frequency:
        """
        self.history = None                 # Object set by the children classes
        self.active = active
        self.capture_frequency = frequency

        self.restrict_chemicals = chem_labels

        self.step_of_last_snapshot = -1     # An out-of-range value

        self.log_file = None



    def get_history(self):
        return self.history



    def enable_history_common(self, frequency=1, chem_labels=None) -> None:
        """

        :param frequency:
        :param chem_labels:
        :return:            None
        """
        self.active = True
        self.capture_frequency = frequency
        self.restrict_chemicals = chem_labels



    def disable_history(self):
        self.active = False



    def to_capture(self, step_count=None) -> bool:
        """
        Return True if various criteria indicate that a snapshot should be taken at the give step_count;
        False otherwise

        :param step_count:  None means always recapture
        :return:
        """
        if step_count is None:
            return True

        if not self.active:
            return False

        if step_count == self.step_of_last_snapshot:
            return False     # Don't re-take snapshot, if still at the same step count

        if self.capture_frequency == 1:
            return True

        if (step_count+1) % self.capture_frequency == 0:
            return True

        return False



    def start_csv_log(self, log_file :str) -> None:
        """
        Register the specified filename for all future CSV logs.
        Any existing file by that name will get over-written.

        :param log_file:Name of a file to use for the CSV logs
        :return:        None
        """
        # TODO: validate log_file
        self.log_file = log_file

        if os.path.exists(log_file):
            os.remove(log_file)
            print(f"-> CSV-format output will be LOGGED into the file '{log_file}' . An existing file by that name was over-written")
        else:
            print(f"-> CSV-format output will be LOGGED into the file '{log_file}'")



    def set_caption_last_snapshot(self, caption) -> None:
        """
        Add a caption to the very last entry in the system history.
        No action taken is caption is None

        :param caption:
        :return:
        """
        if caption is not None:
            self.history.set_caption_last_snapshot(caption)




##########################################################################################################################

class HistoryUniformConcentration(History):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args
        self.history = CollectionTabular(parameter_name="SYSTEM TIME")  # To store user-selected snapshots of quantities of interest



    def enable_history(self, take_snapshot=True, frequency=1, chem_labels=None):
        """

        :param frequency:
        :param chem_labels:
        :param take_snapshot:   If True, a snapshot of the system's current configuration is added to the history
        :return:
        """
        self.enable_history_common(frequency=frequency, chem_labels=chem_labels)

        #if take_snapshot:
            #self.capture_snapshot()



    def save_snapshot(self, system_time, data_snapshot, step_count=None, caption="", initial_caption="") -> None:
        """

        :param system_time:
        :param data_snapshot:   EXAMPLE: {"A": 1.3, "B": 4.9}
        :param step_count:
        :param caption:         [OPTIONAL] String to save alongside this snapshot
        :return:                None
        """
        if not self.to_capture(step_count):
            return

        if (not caption) and initial_caption and (step_count == 0):
            caption = initial_caption

        self.history.store(par=system_time,
                           data_snapshot=data_snapshot, caption=caption)

        if step_count is not None:
            self.step_of_last_snapshot = step_count


        if self.log_file is None:
            return

        # If we get thus far, we have a log file to manage
        if os.path.exists(self.log_file):
            new_file = False
        else:
            new_file = True

        # Open a CSV file in write mode
        with open(self.log_file, mode="a", newline="") as fh:
            fieldnames = ["SYSTEM TIME"] + list(data_snapshot.keys()) + ["caption"]
            writer = csv.DictWriter(fh, fieldnames = fieldnames)

            data_snapshot["SYSTEM TIME"] = system_time
            data_snapshot["caption"] = caption

            # Write the header (only of a newly-created file)
            if new_file:
                writer.writeheader()

            # Write the dictionary as a row
            writer.writerow(data_snapshot)





##########################################################################################################################

class HistoryBinConcentration(History):

    def __init__(self, bins=None, *args, **kwargs):
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args
        self.restrict_bins = bins
        self.history = Collection(parameter_name="SYSTEM TIME")  # To store user-selected snapshots of quantities of interest



    def enable_history(self, frequency=1, chem_labels=None, bins=None):
        """

        :param frequency:
        :param chem_labels:
        :param bins:
        :return:
        """
        self.enable_history_common(frequency=frequency, chem_labels=chem_labels)
        self.restrict_bins = bins



    def save_snapshot(self, system_time, data_snapshot, step_count=None) -> None:
        """

           EXAMPLE of data_snapshot:
                { (0,0): {"A": 1.3, "B": 3.9},
                  (3,5): {"A": 4.6, "B": 2.7}
                }
        :param system_time:
        :param data_snapshot:
        :param step_count:
        :return:            None
        """
        if not self.to_capture(step_count):
            return

        self.history.store(par=system_time,
                           data_snapshot=data_snapshot)

        if step_count is not None:
            self.step_of_last_snapshot = step_count



    def bin_history(self, bin_address, include_captions=False):
        """
        Return the history at the given bin, as a Pandas dataframe.
        The first column is "SYSTEM TIME" and the other ones are the various chemicals for which
        history had been enabled.

        :param bin_address:         EXAMPLE, in 2D : (3,3)
        :param include_captions:    If True, the captions are returned as an extra "caption" column at the end
        :return:                    A Pandas data frame
        """
        first_snapshot = self.history.get_data()[0]         # EXAMPLE: {(3, 3): {'A': 0.0, 'B': 0.0, 'C': 0.0}}
        first_conc_data = first_snapshot.get(bin_address)   # EXAMPLE: {'A': 0.0, 'B': 0.0, 'C': 0.0}
        assert first_conc_data is not None, \
            f"bin_history(): no concentration historical data available for bin {bin_address}"

        chem_labels = list(first_conc_data)                 # Turn the dict's keys into list.  EXAMPLE: ['A', 'B', 'C']

        if include_captions:
            cols = ['SYSTEM TIME'] + chem_labels + ['caption']
        else:
            cols = ['SYSTEM TIME'] + chem_labels

        df = pd.DataFrame(columns = cols) # EXAMPLE: ['SYSTEM TIME', 'A', 'B', 'C']

        for row in self.history.get_collection():
            t , data, caption = row           # Unpack
            conc_dict = data[bin_address]
            l = len(df)
            df.loc[l] = conc_dict           # Add a new row, with all columns except "SYSTEM TIME"
            df.loc[l, "SYSTEM TIME"] = t    # Set the "SYSTEM TIME" cell for the row just added
            if include_captions:
                df.loc[l, "caption"] = caption

        return df




##########################################################################################################################

class HistoryReactionRate(History):

    def __init__(self, rxns=None, *args, **kwargs):
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args
        self.restrict_rxns = rxns
        self.history = CollectionTabular(parameter_name="SYSTEM TIME")  # To store user-selected snapshots of quantities of interest



    def enable_history(self, take_snapshot=True, frequency=1, chem_labels=None, rxns=None):
        """

        :param frequency:
        :param chem_labels:
        :param rxns:
        :param take_snapshot:   If True, a snapshot of the system's current configuration is added to the history
        :return:
        """
        self.enable_history_common(frequency=frequency, chem_labels=chem_labels)
        self.restrict_rxns = rxns
        if take_snapshot:
            self.capture_snapshot()



    def capture_snapshot(self, step_count=None, capture_initial_caption=None) -> None:
        """

        :param step_count:
        :param capture_initial_caption: [OPTIONAL]
        :return:                        None
        """
        if not self.to_capture(step_count):
            return

        system_time = self.obj.get_system_time()

        if (step_count == 0) or (step_count == self.step_of_last_snapshot + 1):
            data_snapshot = {}     # Reaction-rates snapshot
            for k, v in self.obj.system_rxn_rates.items():
                data_snapshot[f"rxn{k}_rate"] = v      # EXAMPLE:  "rxn4_rate" = 18.2

            self.history.store(par=system_time, data_snapshot=data_snapshot)

        self.step_of_last_snapshot = step_count



    def add_caption_to_last_snapshot(self, caption) -> None:
        """
        Add a caption to the very last entry in the system history.
        No action take if the caption is blank or None

        :param caption:
        :return:        None
        """
        if caption:
            self.history.set_caption_last_snapshot(caption)
