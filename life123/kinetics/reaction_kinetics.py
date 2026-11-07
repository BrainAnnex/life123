import math
import numpy as np
import plotly.graph_objects as pgo
from life123.numerical import Numerical
from life123.visualization.plotly_helper import PlotlyHelper




class ReactionKinetics:
    """
    TODO: absorb into new modules

    Static methods about reactions kinetics

    For background, see https://life123.science/reactions
    """


    #####################################################################################################

    '''                                 ~   HALF-TIME RELAXATION   ~                                  '''

    def ________HALF_TIME_RELAXATION________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    @staticmethod
    def half_time_unimolecular_irreversible(kF) -> float:
        """
        Return the time taken for the reactant concentration in an irreversible unimolecular reaction
        to decrease by half (halfway to the asymptotic state, which happens to be zero)

        :param kF:  Forward reaction rate constant
        :return:    Reaction time taken for reactant concentration to decrease by half
        """
        return 0.6931471805599 / kF     # ln 2 / kF



    @staticmethod
    def half_time_relaxation_unimolecular_reversible(kF, kR) -> float:
        """
        Return the time taken for the reactant concentration in a reversible unimolecular reaction R <-> P
        to decrease halfway to its asymptotic equilibrium state

        :param kF:  Forward reaction rate constant
        :param kR:  Reverse reaction rate constant
        :return:    Reaction time taken for reactant concentration
                        to decrease halfway to its asymptotic equilibrium state
        """
        # TODO: maybe simply combine with half_time_unimolecular_irreversible, with kR=0
        return 0.6931471805599 / (kF + kR)      # ln 2 / (kF + kR)


    @staticmethod
    def relaxation_time_unimolecular_reversible(kF, kR) -> float:
        """
        This is the same as the half-time of relaxation, within a factor of ln 2

        :param kF:  Forward reaction rate constant
        :param kR:  Reverse reaction rate constant
        :return:    A quantity named "relaxation time"
        """
        # TODO: is this method needed?
        return 1. / (kF + kR)



    @staticmethod
    def half_time_to_equilibrium_irreversible_synthesis(kF, A0, B0) -> float:
        """
        Return the time taken for the reactant concentration
        in a synthesis reaction A + B -> P
        to decrease halfway to their asymptotic equilibrium state.
        The time taken for 50% of the limiting reactant to be consumed

        :param kF:  Forward reaction rate constant
        :return:    Reaction time taken for reactant concentration
                        to decrease halfway to its asymptotic equilibrium state
        """
        #TODO: test
        #equil_concs = ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, b=1, B0=B0, p=1, P0=P0)

        delta_conc = A0 - B0
        if np.allclose(delta_conc, 0):
            return 1 / (kF * A0)

        if delta_conc > 0:      # When A0 > B0
            return math.log((2 * A0 - B0) / A0) / (kF * delta_conc)
        else:
            return - math.log((2 * B0 - A0) / B0) / (kF * delta_conc)   #TODO: test
