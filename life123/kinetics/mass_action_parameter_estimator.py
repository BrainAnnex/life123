import numpy as np
from life123.numerical import Numerical
import plotly.graph_objects as pgo
from life123.visualization.plotly_helper import PlotlyHelper



class MassAction_Parameter_Estimator:
    """
    Estimate rate constants for concentration trajectory data
    """

    @staticmethod
    def estimate_rate_constants_simple(t :np.ndarray,
                                       A_conc :np.ndarray, B_conc :np.ndarray,
                                       reactant_name="Reactant", product_name="Product"):
        """
        Estimate the rate constants for a 1-st order reaction of the type A <-> B,
        given time evolution of [A] and [B] on a grid of time points (the points don't need to be equally spaced),
        and create a plot to show the fit

        IMPORTANT : This is for reactions with a 1:1 stoichiometry between the given reactant and product

        :param t:               A numpy array of time grid points at which the other functions are specified
                                    (the points do not need to be equally-spaced)
        :param A_conc:          A numpy array of the concentrations of the reactant, at the times in the array t
        :param B_conc:          A numpy array of the concentrations of the product, at the times in the array t
        :param reactant_name:   [OPTIONAL] The name of the reactant (for display purposes)
        :param product_name:    [OPTIONAL] The name of the product (for display purposes)
        :return:                A plotly "Figure" object.  The estimated rate constants are printed out
        """
        #TODO: pytest

        total_conc_arr = A_conc + B_conc
        total_conc = np.median(total_conc_arr)    # TODO: give warning or abort if there's too much variance
        sd = np.std(total_conc_arr)

        print(f"Reaction {reactant_name} <-> {product_name}")
        print(f"Total REACTANT + PRODUCT has a median of {total_conc:,.4g}, "
              f"\n    with standard deviation {sd:,.4g} (ideally should be zero)")


        # The rate of change of reactant concentration with time
        A_prime = np.gradient(A_conc, t, edge_order=2)
        # The rate of change of product concentration with time
        B_prime = np.gradient(B_conc, t, edge_order=2)

        median_sum_derivs = np.median(A_prime + B_prime)
        print(f"The sum of the time derivatives of the reactant and the product "
              f"\n    has a median of {median_sum_derivs:,.4g} (ideally should be zero)")


        # Do a least-square fit
        kF, kR = Numerical.two_vector_least_square(V = A_conc, W = -B_conc, Y = B_prime)

        print(f"Least square fit to model as elementary reaction: {product_name}'(t) = kF * {reactant_name}(t) - kR * {product_name}(t)")

        # Plot both Y and its least-square fit, as functions of t
        fig_main = PlotlyHelper.plot_curves(x=t, y=[B_prime , kF * A_conc - kR * B_conc],
                                       title=f"d/dt {product_name}(t) and its least-square fit",
                                       x_label="t", y_label=f"d/dt  {product_name}(t)",
                                       legend_title="Curves",
                                       curve_labels=[f"d/dt {product_name}(t): exact", f"d/dt {product_name}(t): least-square fit"],
                                       colors=['green', 'red'])

        '''
        # Plot both Y and its least-square fit, as functions of X
        fig = PlotlyHelper.plot_curves(x=A_conc, y=[B_prime , kF * A_conc - kR * B_conc],
                                       title=f"d/dt {product_name}(t) as a function of {reactant_name}(t), alongside its least-square fit",
                                       x_label=f"{reactant_name}(t)", y_label=f"{product_name}'(t)",
                                       curve_labels=[f"{product_name}'(t)", "Linear Fit"], legend_title="Curve vs Fit:",
                                       colors=['green', 'red'])
        '''

        fig_side = PlotlyHelper.plot_curves(x=A_conc, y=B_prime,
                                         title=f"d/dt {product_name}(t) as a function of {reactant_name}(t)",
                                         x_label=f"{reactant_name}(t)", y_label=f"d/dt  {product_name}(t)",
                                         colors="purple")

        print(f"\n-> ESTIMATED RATE CONSTANTS: kF = {kF:,.4g} , kR = {kR:,.4g}")

        return PlotlyHelper.combine_in_vertical_grid(fig1=fig_main, fig2=fig_side,
                                                     title1=f"d/dt {product_name}(t) and its least-square fit",
                                                     title2=f"d/dt {product_name}(t) as a function of {reactant_name}(t)",
                                                     title_combined="LEAST-SQUARE FIT ANALYSIS:")



    @staticmethod
    def estimate_rate_constants_synthesis(t :np.ndarray,
                                          A_conc :np.ndarray, B_conc :np.ndarray, C_conc :np.ndarray,
                                          reactants :[str, str], product :str) -> pgo.Figure:
        """
        Estimate the rate constants for a 1-st order association (synthesis) reaction of the type A + B <-> C,
        given time evolution of [A], [B] and [C] on a grid of time points (don't need to be equally spaced)

        IMPORTANT : This is for reactions with a 1:1:1 stoichiometry

        :param t:           A numpy array of time grid points where the other functions are specified
        :param A_conc:      A numpy array of the concentrations of the reactant, at the times in the array t
        :param B_conc:
        :param C_conc:
        :param reactants:   [OPTIONAL] A list with the names of the 2 reactants, in order (for display purposes)
        :param product:     [OPTIONAL] The name of the product (for display purposes)
        :return:            A plotly "Figure" object.  The estimated rate constants are printed out
        """
        #TODO: pytest

        # The rate of change of [product] with time
        Deriv_C = np.gradient(C_conc, t, edge_order=2)

        # Do a least-square fit
        kF, kR = Numerical.two_vector_least_square(V = A_conc * B_conc, W = - C_conc, Y = Deriv_C)

        print(f"Least square fit to {product}'(t) = kF * {reactants[0]}(t) * {reactants[1]}(t) + kR * (- {product}(t) )")

        # Plot both Y and its least-square fit, as functions of X
        fig = PlotlyHelper.plot_curves(x=A_conc, y=[Deriv_C , kF * A_conc * B_conc - kR * C_conc],
                                       title=f"d/dt {product}(t) as a function of {reactants[0]}(t), alongside its least-square fit",
                                       x_label=f"{reactants[0]}(t)", y_label=f"{product}'(t)",
                                       curve_labels=[f"{product}'(t)", "Linear Fit"], legend_title="Curve vs Fit:",
                                       colors=['green', 'red'])

        print(f"\n-> ESTIMATED RATE CONSTANTS: kF = {kF:,.4g} , kR = {kR:,.4g}")

        return fig
