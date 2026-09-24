import math
import cmath
import numpy as np
import plotly.graph_objects as pgo
from life123.numerical import Numerical
from life123.visualization.plotly_helper import PlotlyHelper
from typing import Tuple



class ReactionKinetics:
    """
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
        Return the time taken for the reactant concentration in a reversible unimolecular reaction
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




    #####################################################################################################

    '''                                      ~   RATES   ~                                            '''

    def ________RATES________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

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



    @staticmethod
    def compute_rate_elementary(reactants :list[str], products :list[str],
                                kF :float, kR :float, reversible :bool,
                                conc_dict :dict) -> float:
        """
        Given a SINGLE elementary reaction, in 1st order to all its chemical species,
        and the specified concentrations of chemicals,
        compute its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

        CAUTION: even though arbitrary lists of reactants and products are accepted as arguments,
                 this kinetic model will generally only hold for elementary reactions

        :param reactants:   List of the species id's of the reactants
        :param products:    List of the species id's of the products
        :param kF:          Forward reaction rate
        :param kR:          Reverse reaction rate; ignored if irreversible
        :param reversible:  True if the reaction is reversible; False otherwise
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}

        :return:            The differences between the reaction's forward and reverse rates
        """
        """
        # TODO: this warning doesn't belong here
        issue_warning = False
        if len(reactants) > 2 or len(products) > 2:
            issue_warning = True

        if len(reactants) == 2 and len(products) > 1:
            issue_warning = True

        if len(products) == 2 and len(reactants) > 1:
            issue_warning = True

        if issue_warning:
            print("compute_rate_elementary(): WARNING - using 1st order kinetic modeling "
                  "for a reaction that's probably not elementary")
        """

        forward_rate = kF        # The initial multiplicative factor
        for r in reactants:
            # Process all the reactants
            conc = conc_dict[r]
            forward_rate *= conc

        if not reversible:
            return forward_rate

        reverse_rate = kR        # The initial multiplicative factor
        for p in products:
            # Process all the reaction products
            conc = conc_dict[p]
            reverse_rate *= conc

        return forward_rate - reverse_rate



    @staticmethod
    def compute_rate_mass_action_kinetics(reactant_terms :[(int, str)], product_terms :[(int, str)],
                                       kF :float, kR :float,
                                       conc_dict :dict) -> float:
        """
        Given a SINGLE arbitrary complex reaction,
        and the specified concentrations of its chemicals,
        compute its initial reaction "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step,
        when the reaction kinetics can be modeled as if it were an elementary reaction.

        This function is largely a convenient default function for testing,
        in scenarios (and hypothetical scenarios) where the the reaction's kinetics
        follow the familiar "Rate Laws",
        with the order of the reaction with respect to its reactant and products
        equals their respective stoichiometric coefficients -
        in other words an elementary reaction or one (hypothetically) modeled as such.

        For example, if the reaction is aA + bB <-> pP + qQ,
        then this function returns:  kF [P]^p [Q]^q - kR [A]^a [B]^b

        Warning: generally speaking, this is NOT a valid kinetic modeling
        of any reaction that isn't elementary

        :param reactant_terms:  A list of pairs (stoichiometry coefficient , species id) for the reactants
        :param product_terms:   A list of pairs (stoichiometry coefficient , species id) for the products
        :param kF:              Forward reaction rate
        :param kR:              Reverse reaction rate; zero if the reaction is irreversible

        :param conc_dict:       A dict mapping chemical labels to their concentrations,
                                    for all the chemicals involved in this reaction
                                    EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}

        :return:            The differences between the reaction's forward and reverse rates
        """
        forward_rate = kF        # The initial multiplicative factor
        for order, reactant_label in reactant_terms:     # The stoichiometry coeff. of each reactant is taken to be its reaction order
            conc = conc_dict.get(reactant_label)
            assert conc is not None, \
                f"compute_rate_mass_action_kinetics(): missing concentration value for reactant chemical with label `{reactant_label}`"
            forward_rate *= conc ** order       # Raise to power


        if kR == 0:
            return forward_rate                 # If there's no reverse reaction (i.e., if reaction is irreversible)


        reverse_rate = kR        # The initial multiplicative factor
        for order, product_label in product_terms:      # The stoichiometry coeff. of each product is taken to be its reaction order
            conc = conc_dict.get(product_label)
            assert conc is not None, \
                f"compute_rate_mass_action_kinetics(): missing concentration value for product chemical `{product_label}`"
            reverse_rate *= conc ** order       # Raise to power

        return forward_rate - reverse_rate



    @staticmethod
    def kinetic_rate_first_order(stoichiometry,
                                 kinetic_parameters :dict,
                                 conc_dict :dict) -> float:
        """
        If the reactions isn't elementary, this is a HYPOTHETICAL scenario (mostly for testing and analysis)
        where the reaction is first order in EACH of the reactants and EACH of products
        """
        kF = kinetic_parameters.get("kF")
        kR = kinetic_parameters.get("kR")
        reversible = True if kR else False

        # Pretend that the reaction is an elementary one
        reactants = stoichiometry.get_reactant_ids()
        products  = stoichiometry.get_product_ids()
        return ReactionKinetics.compute_rate_elementary(reactants=reactants, products=products,
                                                        kF=kF, kR=kR, reversible=reversible,
                                                        conc_dict=conc_dict)





    #####################################################################################################

    '''                           ~   EQUILIBRIUM CONCENTRATIONS   ~                                  '''

    def ________EQUILIBRIUM_CONCENTRATIONS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    @staticmethod
    def compute_equilibrium_conc_elementary_decomposition(kF, kR, A0, P0) -> dict:
        """
        Given a reversible reaction of the form:  A <-> 2P
        whose kinetics follow the "mass-action" law,
        determine the equilibrium concentrations that would be eventually reached by all its species,
        given the specified initial concentrations,
        IN THE ABSENCE of any other reaction.

        :param kF:  The reaction's forward rate constant
        :param kR:  The reaction's reverse rate constant
        :param A0:  Initial concentration of the reactant species `A`
        :param P0:  Initial concentration of the product species `P`
        :return:    A dictionary with two keys, `A` and `P`, containing their equilibrium concentrations.
                        EXAMPLE:    {'A': 24.0, 'P': 1.8}
        """
        # Reverse kF and kR, to obtain the reversed reaction 2 P <-> A
        result = ReactionKinetics.compute_equilibrium_conc_elementary_synthesis(kF=kR, kR=kF, A0=P0, P0=A0)

        return {"A": result["P"], "P": result["A"]}     # Reverse `A` and `P`, because the above function
                                                        # was in terms of 2 A <-> P
        """
        # ALTERNATIVE WAY (successfully tested)
        # Note: the confusing arguments "Q0=P0, p=2, q=2" are a hack to use "Generalization 2" of _compute_equilibrium_conc_first_order()
        # TODO: ditch the helper function, and do a direct computation
        result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR,
                                                                        A0=A0, a=1,
                                                                        P0=P0, Q0=P0, p=2, q=2)
        
        del result["Q"]
        """
        return result


    @staticmethod
    def compute_equilibrium_conc_elementary_synthesis(kF, kR, A0, P0) -> dict:
        """
        Given a reversible reaction of the form:  2 A <-> P
        whose kinetics follow the "mass-action" law,
        determine the equilibrium concentrations that would be eventually reached by all its species,
        given the specified initial concentrations,
        IN THE ABSENCE of any other reaction.

        :param kF:  The reaction's forward rate constant
        :param kR:  The reaction's reverse rate constant
        :param A0:  Initial concentration of the reactant species `A`
        :param P0:  Initial concentration of the product species `P`
        :return:    A dictionary with two keys, `A` and `P`, containing their equilibrium concentrations.
                        EXAMPLE:    {'A': 24.0, 'P': 1.8}
        """
        # Note: the confusing arguments "B0=A0, a=2, b=2" are a hack to use "Generalization 2" of _compute_equilibrium_conc_first_order()
        # TODO: ditch the helper function, and do a direct computation
        result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, A0=A0, B0=A0, a=2, b=2, P0=P0, p=1)
        del result["B"]
        return result



    @staticmethod
    def compute_equilibrium_conc_mass_action(kF, kR, A0, P0, B0=None, Q0=None) -> dict:
        """
        Given a reversible reaction of the form:  A + B <-> P + Q (some terms may be missing)
        whose kinetics follow the "mass-action" law,
        determine the equilibrium concentrations that would be eventually reached by all its species,
        given the specified initial concentrations,
        IN THE ABSENCE of any other reaction.

        Note: for scenarios where one of stoichiometric coefficients is 2,
              use compute_equilibrium_conc_elementary_synthesis() instead

        :param kF:  The reaction's forward rate constant (None will be interpreted as a zero, i.e. no forward reaction)
        :param kR:  The reaction's reverse rate constant (None will be interpreted as a zero, i.e. no reverse reaction)
        :param A0:  The initial concentration of species `A` (i.e. the 1st reactant term)
        :param P0:  The initial concentration of species `P` (i.e. the 1st product term)
        :param B0:  Use None to indicate the absence of a `B` species (i.e. a 2nd term) among the reactants
        :param Q0:  Use None to indicate the absence of a `Q` species (i.e. a 2nd term) among the products

        :return:    A dictionary of the equilibrium concentrations of the
                        species involved in the specified reaction;
                        only the applicable entries will be present
                        EXAMPLES:
                            {'A': 24.0, 'B': 36.0, 'P': 1.8, 'Q': 0.3}
                            {'A': 35.0, 'P': 18.3}
        """
        b = 0 if B0 is None else 1
        q = 0 if Q0 is None else 1

        kF = 0 if kF is None else kF
        kR = 0 if kR is None else kR

        return ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=b, B0=B0, q=q, Q0=Q0)



    @staticmethod
    def _compute_equilibrium_conc_first_order(kF, kR, a, A0, p, P0, b=0, B0=None, q=0, Q0=None) -> dict:
        """
        Determine the equilibrium concentrations that would be eventually reached
        by the chemicals participating in a reversible reaction of the form:
            a A + b B <-> p P + q Q
        whose kinetics are (hypothetically) FIRST ORDER in EACH of the species
        (some of the terms may be missing),
        given the specified initial concentrations,
        IN THE ABSENCE of any other reaction.

        If both a and b are non-zero, or if both p and q are non-zero, the overall order will be 2;
        otherwise, it will be 1.

        In other words, this applies to reactions that can be modeled kinetically as:
            vF = kF [A] [B]  ,  vR = kR [C] [D]    (some of the terms may be missing)
        where vF and vR are, respectively, the forward and reverse reaction velocities (rates).

        IMPORTANT:  cannot be used for an elementary reaction such as 2 A <-> P,
                    because it won't be first-order with respect to `A`

        CAUTION: since all reaction orders are expected to be 1,
                 if any of stoichiometry coefficients of the involved species aren't 1,
                 then it's a scenario of a particular scenario of "generalized kinetics"
                 (non-mass-action with respect to the stoichiometry)

        Note:  we're using "concentrations" instead of "chemical activities";
               concentrations approximate the activities of ideal dilute solutions

        :param kF:  The reaction's forward rate constant.
                        If the reaction isn't elementary (but can be modeled as being of 1st order in each of the species),
                        then this value will be a composite "kF effective"
        :param kR:  The reaction's reverse rate constant.
                        If the reaction isn't elementary (but can be modeled as being of 1st order in each of the species),
                        then this value will be a composite "kR effective"

        :param a:   The stoichiometry coefficient of species A in the reaction (typically 1)
        :param A0:  The initial concentration of species A
        :param b:   [OPTIONAL] The stoichiometry coefficient of species B in the reaction;
                        if not present, accept the default value of zero
        :param B0:  [OPTIONAL] The initial concentration of species B;
                        if not present, accept the default None
        :param p:   The stoichiometry coefficient of product species `P` in the reaction (typically 1)
        :param P0:  The initial concentration of the product species `P`
        :param q:   [OPTIONAL] The stoichiometry coefficient of species D in the reaction;
                        if not present, accept the default value of zero
        :param Q0:  [OPTIONAL] The initial concentration of species `Q``;
                        if not present, accept the default value

        :return:    A dictionary of the equilibrium concentrations of the
                        species involved in the specified reaction;
                        only the applicable entries will be present
                        EXAMPLES:
                            {'A': 24.0, 'B': 36.0, 'P': 1.8, 'Q': 0.3}
                            {'A': 35.0, 'P': 18.3}
        """
        '''
        For (hypothetical) reactions of the form aA + bB <-> cC + dD  
        that are FIRST-ORDER in all species, and with the given initial condition,
        the equilibrium equation needs to equate the reaction's forward and reverse rates 
        after the reaction has advanced by m moles 
        (which consumes reagents in proportion to their stoichiometric coefficients, and 
        generates products similarly):
                  
            kF [(A0 - a*m) (B0 - b*m)] = kR [(C0 + c*m) (D0 + d*m)]     # Forward rate = Reverse rate
        
        where m (to be solved for) is the number of "moles/liter of forward reaction" 

        Generalization 1: 
        The above equilibrium equation can also be made to handle reaction terms that aren't present 
        (for example the "D" part in the reaction A + B <-> C), 
        by simply using 1 for the "initial concentration" and 0 for the "stoichiometry coefficient"
        of any missing term;
        such a trick will make its corresponding multiplicative term to become (1 + 0*m) ,
        i.e. identically equal to 1 , and thus have no effect on the solution of the equation.
        
        Generalization 2:
        The equilibrium equation can also be made to handle scenarios where 2 terms refer to the same chemical,
        such as in the reaction 2 A <-> C , with 2nd order with respect to A, 
        by setting B0 = A0, a=2, b=2,
        which will produce the multiplicative term that we need, namely:
            kF [(A0 - 2*m) **2]
        
        
        SOLUTION -
        Our equilibrium equation can be expanded into a standard quadratic form for the unknown m :
        
                alpha * m**2 + beta * m + gamma = 0
        
        (where alpha, beta and gamma are as computed below), and then solved for m.  
        In case of two solution for the quadratic, we'll pick the one that
        leads to physically-possible results (non-negative concentrations of all the species.)
        '''
        assert kR is not None, \
            "_compute_equilibrium_conc_first_order(): a value must be passed for argument `kR` (currently None)"

        assert kF is not None, \
            "_compute_equilibrium_conc_first_order(): a value must be passed for argument `kF` (currently None)"

        if b == 0:
            assert B0 is None, \
                "_compute_equilibrium_conc_first_order(): unexpected concentration value for a chemical not in reaction (B)"
            B0 = 1      # A trick to reduce to the general equation when b=0

        if q == 0:
            assert Q0 is None, \
                "_compute_equilibrium_conc_first_order(): unexpected concentration value for a chemical not in reaction (D)"
            Q0 = 1      # A trick to reduce to the general equation when d=0

        alpha = (kF * a * b - kR * p * q)
        beta = -kF * (b*A0 + a*B0) -kR * (q * P0 + p * Q0)
        gamma = kF * A0 * B0 - kR * P0 * Q0

        #print("_compute_equilibrium_conc_first_order() - alpha, beta, gamma : ", alpha, beta, gamma)

        if np.allclose(alpha, 0):
            # The quadratic reduces to the linear equation:  beta * m + gamma = 0
            m1 = -gamma / beta
            m2 = m1
        else:
            sqrt_discriminant = math.sqrt(beta**2 - 4 * alpha * gamma)
            m1 = (-beta - sqrt_discriminant) / (2 * alpha)
            m2 = (-beta + sqrt_discriminant) / (2 * alpha)
            #print("m1, m2 : ", m1, m2)

        m = m1  # Let's start with one of the 2 possible solutions of the quadratic equation

        # After m "moles of forward reaction", the concentration of the reactant "A"
        # in aA + bB <-> cC + dD gets reduced by a*m . Likewise for the other terms.
        # Reaction products get increased.  Values for missing terms will be meaningless
        std_result = {"A" : A0 - a*m, "B" : B0 - b*m, "P" : P0 + p * m, "Q" : Q0 + q * m}   # TENTATIVE values!

        if min(std_result.values()) < 0:    # If there's any negative value in the concentrations...
            # ...then repeat the computation using the other solution of the quadratic
            m = m2
            print("Using 2nd solution: m = ", m)
            std_result = {"A" : A0 - a*m, "B" : B0 - b*m, "C" : P0 + p * m, "D" : Q0 + q * m}

        # Eliminate any entries that aren't applicable to the reaction
        if b == 0:
            del std_result["B"]
        if q == 0:
            del std_result["Q"]

        return std_result



    @classmethod
    def compute_reaction_quotient(cls, reactant_data :[str|Tuple[int, str]], product_data :[str|Tuple[int, str]],
                                  conc :dict, explain=False) -> np.double | Tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" Q (aka "Mass–action Ratio"),
        for the reaction with the specified parameters,
        given the concentrations of species involved in the reaction.

        EXAMPLE: use reactant_data=[(2, "A"), (1, "B")] and product_data=["C"] ,
                 for a reaction of the form 2 A + B <-> C ,
                 alongside a dictionary with the concentrations (activities) of A, B and C;
                 the result will be given by the formula:  [C] / ( [A]^2 [B])    , where ^2 represents squaring

        Note: in a heterogeneous mixture, solids, pure liquids and solvents have an activity that has a fixed value of 1,
              and should be omitted from the parameters passed to this function.
              We're using the term "concentrations" instead of "chemical activities";
              concentrations approximate the activities of ideal dilute solutions

        :param reactant_data:   List whose elements can be either STRINGS with the labels of the reactants,
                                    or PAIRS of the form (stoichiometry coefficient, label) of the reactants.

        :param product_data:    List whose elements can be either STRINGS with the labels of the products of the reactions,
                                    or PAIRS of the form (stoichiometry coefficient, label) of the products.

        :param conc:            Dictionary with the concentrations (activities) of the species involved in the reaction.
                                The keys are the chemical labels
                                    EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param explain:         If True, it also returns the math formula being used for the computation
                                    EXAMPLES:   "([C][D]) / ([A][B])"
                                                "[B] /  [A]^2 "

        :return:                If explain is False, return a value for the "Reaction Quotient" (aka "Mass–action Ratio");
                                    if True, return a pair with that quotient and a string with the math formula that was used.
                                    Note that the reaction quotient is a Numpy scalar that might be np.inf or np.nan
        """
        # TODO: could be tidier in avoiding unnecessary blanks in the explanations
        numerator = np.double(1)    # The product of all the concentrations of the reaction products (adjusted for reaction order)
        denominator = np.double(1)  # The product of all the concentrations of the reactants (also adjusted for reaction order)

        numerator_text = ""      # First part of the textual explanation
        denominator_text = ""    # Second part of the textual explanation


        # Compute the numerator of the "Reaction Quotient"
        for term in product_data:
            # Loop over the reaction products
            if type(term) == str:
                stoich_coeff = 1
                p = term
            else:
                (stoich_coeff, p) = term
                assert type(stoich_coeff) == int, f"compute_reaction_quotient(): the argument `product_data` " \
                                                  f"must be a list of pairs (integer and string).  `{stoich_coeff}` is not an integer"
                assert type(p) == str, f"compute_reaction_quotient(): the argument `product_data` " \
                                       f"must be a list of pairs (integer and string).  {p} is not a string"

            species_name = p
            # TODO: Maybe turn the several next lines into a helper function
            species_conc = conc.get(species_name)
            assert species_conc is not None, f"compute_reaction_quotient(): unable to proceed because the " \
                                             f"concentration of product `{species_name}` was not provided"

            numerator *= (species_conc ** stoich_coeff)
            if explain:
                if stoich_coeff > 1:
                    numerator_text += f" [{species_name}]^{stoich_coeff} "
                else:
                    numerator_text += f"[{species_name}]"

        if explain and len(product_data) > 1:
            numerator_text = f"({numerator_text})"  # In case of multiple terms, enclose them in parenthesis


        # Compute the denominator of the "Reaction Quotient"
        for term in reactant_data:
            # Loop over the reactants
            if type(term) == str:
                stoich_coeff = 1
                r = term
            else:
                (stoich_coeff, r) = term
                assert type(stoich_coeff) == int, f"compute_reaction_quotient(): the argument `reactant_data` " \
                                                  f"must be a list of pairs (integer and string).  `{stoich_coeff}` is not an integer"
                assert type(r) == str, f"compute_reaction_quotient(): the argument `reactant_data` " \
                                       f"must be a list of pairs (integer and string).  {r} is not a string"

            species_name =  r
            # TODO: Maybe turn the several next lines into a helper function
            species_conc = conc.get(species_name)
            assert species_conc is not None, f"compute_reaction_quotient(): unable to proceed because the " \
                                             f"concentration of reactant `{species_name}` was not provided"

            denominator *= (species_conc ** stoich_coeff)
            if explain:
                if stoich_coeff > 1:
                    denominator_text += f" [{species_name}]^{stoich_coeff} "
                else:
                    denominator_text += f"[{species_name}]"

        if explain and len(reactant_data) > 1:
            denominator_text = f"({denominator_text})"  # In case of multiple terms, enclose them in parenthesis


        with np.errstate(divide='ignore', invalid='ignore'):
            # It might be np.inf (if just the denominator is zero) or np.nan (if both are zero)
            quotient = numerator / denominator

        if explain:
            formula = f"{numerator_text} / {denominator_text}"
            return (quotient, formula)

        return quotient
