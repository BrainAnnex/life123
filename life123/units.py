from dataclasses import dataclass

"""
    EARLY implementation of a unified, standardized approach to deal with units in Life123
"""


@dataclass(frozen=True)
class Unit:
    name: str       # EXAMPLE: "Celsius"
    symbol: str     # EXAMPLE: "°C"

    dimension: str  # The type of quantity.  EXAMPLES: "temperature", "molar energy"

    factor: float   # The 1st (multiplicative) step in converting to the Life123 standard internal units
    offset: float   # The 2nd (additive) step in converting to the Life123 standard internal units



def standardize_units(value :float | tuple[float, Unit], dimension=None) -> float:
    """
    If a numerical value is given, simply pass it thru.
    Otherwise, if a pair such as (25, C) is passed, convert the numeric part
    into the Life123 standard internal units

    :param value:   Either a number, or a pair such as (25, C)
    :param dimension:[OPTIONAL] If passed, it must be a standard string
                        such as "temperature" or "molar energy" - and the units (if specified)
                        are validated against it
    :return:        The value converted to the internal standard units;
                        if no units were specified, the value is assumed to already be in those units
    """
    if type(value) != tuple:
        return value

    value, units = value    # Unpack

    if dimension is not None:
        assert units.dimension == dimension, \
            f"standardize_units(): {units.name} ({units.symbol}) is not a valid unit for `{dimension}`"

    return value * units.factor + units.offset



def convert(quantity :float, from_unit :Unit, to_unit :Unit):
    """
    EXAMPLE:  convert(quantity=37, from_unit=C, to_unit=K)

    :param quantity:    The numerical value of the quantity of interest, in the given units
    :param from_unit:   Name of a "Unit" dataclass.  EXAMPLE:  C
    :param to_unit:     Name of a "Unit" dataclass.  EXAMPLE:  K
    :return:            The numerical value, converted to the new units
    """
    assert from_unit.dimension == to_unit.dimension, \
        f"convert(): cannot convert a `{from_unit.dimension}` unit ({from_unit.symbol}) " \
        f"to a `{to_unit.dimension}` unit ({to_unit.symbol})"

    standard_value = standardize_units(value=(quantity, from_unit))
    return (standard_value - to_unit.offset) / to_unit.factor





##################################################################################

# This is our internal standard unit for temperatures
K = Unit(
    name="Kelvin",
    symbol="K",

    dimension="temperature",

    # To standardize to degree Kelvin
    factor=1,
    offset=0
)

C = Unit(
    name="Celsius",
    symbol="°C",

    dimension="temperature",

    # To standardize to degree Kelvin
    factor=1,
    offset=273.15
)

F = Unit(
    name="Fahrenheit",
    symbol="°F",

    dimension="temperature",

    # To standardize to degree Kelvin
    factor=5/9,
    offset=255.3722222222222
)


# This is our internal standard unit for molar energy
KJ_PER_MOL = Unit(
    name="kiloJules per mol",
    symbol="kJ/mol",

    dimension="molar energy",

    # To standardize to kJ/mol
    factor=1,
    offset=0
)

J_PER_MOL = Unit(
    name="kiloJules per mol",
    symbol="kJ/mol",

    dimension="molar energy",

    # To standardize to kJ/mol
    factor=1/1000,
    offset=0
)
