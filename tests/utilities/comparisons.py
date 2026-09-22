# Utilities for comparisons

import math
import collections


def compare_recordsets(rs1: [{}], rs2: [{}]) -> bool:
    """
    We define "recordsets" as "lists of dictionaries".  Each element of the lists is regarded as a "record".

    EXAMPLE of recordset:  [{'Field_A': 1},
                            {'Field_A': 1},
                            {'Field_A': 99, 'Field_B': 'hello'}]

    Compare 2 recordsets WITHOUT REGARD to the position of the dictionaries within the lists,
    and also WITHOUT REGARD to the position of the key:value pairs within the dictionaries.
    Duplicates records, if present, are treated as completely separate.

    Return True if the given recordsets match, as defined above; False, otherwise.

    WARNING: this function is meant for comparing SMALL datasets, because it's Order n square!

    :param rs1: A (possibly empty) list of dictionaries
    :param rs2: A (possibly empty) list of dictionaries

    :return:    True if there's a match, or False otherwise
    """

    # Verify the type of the arguments
    assert isinstance(rs1, list), "compare_recordsets() : The 1st argument is not a list!  Value = " + str(rs1)
    assert isinstance(rs2, list), "compare_recordsets() : The 2nd argument is not a list!  Value = " + str(rs2)

    if len(rs1) != len(rs2):
        return False    # Datasets of different sizes will never match

    rs2_clone = rs2.copy()  # IMPORTANT: this is done to avoid altering the list that was passed as argument

    # Consider each element (i.e. a dictionary) in turn in the first list:
    #   attempt to remove it from the other list; if the removal fails, then it means
    #   that we have an element in the 1st list that is not present in the 2nd one (hence a mismatch)
    for rec1 in rs1:
        # Note: since Python 3.7 dictionaries are order-preserving, but
        #       built-in Python functions such as "remove"
        #       do not distinguish dictionaries based on order:
        #           {'a': 1, 'b': 2} will match {'b': 2, 'a': 1}

        try:
            rs2_clone.remove(rec1)    # Remove (the first instance of) the element rec1 from the list rs2
        except Exception:
            return False        # The remove failed - i.e. the first list contains an element not in the 2nd one

    return True



def compare_unordered_lists(l1 :dict, l2 :dict) -> bool:
    """
    Compare two lists regardless of order of elements.
    Duplicates elements, if present, are treated as completely separate.
    IMPORTANT:  the elements of the list must match exactly,
                and *must* be HASHABLE Python entities, such as
                strings, numbers or tuples; they can NOT be, for example, dictionaries

    Return True if the given lists match, as defined above; False, otherwise.

    EXAMPLES:   [1, 2, 3] will match [3, 2, 1]
                [] and [] will match
                ["x", (1, 2)] will match [(1, 2) , "x"] but NOT ["x", (2, 1)]
                ["a", "a"] will NOT match ["a"]

    :param l1:  A list of HASHABLE Python entities (e.g. strings or numbers)
    :param l2:  Same as above
    :return:    True if there's a match, or False otherwise
    """
    return collections.Counter(l1) == collections.Counter(l2)



def compare_dicts(d1 :dict, d2 :dict, **kwargs) -> bool:
    """
    Return True if the two dictionaries have the same keys,
    and approximately equality of their corresponding values

    :param d1:
    :param d2:
    :param kwargs:  [OPTIONAL] Extra arguments of math.isclose(),
                        such as rel_tol or abs_tol
    :return:
    """
    if d1.keys() != d2.keys():
        return False

    for key, value in d1.items():
        other_value = d2[key]
        if (value is None) or (other_value is None):
            if value != other_value:
                return False

        elif not math.isclose(value, other_value, **kwargs):
            return  False

    return True
    """
    return d1.keys() == d2.keys() and all(
        math.isclose(d1[k], d2[k], **kwargs) for k in d1
    )
    """
