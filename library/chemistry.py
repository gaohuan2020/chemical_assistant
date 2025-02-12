class Substance(object):
    """ Class representing a chemical substance

    Parameters
    ----------
    name : str
    charge : int (optional, default: None)
        Will be stored in composition[0], prefer composition when possible.
    latex_name : str
    unicode_name : str
    html_name : str
    composition : dict or None (default)
        Dictionary (int -> number) e.g. {atomic number: count}, zero has special
        meaning (net charge). Avoid using the key 0 unless you specifically mean
        net charge. The motivation behind this is that it is easier to track a
        net-charge of e.g. 6 for U(VI) than it is to remember that uranium has 92
        electrons and use 86 as the value).
    data : dict
        Free form dictionary. Could be simple such as ``{'mp': 0, 'bp': 100}``
        or considerably more involved, e.g.: ``{'diffusion_coefficient': {\
 'water': lambda T: 2.1*m**2/s/K*(T - 273.15*K)}}``.

    Attributes
    ----------
    mass
        Maps to data['mass'], and when unavailable looks for ``formula.mass``.
    attrs
        A tuple of attribute names for serialization.
    composition : dict or None
        Dictionary mapping fragment key (str) to amount (int).
    data
        Free form dictionary.

    Examples
    --------
    >>> ammonium = Substance('NH4+', 1, 'NH_4^+', composition={7: 1, 1: 4},
    ...     data={'mass': 18.0385, 'pKa': 9.24})
    >>> ammonium.name
    'NH4+'
    >>> ammonium.composition == {0: 1, 1: 4, 7: 1}  # charge represented by key '0'
    True
    >>> ammonium.data['mass']
    18.0385
    >>> ammonium.data['pKa']
    9.24
    >>> ammonium.mass  # mass is a special case (also attribute)
    18.0385
    >>> ammonium.pKa
    Traceback (most recent call last):
        ...
    AttributeError: 'Substance' object has no attribute 'pKa'
    >>> nh4p = Substance.from_formula('NH4+')  # simpler
    >>> nh4p.composition == {7: 1, 1: 4, 0: 1}
    True
    >>> nh4p.latex_name
    'NH_{4}^{+}'

    """

    attrs = ("name", "latex_name", "unicode_name", "html_name", "composition", "data")
