"""
Thermodynamic formulas and calculations for gas phase species.

This module provides functions for calculating standard thermodynamic properties
of gas phase species using polynomial fits and numerical integration.
"""

from typing import List, Callable
from scipy import integrate


def calculate_standard_cp(coefficients: List[float], temperature: float,
                        R: float) -> float:
    """
    Calculate standard specific heat capacity using polynomial fit.

    The specific heat capacity is calculated using a fourth-order polynomial:
    cp°/R = Σ(aₙᵢT^(n-1)) from n=1 to 5
          = a₁ᵢ + a₂ᵢT + a₃ᵢT² + a₄ᵢT³ + a₅ᵢT⁴

    Args:
        coefficients: Polynomial coefficients [a₁, a₂, a₃, a₄, a₅]
        temperature: Temperature in K
        R: Gas constant in J/(mol·K)

    Returns:
        Standard specific heat capacity in J/(mol·K)
    """
    cp_over_R = sum(
        coeff * temperature**power
        for power, coeff in enumerate(coefficients)
    )
    return R * cp_over_R


def _create_cp_function(coefficients: List[float], R: float) -> Callable[[float],
                                                                       float]:
    """
    Create a function that calculates cp° at a given temperature.

    Args:
        coefficients: Polynomial coefficients [a₁, a₂, a₃, a₄, a₅]
        R: Gas constant in J/(mol·K)

    Returns:
        Function that takes temperature and returns cp°
    """
    return lambda T: calculate_standard_cp(coefficients, T, R)


def calculate_standard_enthalpy(coefficients: List[float], T0: float, T: float,
                              R: float) -> float:
    """
    Calculate standard state enthalpy by integrating cp° from T₀ to T.

    The enthalpy is calculated using the fundamental relation:
    H°(T) = ∫[T₀ to T] cp°(T') dT'

    Args:
        coefficients: Polynomial coefficients [a₁, a₂, a₃, a₄, a₅]
        T0: Reference temperature in K
        T: Final temperature in K
        R: Gas constant in J/(mol·K)

    Returns:
        Standard state enthalpy in J/mol
    """
    cp_function = _create_cp_function(coefficients, R)
    enthalpy, _ = integrate.quad(cp_function, T0, T)
    return enthalpy


def calculate_standard_entropy(coefficients: List[float], T0: float, T: float,
                             R: float) -> float:
    """
    Calculate standard state entropy by integrating cp°/T from T₀ to T.

    The entropy is calculated using the fundamental relation:
    S°(T) = ∫[T₀ to T] (cp°(T')/T') dT'

    Args:
        coefficients: Polynomial coefficients [a₁, a₂, a₃, a₄, a₅]
        T0: Reference temperature in K
        T: Final temperature in K
        R: Gas constant in J/(mol·K)

    Returns:
        Standard state entropy in J/(mol·K)
    """
    cp_function = _create_cp_function(coefficients, R)
    cp_over_T_function = lambda T: cp_function(T) / T
    entropy, _ = integrate.quad(cp_over_T_function, T0, T)
    return entropy 