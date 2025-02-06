"""
Thermodynamic formulas and calculations for gas phase species.

This module provides functions for calculating standard thermodynamic properties
of gas phase species using polynomial fits and numerical integration.
"""

from typing import List, Callable
from scipy import integrate
from math import log


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


def calculate_standard_enthalpy(coefficients: List[float], T: float,
                              R: float, H_298: float) -> float:
    """
    Calculate standard state enthalpy using polynomial integration.

    The enthalpy is calculated using the analytical form:
    H°ᵢ = R(a₁ᵢT + a₂ᵢT²/2 + a₃ᵢT³/3 + a₄ᵢT⁴/4 + a₅ᵢT⁵/5 + a₆ᵢ)
    where a₆ᵢ is determined from H°₂₉₈

    Args:
        coefficients: Polynomial coefficients [a₁, a₂, a₃, a₄, a₅]
        T: Temperature in K
        R: Gas constant in J/(mol·K)
        H_298: Standard heat of formation at 298K in J/mol

    Returns:
        Standard state enthalpy in J/mol
    """
    # Calculate H°(T)/R without a₆ᵢ
    H_over_R = sum(
        coeff * T**(power + 1) / (power + 1)
        for power, coeff in enumerate(coefficients)
    )
    
    # Calculate H°(298)/R without a₆ᵢ
    H_298_over_R = sum(
        coeff * 298.15**(power + 1) / (power + 1)
        for power, coeff in enumerate(coefficients)
    )
    
    # Solve for a₆ᵢ using H°₂₉₈
    a6 = (H_298 - R * H_298_over_R) / R
    
    return R * (H_over_R + a6)


def calculate_standard_entropy(coefficients: List[float], T: float,
                             R: float, S_298: float) -> float:
    """
    Calculate standard state entropy using polynomial integration.

    The entropy is calculated using the analytical form:
    S°ᵢ = R(a₁ᵢlnT + a₂ᵢT + a₃ᵢT²/2 + a₄ᵢT³/3 + a₅ᵢT⁴/4 + a₇ᵢ)
    where a₇ᵢ is determined from S°₂₉₈

    Args:
        coefficients: Polynomial coefficients [a₁, a₂, a₃, a₄, a₅]
        T: Temperature in K
        R: Gas constant in J/(mol·K)
        S_298: Standard entropy at 298K in J/(mol·K)

    Returns:
        Standard state entropy in J/(mol·K)
    """
    # Calculate S°(T)/R without a₇ᵢ
    S_over_R = coefficients[0] * log(T)  # First term (a₁lnT)
    S_over_R += sum(
        coeff * T**power / power
        for power, coeff in enumerate(coefficients[1:], 1)
    )
    
    # Calculate S°(298)/R without a₇ᵢ
    S_298_over_R = coefficients[0] * log(298.15)  # First term at 298.15K
    S_298_over_R += sum(
        coeff * 298.15**power / power
        for power, coeff in enumerate(coefficients[1:], 1)
    )
    
    # Solve for a₇ᵢ using S°₂₉₈
    a7 = (S_298 - R * S_298_over_R) / R
    
    return R * (S_over_R + a7) 

