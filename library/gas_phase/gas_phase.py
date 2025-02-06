from typing import Dict, List
import yaml
import os
from .data_types import Species, ParticleType
from . import thermodynamics


class GasPhase:
    """Represents the gas phase (Sg) containing multiple species"""

    # Load constants from config file
    config_path = os.path.join(os.path.dirname(__file__), 'config.yaml')
    with open(config_path, 'r', encoding='utf-8') as f:
        _config = yaml.safe_load(f)
        _constants = _config['constants']
        R = _constants['gas_constant']
        AMU_TO_KG_MOL = _constants['amu_to_kg_mol']
        _STANDARD_TEMP = _constants['standard_temperature']
        _STANDARD_PRESSURE = _constants['standard_pressure']
        _CP_COEFFICIENTS = _config['cp_coefficients']

    def __init__(self,
                 temperature: float = _STANDARD_TEMP,
                 pressure: float = _STANDARD_PRESSURE):
        self.species: Dict[str, Species] = {}
        self.temperature = temperature  # K
        self.pressure = pressure  # Pa

    def add_species(self, species: Species) -> None:
        """Add a species to the gas phase"""
        self.species[species.name] = species

    def get_total_concentration(self) -> float:
        """Calculate total concentration of all species in the gas phase"""
        return sum(species.concentration for species in self.species.values())

    def get_species_by_type(self,
                            particle_type: ParticleType) -> List[Species]:
        """Get all species of a specific type"""
        return [
            s for s in self.species.values()
            if s.particle_type == particle_type
        ]

    def get_mole_fraction(self, species_name: str) -> float:
        """Calculate mole fraction of a specific species"""
        if species_name not in self.species:
            raise KeyError(f"Species {species_name} not found")
        return (self.species[species_name].concentration /
                self.get_total_concentration())

    def get_mean_molar_mass(self) -> float:
        """Calculate the mean molar mass of the gas mixture"""
        numerator = sum(species.concentration * species.molecular_weight
                        for species in self.species.values())
        return numerator / self.get_total_concentration()

    def get_mass_fraction(self, species_name: str) -> float:
        """Calculate mass fraction of a specific species"""
        if species_name not in self.species:
            raise KeyError(f"Species {species_name} not found")

        species = self.species[species_name]
        mean_molar_mass = self.get_mean_molar_mass()
        mole_fraction = self.get_mole_fraction(species_name)

        return (species.molecular_weight / mean_molar_mass) * mole_fraction

    def get_density(self) -> float:
        """Calculate the density of the gas mixture using the ideal gas law"""
        mean_molar_mass = self.get_mean_molar_mass()
        mean_molar_mass_kg = mean_molar_mass * self.AMU_TO_KG_MOL

        return (self.pressure * mean_molar_mass_kg) / (self.R *
                                                       self.temperature)

    def get_heat_capacity(self, delta_t: float = 0.1) -> float:
        """
        Calculate the heat capacity (C) of the gas mixture, defined as δQ/dT.
        Uses numerical differentiation to calculate the rate of heat change 
        with temperature.
        
        Args:
            delta_t: Temperature difference for numerical derivative (default 0.1 K)
            
        Returns:
            Heat capacity in J/K
        """
        # Store original temperature
        original_temp = self.temperature

        # Calculate heat energy at T + ΔT
        self.temperature = original_temp + delta_t
        q1 = sum(species.concentration * species.heat_capacity *
                 self.temperature for species in self.species.values())

        # Calculate heat energy at T - ΔT
        self.temperature = original_temp - delta_t
        q2 = sum(species.concentration * species.heat_capacity *
                 self.temperature for species in self.species.values())

        # Restore original temperature
        self.temperature = original_temp

        # Calculate C using central difference approximation
        # C = δQ/dT ≈ (Q(T+ΔT) - Q(T-ΔT))/(2ΔT)
        return (q1 - q2) / (2 * delta_t)

    def get_enthalpy(self) -> float:
        """
        Calculate the molar enthalpy of the gas mixture.
        Returns enthalpy in J/mol
        """
        total_concentration = self.get_total_concentration()
        if total_concentration == 0:
            return 0.0

        # Calculate molar enthalpy as weighted sum of species enthalpies
        return sum(species.heat_capacity * self.temperature *
                   (species.concentration / total_concentration)
                   for species in self.species.values())

    def get_specific_heat_capacity(self, delta_t: float = 0.1) -> float:
        """
        Calculate the specific heat capacity at constant pressure (cp)
        using the derivative of enthalpy with respect to temperature: cp = (∂h/∂T)_p
        
        Args:
            delta_t: Temperature difference for numerical derivative (default 0.1 K)
            
        Returns:
            Specific heat capacity in J/(mol·K)
        """
        # Store original temperature
        original_temp = self.temperature

        # Calculate enthalpy at T + ΔT
        self.temperature = original_temp + delta_t
        h1 = self.get_enthalpy()

        # Calculate enthalpy at T - ΔT
        self.temperature = original_temp - delta_t
        h2 = self.get_enthalpy()

        # Restore original temperature
        self.temperature = original_temp

        # Calculate cp using central difference approximation
        # cp = (∂h/∂T)_p ≈ (h(T+ΔT) - h(T-ΔT))/(2ΔT)
        return (h1 - h2) / (2 * delta_t)

    def get_specific_heat_capacity_v(self) -> float:
        """
        Calculate the specific heat capacity at constant volume (cv)
        as sum of translational, rotational, and vibrational contributions:
        cv = cv,tr + cv,rot + cv,vib
        
        Returns:
            Specific heat capacity at constant volume in J/(mol·K)
        """
        cv_tr = self.get_translational_heat_capacity_v()
        cv_rot = self.get_rotational_heat_capacity_v()
        cv_vib = self.get_vibrational_heat_capacity_v()

        return cv_tr + cv_rot + cv_vib

    def get_translational_heat_capacity_v(self) -> float:
        """
        Calculate translational contribution to cv.
        For ideal gas: cv,tr = 3/2 R
        
        Returns:
            Translational heat capacity in J/(mol·K)
        """
        return 1.5 * self.R

    def get_rotational_heat_capacity_v(self) -> float:
        """
        Calculate rotational contribution to cv.
        cv,rot = f/2 R where f is the number of excited rotational degrees of freedom:
        f = 0 for atoms
        f = 2 for linear molecules
        f = 3 for non-linear molecules
        
        Returns:
            Rotational heat capacity in J/(mol·K)
        """
        total_concentration = self.get_total_concentration()
        if total_concentration == 0:
            return 0.0

        # Sum up rotational contributions weighted by mole fraction
        cv_rot = 0.0
        for species in self.species.values():
            mole_fraction = species.concentration / total_concentration

            # Determine degrees of freedom based on particle type
            if species.particle_type == ParticleType.ATOMIC:
                f = 0  # Atoms have no rotational modes
            elif species.particle_type == ParticleType.MOLECULAR:
                f = 2  # Linear molecules have 2 rotational degrees of freedom
            else:
                f = 3  # Non-linear molecules have 3 rotational degrees of freedom

            cv_rot += (f / 2) * self.R * mole_fraction

        return cv_rot

    def get_vibrational_heat_capacity_v(self) -> float:
        """
        Calculate vibrational contribution to cv.
        This is the remainder after subtracting translational and rotational parts
        from the total cv.
        
        Returns:
            Vibrational heat capacity in J/(mol·K)
        """
        cv_total = self.get_specific_heat_capacity() - self.R  # cv = cp - R
        cv_tr = self.get_translational_heat_capacity_v()
        cv_rot = self.get_rotational_heat_capacity_v()

        return cv_total - cv_tr - cv_rot

    def get_entropy_change(self, delta_t: float = 0.1) -> float:
        """
        Calculate entropy change (dS) using the second law of thermodynamics:
        dS = δQ/T
        
        Args:
            delta_t: Temperature difference for heat calculation (default 0.1 K)
            
        Returns:
            Entropy change in J/K
        """
        # Store original temperature
        original_temp = self.temperature

        # Calculate heat change (δQ)
        self.temperature = original_temp + delta_t
        q1 = sum(species.concentration * species.heat_capacity *
                 self.temperature for species in self.species.values())

        self.temperature = original_temp - delta_t
        q2 = sum(species.concentration * species.heat_capacity *
                 self.temperature for species in self.species.values())

        # Restore original temperature
        self.temperature = original_temp

        # Calculate δQ
        delta_q = q1 - q2

        # Calculate entropy change: dS = δQ/T
        return delta_q / self.temperature

    def get_standard_cp(self, species_name: str, temperature: float) -> float:
        """
        Calculate standard specific heat capacity at constant pressure.
        
        Args:
            species_name: Name of the species
            temperature: Temperature in K
            
        Returns:
            Standard specific heat capacity in J/(mol·K)
        """
        if species_name not in self.species:
            raise KeyError(f"Species {species_name} not found")
            
        if species_name not in self._CP_COEFFICIENTS:
            raise ValueError(f"No polynomial coefficients available for {species_name}")
            
        coefficients = self._CP_COEFFICIENTS[species_name]
        return thermodynamics.calculate_standard_cp(coefficients, temperature, self.R)

    def get_standard_enthalpy(self, species_name: str, temperature: float) -> float:
        """
        Calculate standard state enthalpy.
        
        Args:
            species_name: Name of the species
            temperature: Temperature in K
            
        Returns:
            Standard state enthalpy in J/mol
        """
        if species_name not in self.species:
            raise KeyError(f"Species {species_name} not found")
            
        if species_name not in self._CP_COEFFICIENTS:
            raise ValueError(f"No polynomial coefficients available for {species_name}")
            
        coefficients = self._CP_COEFFICIENTS[species_name]
        return thermodynamics.calculate_standard_enthalpy(
            coefficients, self._STANDARD_TEMP, temperature, self.R)

    def get_standard_entropy(self, species_name: str, temperature: float) -> float:
        """
        Calculate standard state entropy.
        
        Args:
            species_name: Name of the species
            temperature: Temperature in K
            
        Returns:
            Standard state entropy in J/(mol·K)
        """
        if species_name not in self.species:
            raise KeyError(f"Species {species_name} not found")
            
        if species_name not in self._CP_COEFFICIENTS:
            raise ValueError(f"No polynomial coefficients available for {species_name}")
            
        coefficients = self._CP_COEFFICIENTS[species_name]
        return thermodynamics.calculate_standard_entropy(
            coefficients, self._STANDARD_TEMP, temperature, self.R)
