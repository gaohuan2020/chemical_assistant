from dataclasses import dataclass
from typing import Dict
from enum import Enum


class ParticleType(Enum):
    """Classification of gas species"""
    ATOMIC = "atomic"
    MOLECULAR = "molecular"
    RADICAL = "radical"
    ION = "ion"


class ElectronicState(Enum):
    """Common electronic states"""
    GROUND = "ground"
    EXCITED = "excited"
    IONIZED = "ionized"


@dataclass(frozen=True)
class Atom:
    """Basic atomic properties"""
    symbol: str
    atomic_number: int
    mass_number: int


@dataclass
class Species:
    """Represents a gas-phase species with its composition and properties"""
    name: str
    molecular_weight: float
    particle_type: ParticleType
    concentration: float = 0.0
    heat_capacity: float = 0.0      # J/(mol·K)
    formation_enthalpy: float = 0.0 # J/mol at standard conditions

    @property
    def atom_count(self) -> int:
        """Get total number of atoms in the species"""
        return sum(self.composition.values()) 