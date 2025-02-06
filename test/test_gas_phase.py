from library.gas_phase.data_types import Atom, Species, ParticleType, ElectronicState
from library.gas_phase import GasPhase


def create_test_species():
    """Create test species for the examples"""
    # Define atoms
    H = Atom("H", atomic_number=1, mass_number=1)
    O = Atom("O", atomic_number=8, mass_number=16)
    N = Atom("N", atomic_number=7, mass_number=14)

    # Create species
    species = [
        Species(name="O",
                molecular_weight=16.0,  # AMU
                particle_type=ParticleType.ATOMIC,
                concentration=0.001,
                heat_capacity=20.8),  # J/(mol·K)
        Species(name="O2",
                molecular_weight=32.0,  # AMU
                particle_type=ParticleType.MOLECULAR,
                concentration=0.20946,
                heat_capacity=29.4),  # J/(mol·K)
        Species(name="OH",
                molecular_weight=17.0,  # AMU
                particle_type=ParticleType.RADICAL,
                concentration=0.0001,
                heat_capacity=29.9)   # J/(mol·K)
    ]
    return species


def test_gas_phase():
    """Run all tests for the gas phase calculations"""
    # Initialize
    species_list = create_test_species()
    gas_phase = GasPhase()
    for species in species_list:
        gas_phase.add_species(species)

    # Test total concentration
    test_total_concentration(gas_phase)

    # Test mean molar mass
    test_mean_molar_mass(gas_phase)

    # Test mole fractions
    test_mole_fractions(gas_phase)

    # Test mass fractions
    test_mass_fractions(gas_phase)

    # Test density
    test_density(gas_phase)

    # Test heat capacity
    test_heat_capacity(gas_phase)

    # Test specific heat capacity
    test_specific_heat_capacity(gas_phase)

    # Test specific heat capacity at constant volume
    test_specific_heat_capacity_v(gas_phase)
    
    # Test cv components
    test_heat_capacity_v_components(gas_phase)
    
    # Test entropy change
    test_entropy_change(gas_phase)
    
    # Test standard cp polynomial
    test_standard_cp_polynomial(gas_phase)
    
    # Test standard enthalpy
    test_standard_enthalpy(gas_phase)
    
    # Test standard entropy
    test_standard_entropy(gas_phase)
    
    # Print composition summary
    print_composition_summary(gas_phase)


def test_total_concentration(gas_phase):
    """Test the total concentration calculation"""
    expected_total = sum(s.concentration for s in gas_phase.species.values())
    calculated_total = gas_phase.get_total_concentration()

    print("\nTesting total concentration calculation:")
    print(f"Expected total concentration: {expected_total:.6f}")
    print(f"Calculated total concentration: {calculated_total:.6f}")
    print(f"Difference: {abs(expected_total - calculated_total):.10f}")

    assert abs(expected_total - calculated_total) < 1e-10


def test_mean_molar_mass(gas_phase):
    """Test the mean molar mass calculation"""
    # Calculate expected mean molar mass manually
    total_concentration = gas_phase.get_total_concentration()
    expected_mean = sum(
        species.concentration * species.molecular_weight
        for species in gas_phase.species.values()) / total_concentration

    calculated_mean = gas_phase.get_mean_molar_mass()

    print("\nTesting mean molar mass calculation:")
    print(f"Expected mean molar mass: {expected_mean:.6f}")
    print(f"Calculated mean molar mass: {calculated_mean:.6f}")
    print(f"Difference: {abs(expected_mean - calculated_mean):.10f}")

    assert abs(expected_mean - calculated_mean) < 1e-10


def test_mole_fractions(gas_phase):
    """Test the mole fraction calculations"""
    total_concentration = gas_phase.get_total_concentration()

    print("\nTesting mole fraction calculations:")
    for species_name, species in gas_phase.species.items():
        expected_fraction = species.concentration / total_concentration
        calculated_fraction = gas_phase.get_mole_fraction(species_name)

        print(f"{species_name}:")
        print(f"Expected mole fraction: {expected_fraction:.6f}")
        print(f"Calculated mole fraction: {calculated_fraction:.6f}")
        print(
            f"Difference: {abs(expected_fraction - calculated_fraction):.10f}")

        assert abs(expected_fraction - calculated_fraction) < 1e-10


def test_mass_fractions(gas_phase):
    """Test the mass fraction calculations"""
    mean_molar_mass = gas_phase.get_mean_molar_mass()

    print("\nTesting mass fraction calculations:")
    for species_name, species in gas_phase.species.items():
        mole_fraction = gas_phase.get_mole_fraction(species_name)
        expected_fraction = (species.molecular_weight /
                             mean_molar_mass) * mole_fraction
        calculated_fraction = gas_phase.get_mass_fraction(species_name)

        print(f"{species_name}:")
        print(f"Expected mass fraction: {expected_fraction:.6f}")
        print(f"Calculated mass fraction: {calculated_fraction:.6f}")
        print(
            f"Difference: {abs(expected_fraction - calculated_fraction):.10f}")

        assert abs(expected_fraction - calculated_fraction) < 1e-10


def test_density(gas_phase):
    """Test the density calculation"""
    mean_molar_mass = gas_phase.get_mean_molar_mass()
    mean_molar_mass_kg = mean_molar_mass * gas_phase.AMU_TO_KG_MOL
    expected_density = (gas_phase.pressure * mean_molar_mass_kg) / (
        gas_phase.R * gas_phase.temperature)

    calculated_density = gas_phase.get_density()

    print("\nTesting density calculation:")
    print(f"Expected density: {expected_density:.6f}")
    print(f"Calculated density: {calculated_density:.6f}")
    print(f"Difference: {abs(expected_density - calculated_density):.10f}")

    assert abs(expected_density - calculated_density) < 1e-10


def test_heat_capacity(gas_phase):
    """Test the heat capacity calculation"""
    print("\nTesting heat capacity calculation:")
    
    # For a small temperature change dT
    dT = 0.1  # K
    original_temp = gas_phase.temperature
    
    # Calculate heat change (δQ) directly
    # Q = sum(n * cp * T) for each species
    gas_phase.temperature = original_temp + dT
    q1 = sum(
        species.concentration * species.heat_capacity * gas_phase.temperature
        for species in gas_phase.species.values()
    )
    
    gas_phase.temperature = original_temp - dT
    q2 = sum(
        species.concentration * species.heat_capacity * gas_phase.temperature
        for species in gas_phase.species.values()
    )
    
    # Reset temperature
    gas_phase.temperature = original_temp
    
    # Expected heat capacity from definition C = δQ/dT
    expected_capacity = (q1 - q2) / (2 * dT)
    
    # Calculate using the method
    calculated_capacity = gas_phase.get_heat_capacity()
    
    print(f"Expected heat capacity: {expected_capacity:.6f} J/K")
    print(f"Calculated heat capacity: {calculated_capacity:.6f} J/K")
    print(f"Difference: {abs(expected_capacity - calculated_capacity):.10f}")
    
    assert abs(expected_capacity - calculated_capacity) < 1e-10


def test_specific_heat_capacity(gas_phase):
    """Test the specific heat capacity (cp) calculation"""
    print("\nTesting specific heat capacity calculation:")
    
    # For a small temperature change dT
    dT = 0.1  # K
    original_temp = gas_phase.temperature
    
    # Calculate enthalpy change directly
    gas_phase.temperature = original_temp + dT
    h1 = gas_phase.get_enthalpy()
    
    gas_phase.temperature = original_temp - dT
    h2 = gas_phase.get_enthalpy()
    
    # Reset temperature
    gas_phase.temperature = original_temp
    
    # Expected cp from definition cp = (∂h/∂T)_p
    expected_cp = (h1 - h2) / (2 * dT)
    
    # Calculate using the method
    calculated_cp = gas_phase.get_specific_heat_capacity()
    
    print(f"Expected cp: {expected_cp:.6f} J/(mol·K)")
    print(f"Calculated cp: {calculated_cp:.6f} J/(mol·K)")
    print(f"Difference: {abs(expected_cp - calculated_cp):.10f}")
    
    assert abs(expected_cp - calculated_cp) < 1e-10


def test_specific_heat_capacity_v(gas_phase):
    """Test the specific heat capacity at constant volume (cv) calculation"""
    print("\nTesting specific heat capacity at constant volume calculation:")
    
    # Get cp
    cp = gas_phase.get_specific_heat_capacity()
    
    # Expected cv from ideal gas relation cv = cp - R
    expected_cv = cp - gas_phase.R
    
    # Calculate using the method
    calculated_cv = gas_phase.get_specific_heat_capacity_v()
    
    print(f"cp: {cp:.6f} J/(mol·K)")
    print(f"R: {gas_phase.R:.6f} J/(mol·K)")
    print(f"Expected cv: {expected_cv:.6f} J/(mol·K)")
    print(f"Calculated cv: {calculated_cv:.6f} J/(mol·K)")
    print(f"Difference: {abs(expected_cv - calculated_cv):.10f}")
    
    assert abs(expected_cv - calculated_cv) < 1e-10


def test_heat_capacity_v_components(gas_phase):
    """Test the components of specific heat capacity at constant volume"""
    print("\nTesting cv components:")
    
    # Get all components
    cv_tr = gas_phase.get_translational_heat_capacity_v()
    cv_rot = gas_phase.get_rotational_heat_capacity_v()
    cv_vib = gas_phase.get_vibrational_heat_capacity_v()
    cv_total = gas_phase.get_specific_heat_capacity_v()
    
    # Print components
    print(f"Translational (cv,tr): {cv_tr:.6f} J/(mol·K)")
    print(f"Rotational (cv,rot): {cv_rot:.6f} J/(mol·K)")
    print(f"Vibrational (cv,vib): {cv_vib:.6f} J/(mol·K)")
    print(f"Total cv: {cv_total:.6f} J/(mol·K)")
    
    # Test 1: Verify cv = cv,tr + cv,rot + cv,vib
    expected_cv = cv_tr + cv_rot + cv_vib
    print(f"\nSum of components: {expected_cv:.6f} J/(mol·K)")
    print(f"Difference from total: {abs(cv_total - expected_cv):.10f}")
    assert abs(cv_total - expected_cv) < 1e-10
    
    # Test 2: Verify cv,tr = 3/2 R
    expected_cv_tr = 1.5 * gas_phase.R
    print(f"\nExpected cv,tr (3/2 R): {expected_cv_tr:.6f} J/(mol·K)")
    print(f"Difference: {abs(cv_tr - expected_cv_tr):.10f}")
    assert abs(cv_tr - expected_cv_tr) < 1e-10
    
    # Test 3: Verify cv,rot for each species type
    print("\nTesting rotational contributions:")
    for species_name, species in gas_phase.species.items():
        if species.particle_type == ParticleType.ATOMIC:
            f = 0
        elif species.particle_type == ParticleType.MOLECULAR:
            f = 2
        else:
            f = 3
        print(f"{species_name} (f={f}): {f/2 * gas_phase.R:.6f} J/(mol·K)")


def test_entropy_change(gas_phase):
    """Test the entropy change calculation"""
    print("\nTesting entropy change calculation:")
    
    # For a small temperature change
    dT = 0.1  # K
    original_temp = gas_phase.temperature
    
    # Calculate heat change directly
    gas_phase.temperature = original_temp + dT
    q1 = sum(
        species.concentration * species.heat_capacity * gas_phase.temperature
        for species in gas_phase.species.values()
    )
    
    gas_phase.temperature = original_temp - dT
    q2 = sum(
        species.concentration * species.heat_capacity * gas_phase.temperature
        for species in gas_phase.species.values()
    )
    
    # Reset temperature
    gas_phase.temperature = original_temp
    
    # Expected entropy change from dS = δQ/T
    delta_q = q1 - q2
    expected_ds = delta_q / gas_phase.temperature
    
    # Calculate using the method
    calculated_ds = gas_phase.get_entropy_change()
    
    print(f"Temperature: {gas_phase.temperature:.2f} K")
    print(f"Heat change (δQ): {delta_q:.6f} J")
    print(f"Expected dS: {expected_ds:.6f} J/K")
    print(f"Calculated dS: {calculated_ds:.6f} J/K")
    print(f"Difference: {abs(expected_ds - calculated_ds):.10f}")
    
    assert abs(expected_ds - calculated_ds) < 1e-10


def test_standard_cp_polynomial(gas_phase):
    """Test the polynomial fit calculation of standard cp"""
    print("\nTesting standard cp polynomial fit:")
    
    # Test temperatures
    temperatures = [300, 500, 1000, 1500, 2000]  # K
    
    # Expected values at T = 300K (calculated from polynomial)
    expected_values = {
        "O": 3.025 * gas_phase.R + 1.639e-3 * gas_phase.R * 300 + 
             -1.138e-6 * gas_phase.R * 300**2 + 3.331e-10 * gas_phase.R * 300**3 + 
             -3.708e-14 * gas_phase.R * 300**4,
        "O2": 3.784 * gas_phase.R + -2.996e-3 * gas_phase.R * 300 + 
              9.847e-6 * gas_phase.R * 300**2 + -9.681e-9 * gas_phase.R * 300**3 + 
              3.243e-12 * gas_phase.R * 300**4,
        "OH": 3.637 * gas_phase.R + -1.794e-3 * gas_phase.R * 300 + 
              3.707e-6 * gas_phase.R * 300**2 + -3.142e-9 * gas_phase.R * 300**3 + 
              9.758e-13 * gas_phase.R * 300**4
    }
    
    # Test at different temperatures
    for T in temperatures:
        print(f"\nAt T = {T} K:")
        for species_name in gas_phase.species:
            cp = gas_phase.get_standard_cp(species_name, T)
            print(f"{species_name}: Cp° = {cp:.3f} J/(mol·K)")
            
            if T == 300:  # Verify against expected values at 300K
                print(f"Expected at 300K: {expected_values[species_name]:.3f} J/(mol·K)")
                assert abs(cp - expected_values[species_name]) < 1e-10


def test_standard_enthalpy(gas_phase):
    """Test the standard state enthalpy calculation"""
    print("\nTesting standard state enthalpy calculation:")
    
    # Test at different temperatures
    temperatures = [300, 500, 1000, 1500, 2000]  # K
    T0 = gas_phase._STANDARD_TEMP
    
    for T in temperatures:
        print(f"\nAt T = {T} K:")
        for species_name in gas_phase.species:
            # Calculate using the method which uses numerical integration
            calculated_H = gas_phase.get_standard_enthalpy(species_name, T)
            
            # Calculate expected value using numerical integration for verification
            from scipy import integrate
            
            def cp_function(T):
                coeffs = gas_phase._CP_COEFFICIENTS[species_name]
                cp_over_R = sum(a * T**n for n, a in enumerate(coeffs))
                return gas_phase.R * cp_over_R
            
            expected_H, _ = integrate.quad(cp_function, T0, T)
            
            print(f"{species_name}:")
            print(f"  H° = {calculated_H:.3f} J/mol")
            print(f"  Expected H° = {expected_H:.3f} J/mol")
            print(f"  Difference: {abs(calculated_H - expected_H):.10f}")
            
            assert abs(calculated_H - expected_H) < 1e-10


def test_standard_entropy(gas_phase):
    """Test the standard state entropy calculation"""
    print("\nTesting standard state entropy calculation:")
    
    # Test at different temperatures
    temperatures = [300, 500, 1000, 1500, 2000]  # K
    T0 = gas_phase._STANDARD_TEMP
    
    for T in temperatures:
        print(f"\nAt T = {T} K:")
        for species_name in gas_phase.species:
            # Calculate using the method which uses numerical integration
            calculated_S = gas_phase.get_standard_entropy(species_name, T)
            
            # Calculate expected value using numerical integration for verification
            from scipy import integrate
            
            def cp_over_T_function(T):
                coeffs = gas_phase._CP_COEFFICIENTS[species_name]
                cp_over_R = sum(a * T**n for n, a in enumerate(coeffs))
                return gas_phase.R * cp_over_R / T
            
            expected_S, _ = integrate.quad(cp_over_T_function, T0, T)
            
            print(f"{species_name}:")
            print(f"  S° = {calculated_S:.3f} J/(mol·K)")
            print(f"  Expected S° = {expected_S:.3f} J/(mol·K)")
            print(f"  Difference: {abs(calculated_S - expected_S):.10f}")
            
            assert abs(calculated_S - expected_S) < 1e-10


def print_composition_summary(gas_phase):
    """Print a summary of the gas phase composition"""
    print("\nGas Phase Composition Summary:")
    print("-" * 50)
    print(
        f"{'Species':<10} {'Conc.':<12} {'Mole Frac.':<12} {'Mass Frac.':<12}")
    print("-" * 50)

    for species_name in gas_phase.species:
        concentration = gas_phase.species[species_name].concentration
        mole_fraction = gas_phase.get_mole_fraction(species_name)
        mass_fraction = gas_phase.get_mass_fraction(species_name)

        print(
            f"{species_name:<10} {concentration:<12.6f} {mole_fraction:<12.6f} "
            f"{mass_fraction:<12.6f}")

    print("-" * 50)
    print(f"Total concentration: {gas_phase.get_total_concentration():.6f}")
    print(f"Mean molar mass: {gas_phase.get_mean_molar_mass():.6f}")
    print(f"Density: {gas_phase.get_density():.6f}")


if __name__ == "__main__":
    test_gas_phase()
