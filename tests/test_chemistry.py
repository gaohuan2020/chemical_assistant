import pytest
from library.chemistry import Substance, Species

def test_substance_basic_properties():
    """Test basic properties of Substance class"""
    # Test creation with basic properties
    nh4 = Substance(
        name="NH4+",
        charge=1,
        latex_name="NH_4^+",
        unicode_name="NH₄⁺",
        html_name="NH<sub>4</sub><sup>+</sup>",
        composition={7: 1, 1: 4},
        data={'mass': 18.0385, 'pKa': 9.24}
    )
    
    # Test name properties
    assert nh4.name == "NH4+"
    assert nh4.latex_name == "NH_4^+"
    assert nh4.unicode_name == "NH₄⁺"
    assert nh4.html_name == "NH<sub>4</sub><sup>+</sup>"
    
    # Test composition
    expected_composition = {0: 1, 1: 4, 7: 1}  # 0 represents charge
    assert nh4.composition == expected_composition
    
    # Test data access
    assert nh4.data['mass'] == 18.0385
    assert nh4.data['pKa'] == 9.24
    assert nh4.mass == 18.0385  # Test mass property

def test_substance_from_formula():
    """Test creation of Substance from formula"""
    nh4 = Substance.from_formula('NH4+')
    
    # Test composition
    expected_composition = {0: 1, 1: 4, 7: 1}
    assert nh4.composition == expected_composition
    
    # Test generated latex name
    assert nh4.latex_name == 'NH_{4}^{+}'

def test_substance_invalid_access():
    """Test invalid property access"""
    substance = Substance(
        name="Test",
        data={'mass': 10.0}
    )
    
    # Should raise AttributeError when accessing non-existent property
    with pytest.raises(AttributeError):
        substance.non_existent_property

    # Should raise KeyError when accessing non-existent data
    with pytest.raises(KeyError):
        substance.data['non_existent_key']

def test_substance_without_optional_params():
    """Test creation of Substance without optional parameters"""
    substance = Substance(name="Simple")
    
    assert substance.name == "Simple"
    assert substance.composition is None
    assert substance.data == {}

def test_species_basic_properties():
    """Test basic properties of Species class"""
    # Test creation with basic properties
    water = Species(
        name="H2O",
        composition={1: 2, 8: 1},
        phase_idx=0,  # default phase (usually aqueous or liquid)
        data={'mass': 18.015}
    )
    
    assert water.name == "H2O"
    assert water.composition == {1: 2, 8: 1}
    assert water.phase_idx == 0
    assert water.data['mass'] == 18.015

def test_species_from_formula():
    """Test creation of Species from formula with different phases"""
    # Test aqueous species (default phase)
    h2o_aq = Species.from_formula('H2O')
    assert h2o_aq.phase_idx == 0
    
    # Test solid phase
    nacl_s = Species.from_formula('NaCl(s)')
    assert nacl_s.phase_idx == 1
    assert nacl_s.composition == {11: 1, 17: 1}
    
    # Test liquid phase
    hg_l = Species.from_formula('Hg(l)')
    assert hg_l.phase_idx == 2
    
    # Test gas phase
    co2_g = Species.from_formula('CO2(g)')
    assert co2_g.phase_idx == 3
    assert co2_g.composition == {6: 1, 8: 2}

def test_species_custom_phases():
    """Test Species creation with custom phase definitions"""
    # Test with custom phase list
    custom_phases = ['(aq)', '(s)', '(l)', '(g)']
    water = Species.from_formula('H2O(aq)', phases=custom_phases)
    assert water.phase_idx == 1  # phase_idx is 1-based for list phases
    
    # Test with phase dictionary
    phase_dict = {'(aq)': 0, '(s)': 1, '(l)': 2, '(g)': 3}
    water = Species.from_formula('H2O(aq)', phases=phase_dict)
    assert water.phase_idx == 0  # phase_idx directly from dictionary

def test_species_formatting():
    """Test string formatting of Species with different phases"""
    # Test LaTeX formatting
    co2_g = Species.from_formula('CO2(g)')
    assert co2_g.latex_name == 'CO_{2}(g)'
    
    # Test Unicode formatting
    h2o_l = Species.from_formula('H2O(l)')
    assert h2o_l.unicode_name == 'H₂O(l)'
    
    # Test HTML formatting
    nacl_s = Species.from_formula('NaCl(s)')
    assert nacl_s.html_name == 'NaCl(s)'

if __name__ == '__main__':
    pytest.main(['-v', __file__])
