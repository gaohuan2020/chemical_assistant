import pytest
from library.chemistry import Substance

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

if __name__ == '__main__':
    pytest.main(['-v', __file__])
