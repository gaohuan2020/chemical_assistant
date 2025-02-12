import pytest
from library.util.parsing import (
    formula_to_composition,
    formula_to_html,
    formula_to_unicode,
    formula_to_latex,
)


def test_simple_molecules():
    """Test parsing of simple molecular formulas"""
    assert formula_to_composition("H2O") == {1: 2, 8: 1}
    assert formula_to_composition("CO2") == {6: 1, 8: 2}
    assert formula_to_composition("NaCl") == {11: 1, 17: 1}


def test_ions():
    """Test parsing of ionic compounds"""
    assert formula_to_composition("NH4+") == {0: 1, 1: 4, 7: 1}
    assert formula_to_composition("SO4-2") == {0: -2, 16: 1, 8: 4}
    assert formula_to_composition("Fe+3") == {0: 3, 26: 1}
    assert formula_to_composition("OH-") == {0: -1, 1: 1, 8: 1}


def test_hydrates():
    """Test parsing of hydrated compounds"""
    assert formula_to_composition("CuSO4..5H2O") == {
        29: 1,  # Cu
        16: 1,  # S
        8: 9,  # O (4 + 5)
        1: 10,  # H (2 * 5)
    }
    assert formula_to_composition("Na2CO3..7H2O") == {
        11: 2,  # Na
        6: 1,  # C
        8: 10,  # O (3 + 7)
        1: 14,  # H (2 * 7)
    }


def test_complex_ions():
    """Test parsing of complex ions"""
    assert formula_to_composition("Fe(CN)6-3") == {
        0: -3,  # charge
        26: 1,  # Fe
        6: 6,  # C
        7: 6,  # N
    }
    assert formula_to_composition("Co(NH3)6+3") == {
        0: 3,  # charge
        27: 1,  # Co
        7: 6,  # N
        1: 18,  # H
    }


def test_state_markers():
    """Test parsing with state markers"""
    assert formula_to_composition("H2O(l)") == {1: 2, 8: 1}
    assert formula_to_composition("NaCl(s)") == {11: 1, 17: 1}
    assert formula_to_composition("CO2(g)") == {6: 1, 8: 2}
    assert formula_to_composition("Na+(aq)") == {0: 1, 11: 1}


def test_greek_prefixes():
    """Test parsing with greek letter prefixes"""
    assert formula_to_composition("alpha-Fe2O3") == {26: 2, 8: 3}
    assert formula_to_composition("beta-HCl") == {1: 1, 17: 1}


def test_decimal_subscripts():
    """Test parsing formulas with decimal subscripts"""
    assert formula_to_composition("UO2.5") == {92: 1, 8: 2.5}
    assert formula_to_composition("Fe0.98O") == {26: 0.98, 8: 1}


def test_error_cases():
    """Test error handling"""
    with pytest.raises(ValueError):
        formula_to_composition("H2O++")  # Multiple charge tokens

    with pytest.raises(ValueError):
        formula_to_composition("Fe/3+")  # Deprecated charge notation

    with pytest.raises(ValueError):
        formula_to_composition("H2O+2-")  # Both + and - present


def test_edge_cases():
    """Test edge cases"""
    assert formula_to_composition("e") == {}  # electron
    assert formula_to_composition("H") == {1: 1}  # single atom
    assert formula_to_composition("H2") == {1: 2}  # diatomic molecule


def test_complex_cases():
    """Test more complex formulas"""
    assert formula_to_composition("[Fe(H2O)6]Cl3") == {
        26: 1,  # Fe
        1: 12,  # H
        8: 6,  # O
        17: 3,  # Cl
    }

    assert formula_to_composition("K4[Fe(CN)6]..3H2O") == {
        19: 4,  # K
        26: 1,  # Fe
        6: 6,  # C
        7: 6,  # N
        1: 6,  # H (from water)
        8: 3,  # O (from water)
    }


def test_formula_to_latex():
    """Test conversion of chemical formulas to LaTeX format"""
    # Test basic molecules
    assert formula_to_latex("H2O") == "H_{2}O"
    assert formula_to_latex("NH4+") == "NH_{4}^{+}"
    assert formula_to_latex("Fe(CN)6+2") == "Fe(CN)_{6}^{2+}"

    # Test with state markers
    assert formula_to_latex("Fe(CN)6+2(aq)") == "Fe(CN)_{6}^{2+}(aq)"

    # Test with prefixes
    assert formula_to_latex(".NHO-(aq)") == "^\\bullet NHO^{-}(aq)"
    assert formula_to_latex("alpha-FeOOH(s)") == "\\alpha-FeOOH(s)"

    # Test complex cases
    assert formula_to_latex("Na2CO3..7H2O") == "Na_{2}CO_{3}\\cdot 7H_{2}O"
    assert formula_to_latex("[Fe(H2O)6]Cl3") == "[Fe(H_{2}O)_{6}]Cl_{3}"
    assert formula_to_latex("K4[Fe(CN)6]..3H2O") == "K_{4}[Fe(CN)_{6}]\\cdot 3H_{2}O"


def test_formula_to_unicode():
    """Test conversion of chemical formulas to Unicode format"""
    # Test basic molecules
    assert formula_to_unicode("H2O") == "H₂O"
    assert formula_to_unicode("NH4+") == "NH₄⁺"
    assert formula_to_unicode("Fe(CN)6+2") == "Fe(CN)₆²⁺"

    # Test with state markers
    assert formula_to_unicode("Fe(CN)6+2(aq)") == "Fe(CN)₆²⁺(aq)"

    # Test with prefixes
    assert formula_to_unicode(".NHO-(aq)") == "⋅NHO⁻(aq)"
    assert formula_to_unicode("alpha-FeOOH(s)") == "α-FeOOH(s)"

    # Test hydrates and complex compounds
    assert formula_to_unicode("Na2CO3..7H2O") == "Na₂CO₃·7H₂O"
    assert formula_to_unicode("[Fe(H2O)6]Cl3") == "[Fe(H₂O)₆]Cl₃"
    assert formula_to_unicode("K4[Fe(CN)6]..3H2O") == "K₄[Fe(CN)₆]·3H₂O"


def test_formula_to_html():
    """Test conversion of chemical formulas to HTML format"""
    # Test basic molecules
    assert formula_to_html("H2O") == "H<sub>2</sub>O"
    assert formula_to_html("NH4+") == "NH<sub>4</sub><sup>+</sup>"
    assert formula_to_html("Fe(CN)6+2") == "Fe(CN)<sub>6</sub><sup>2+</sup>"

    # Test with state markers
    assert formula_to_html("Fe(CN)6+2(aq)") == "Fe(CN)<sub>6</sub><sup>2+</sup>(aq)"

    # Test with prefixes
    assert formula_to_html(".NHO-(aq)") == "&sdot;NHO<sup>-</sup>(aq)"
    assert formula_to_html("alpha-FeOOH(s)") == "&alpha;-FeOOH(s)"

    # Test hydrates and complex compounds
    assert (
        formula_to_html("Na2CO3..7H2O")
        == "Na<sub>2</sub>CO<sub>3</sub>&sdot;7H<sub>2</sub>O"
    )
    assert (
        formula_to_html("[Fe(H2O)6]Cl3")
        == "[Fe(H<sub>2</sub>O)<sub>6</sub>]Cl<sub>3</sub>"
    )
    assert (
        formula_to_html("K4[Fe(CN)6]..3H2O")
        == "K<sub>4</sub>[Fe(CN)<sub>6</sub>]&sdot;3H<sub>2</sub>O"
    )


if __name__ == "__main__":
    pytest.main(["-v", __file__])
