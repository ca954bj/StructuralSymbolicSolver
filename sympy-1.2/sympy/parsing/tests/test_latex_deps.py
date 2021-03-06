from sympy.external import import_module
from sympy.utilities import pytest


antlr4 = import_module("antlr4")

# disable tests if antlr4-python*-runtime is not present
if antlr4:
    disabled = True


def test_no_import():
    from sympy.parsing.latex import parse_latex

    with pytest.raises(ImportError):
        parse_latex('1 + 1')
