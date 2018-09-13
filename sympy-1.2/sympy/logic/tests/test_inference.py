"""For more tests on satisfiability, see test_dimacs"""

from sympy import symbols, Q
from sympy.core.compatibility import range
from sympy.logic.boolalg import And, Implies, Equivalent, true, false
from sympy.logic.inference import literal_symbol, \
     pl_true, satisfiable, valid, entails, PropKB
from sympy.logic.algorithms.dpll import dpll, dpll_satisfiable, \
    find_pure_symbol, find_unit_clause, unit_propagate, \
    find_pure_symbol_int_repr, find_unit_clause_int_repr, \
    unit_propagate_int_repr
from sympy.logic.algorithms.dpll2 import dpll_satisfiable as dpll2_satisfiable
from sympy.utilities.pytest import raises


def test_literal():
    A, B = symbols('A,B')
    assert literal_symbol(True) is True
    assert literal_symbol(False) is False
    assert literal_symbol(A) is A
    assert literal_symbol(~A) is A


def test_find_pure_symbol():
    A, B, C = symbols('A,B,C')
    assert find_pure_symbol([A], [A]) == (A, True)
    assert find_pure_symbol([A, B], [~A | B, ~B | A]) == (None, None)
    assert find_pure_symbol([A, B, C], [ A | ~B, ~B | ~C, C | A]) == (A, True)
    assert find_pure_symbol([A, B, C], [~A | B, B | ~C, C | A]) == (B, True)
    assert find_pure_symbol([A, B, C], [~A | ~B, ~B | ~C, C | A]) == (B, False)
    assert find_pure_symbol(
        [A, B, C], [~A | B, ~B | ~C, C | A]) == (None, None)


def test_find_pure_symbol_int_repr():
    assert find_pure_symbol_int_repr([1], [set([1])]) == (1, True)
    assert find_pure_symbol_int_repr([1, 2],
                [set([-1, 2]), set([-2, 1])]) == (None, None)
    assert find_pure_symbol_int_repr([1, 2, 3],
                [set([1, -2]), set([-2, -3]), set([3, 1])]) == (1, True)
    assert find_pure_symbol_int_repr([1, 2, 3],
                [set([-1, 2]), set([2, -3]), set([3, 1])]) == (2, True)
    assert find_pure_symbol_int_repr([1, 2, 3],
                [set([-1, -2]), set([-2, -3]), set([3, 1])]) == (2, False)
    assert find_pure_symbol_int_repr([1, 2, 3],
                [set([-1, 2]), set([-2, -3]), set([3, 1])]) == (None, None)


def test_unit_clause():
    A, B, C = symbols('A,B,C')
    assert find_unit_clause([A], {}) == (A, True)
    assert find_unit_clause([A, ~A], {}) == (A, True)  # Wrong ??
    assert find_unit_clause([A | B], {A: True}) == (B, True)
    assert find_unit_clause([A | B], {B: True}) == (A, True)
    assert find_unit_clause(
        [A | B | C, B | ~C, A | ~B], {A: True}) == (B, False)
    assert find_unit_clause([A | B | C, B | ~C, A | B], {A: True}) == (B, True)
    assert find_unit_clause([A | B | C, B | ~C, A ], {}) == (A, True)


def test_unit_clause_int_repr():
    assert find_unit_clause_int_repr(map(set, [[1]]), {}) == (1, True)
    assert find_unit_clause_int_repr(map(set, [[1], [-1]]), {}) == (1, True)
    assert find_unit_clause_int_repr([set([1, 2])], {1: True}) == (2, True)
    assert find_unit_clause_int_repr([set([1, 2])], {2: True}) == (1, True)
    assert find_unit_clause_int_repr(map(set,
        [[1, 2, 3], [2, -3], [1, -2]]), {1: True}) == (2, False)
    assert find_unit_clause_int_repr(map(set,
        [[1, 2, 3], [3, -3], [1, 2]]), {1: True}) == (2, True)

    A, B, C = symbols('A,B,C')
    assert find_unit_clause([A | B | C, B | ~C, A ], {}) == (A, True)


def test_unit_propagate():
    A, B, C = symbols('A,B,C')
    assert unit_propagate([A | B], A) == []
    assert unit_propagate([A | B, ~A | C, ~C | B, A], A) == [C, ~C | B, A]


def test_unit_propagate_int_repr():
    assert unit_propagate_int_repr([set([1, 2])], 1) == []
    assert unit_propagate_int_repr(map(set,
        [[1, 2], [-1, 3], [-3, 2], [1]]), 1) == [set([3]), set([-3, 2])]


def test_dpll():
    """This is also tested in test_dimacs"""
    A, B, C = symbols('A,B,C')
    assert dpll([A | B], [A, B], {A: True, B: True}) == {A: True, B: True}


def test_dpll_satisfiable():
    A, B, C = symbols('A,B,C')
    assert dpll_satisfiable( A & ~A ) is False
    assert dpll_satisfiable( A & ~B ) == {A: True, B: False}
    assert dpll_satisfiable(
        A | B ) in ({A: True}, {B: True}, {A: True, B: True})
    assert dpll_satisfiable(
        (~A | B) & (~B | A) ) in ({A: True, B: True}, {A: False, B: False})
    assert dpll_satisfiable( (A | B) & (~B | C) ) in ({A: True, B: False},
            {A: True, C: True}, {B: True, C: True})
    assert dpll_satisfiable( A & B & C  ) == {A: True, B: True, C: True}
    assert dpll_satisfiable( (A | B) & (A >> B) ) == {B: True}
    assert dpll_satisfiable( Equivalent(A, B) & A ) == {A: True, B: True}
    assert dpll_satisfiable( Equivalent(A, B) & ~A ) == {A: False, B: False}


def test_dpll2_satisfiable():
    A, B, C = symbols('A,B,C')
    assert dpll2_satisfiable( A & ~A ) is False
    assert dpll2_satisfiable( A & ~B ) == {A: True, B: False}
    assert dpll2_satisfiable(
        A | B ) in ({A: True}, {B: True}, {A: True, B: True})
    assert dpll2_satisfiable(
        (~A | B) & (~B | A) ) in ({A: True, B: True}, {A: False, B: False})
    assert dpll2_satisfiable( (A | B) & (~B | C) ) in ({A: True, B: False, C: True},
        {A: True, B: True, C: True})
    assert dpll2_satisfiable( A & B & C  ) == {A: True, B: True, C: True}
    assert dpll2_satisfiable( (A | B) & (A >> B) ) in ({B: True, A: False},
        {B: True, A: True})
    assert dpll2_satisfiable( Equivalent(A, B) & A ) == {A: True, B: True}
    assert dpll2_satisfiable( Equivalent(A, B) & ~A ) == {A: False, B: False}


def test_satisfiable():
    A, B, C = symbols('A,B,C')
    assert satisfiable(A & (A >> B) & ~B) is False


def test_valid():
    A, B, C = symbols('A,B,C')
    assert valid(A >> (B >> A)) is True
    assert valid((A >> (B >> C)) >> ((A >> B) >> (A >> C))) is True
    assert valid((~B >> ~A) >> (A >> B)) is True
    assert valid(A | B | C) is False
    assert valid(A >> B) is False


def test_pl_true():
    A, B, C = symbols('A,B,C')
    assert pl_true(True) is True
    assert pl_true( A & B, {A: True, B: True}) is True
    assert pl_true( A | B, {A: True}) is True
    assert pl_true( A | B, {B: True}) is True
    assert pl_true( A | B, {A: None, B: True}) is True
    assert pl_true( A >> B, {A: False}) is True
    assert pl_true( A | B | ~C, {A: False, B: True, C: True}) is True
    assert pl_true(Equivalent(A, B), {A: False, B: False}) is True

    # test for false
    assert pl_true(False) is False
    assert pl_true( A & B, {A: False, B: False}) is False
    assert pl_true( A & B, {A: False}) is False
    assert pl_true( A & B, {B: False}) is False
    assert pl_true( A | B, {A: False, B: False}) is False

    #test for None
    assert pl_true(B, {B: None}) is None
    assert pl_true( A & B, {A: True, B: None}) is None
    assert pl_true( A >> B, {A: True, B: None}) is None
    assert pl_true(Equivalent(A, B), {A: None}) is None
    assert pl_true(Equivalent(A, B), {A: True, B: None}) is None

    # Test for deep
    assert pl_true(A | B, {A: False}, deep=True) is None
    assert pl_true(~A & ~B, {A: False}, deep=True) is None
    assert pl_true(A | B, {A: False, B: False}, deep=True) is False
    assert pl_true(A & B & (~A | ~B), {A: True}, deep=True) is False
    assert pl_true((C >> A) >> (B >> A), {C: True}, deep=True) is True


def test_pl_true_wrong_input():
    from sympy import pi
    raises(ValueError, lambda: pl_true('John Cleese'))
    raises(ValueError, lambda: pl_true(42 + pi + pi ** 2))
    raises(ValueError, lambda: pl_true(42))


def test_entails():
    A, B, C = symbols('A, B, C')
    assert entails(A, [A >> B, ~B]) is False
    assert entails(B, [Equivalent(A, B), A]) is True
    assert entails((A >> B) >> (~A >> ~B)) is False
    assert entails((A >> B) >> (~B >> ~A)) is True


def test_PropKB():
    A, B, C = symbols('A,B,C')
    kb = PropKB()
    assert kb.ask(A >> B) is False
    assert kb.ask(A >> (B >> A)) is True
    kb.tell(A >> B)
    kb.tell(B >> C)
    assert kb.ask(A) is False
    assert kb.ask(B) is False
    assert kb.ask(C) is False
    assert kb.ask(~A) is False
    assert kb.ask(~B) is False
    assert kb.ask(~C) is False
    assert kb.ask(A >> C) is True
    kb.tell(A)
    assert kb.ask(A) is True
    assert kb.ask(B) is True
    assert kb.ask(C) is True
    assert kb.ask(~C) is False
    kb.retract(A)
    assert kb.ask(C) is False


def test_propKB_tolerant():
    """"tolerant to bad input"""
    kb = PropKB()
    A, B, C = symbols('A,B,C')
    assert kb.ask(B) is False

def test_satisfiable_non_symbols():
    x, y = symbols('x y')
    assumptions = Q.zero(x*y)
    facts = Implies(Q.zero(x*y), Q.zero(x) | Q.zero(y))
    query = ~Q.zero(x) & ~Q.zero(y)
    refutations = [
        {Q.zero(x): True, Q.zero(x*y): True},
        {Q.zero(y): True, Q.zero(x*y): True},
        {Q.zero(x): True, Q.zero(y): True, Q.zero(x*y): True},
        {Q.zero(x): True, Q.zero(y): False, Q.zero(x*y): True},
        {Q.zero(x): False, Q.zero(y): True, Q.zero(x*y): True}]
    assert not satisfiable(And(assumptions, facts, query), algorithm='dpll')
    assert satisfiable(And(assumptions, facts, ~query), algorithm='dpll') in refutations
    assert not satisfiable(And(assumptions, facts, query), algorithm='dpll2')
    assert satisfiable(And(assumptions, facts, ~query), algorithm='dpll2') in refutations

def test_satisfiable_bool():
    from sympy.core.singleton import S
    assert satisfiable(true) == {true: true}
    assert satisfiable(S.true) == {true: true}
    assert satisfiable(false) is False
    assert satisfiable(S.false) is False


def test_satisfiable_all_models():
    from sympy.abc import A, B
    assert next(satisfiable(False, all_models=True)) is False
    assert list(satisfiable((A >> ~A) & A , all_models=True)) == [False]
    assert list(satisfiable(True, all_models=True)) == [{true: true}]

    models = [{A: True, B: False}, {A: False, B: True}]
    result = satisfiable(A ^ B, all_models=True)
    models.remove(next(result))
    models.remove(next(result))
    raises(StopIteration, lambda: next(result))
    assert not models

    assert list(satisfiable(Equivalent(A, B), all_models=True)) == \
    [{A: False, B: False}, {A: True, B: True}]

    models = [{A: False, B: False}, {A: False, B: True}, {A: True, B: True}]
    for model in satisfiable(A >> B, all_models=True):
        models.remove(model)
    assert not models

    # This is a santiy test to check that only the required number
    # of solutions are generated. The expr below has 2**100 - 1 models
    # which would time out the test if all are generated at once.
    from sympy import numbered_symbols
    from sympy.logic.boolalg import Or
    sym = numbered_symbols()
    X = [next(sym) for i in range(100)]
    result = satisfiable(Or(*X), all_models=True)
    for i in range(10):
        assert next(result)
