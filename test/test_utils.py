import pytest
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np
from atatutils.utils import (
    str2select_cond, build_linear_problem, get_phase_fractions,
    conc_str2conc_dict
)


@pytest.mark.parametrize('query, want', [
    ('key1=value', [('key1', '=', 'value')]),
    ('key=true', [('key', '=', 1)]),
    ('key=false', [('key', '=', 0)]),
    ('key=false,temp=2', [('key', '=', 0), ('temp', '=', 2)]),
    ('key=false,temp=something', [('key', '=', 0), ('temp', '=', 'something')]),
    ('key>1', [('key', '>', 1)]),
    ('key<1', [('key', '<', 1)]),
    ('key<=1', [('key', '<=', 1)]),
    ('key>=1', [('key', '>=', 1)]),
])
def test_str2select_cond(query, want):
    got = str2select_cond(query)

    assert len(got) == len(want)
    print(got, want)
    for g, w in zip(got, want):
        assert all(x == y for x, y in zip(g, w))


def test_build_linear_problem():
    atoms1 = bulk('Al')
    atoms1.calc = SinglePointCalculator(atoms1, energy=0.5)
    atoms2 = bulk('Mg')
    atoms2.calc = SinglePointCalculator(atoms2, energy=-0.5)
    conc = {'Al': 0.5, 'Mg': 0.5}

    problem = build_linear_problem([atoms1, atoms2], conc)

    c = [0.5, -0.25] # -0.25 because there are two atoms in the Mg cell
    A = np.array([[1.0, 1.0],
                  [1.0, 0.0],
                  [0.0, 1.0]])
    b = np.array([1.0, 0.5, 0.5])
    assert np.allclose(c, problem.cost_vector)
    assert np.allclose(A, problem.A_eq)
    assert np.allclose(b, problem.b_eq)

    fractions = get_phase_fractions(problem)
    assert np.allclose([0.5, 0.5], fractions)

@pytest.mark.parametrize('query, expect', [
    ('Al=0.5,Mg=0.5', {'Al': 0.5, 'Mg': 0.5}),
    ('Al=0.5,Mg=0.2', {'Al': 0.5, 'Mg': 0.2}),
    ('Al=0.5,Mg=0.2,Cu=0.3', {'Al': 0.5, 'Mg': 0.2, 'Cu': 0.3})
])
def test_conc_str2conc_dict(query, expect):
    res = conc_str2conc_dict(query)
    assert len(res) == len(expect)

    for k in expect.keys():
        assert expect[k] == pytest.approx(res[k])