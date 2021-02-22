import pytest
from ase.build import bulk
from atatutils.convex import (
    unique_elements, normalize_energy_data, convex_hull
)
from ase.calculators.singlepoint import SinglePointCalculator


def test_unique_elements():
    atoms = [
        bulk('MgCu', crystalstructure='rocksalt', a=4.0),
        bulk('AlMg', crystalstructure='rocksalt', a=4.0),
        bulk('Cu', crystalstructure='fcc', a=4.0),
        bulk('Al', crystalstructure='fcc', a=4.0)
    ]

    unique = unique_elements(atoms)
    want = set(['Mg', 'Al', 'Cu'])
    print(unique)
    assert unique == want


def test_normalize_energy_data():
    atoms1 = bulk('MgCu', crystalstructure='rocksalt', a=4.0)
    atoms1.calc  = SinglePointCalculator(atoms1, energy=2.0)

    atoms2 = bulk('MgCu', crystalstructure='rocksalt', a=4.0)
    atoms2.calc  = SinglePointCalculator(atoms2, energy=1.0)

    res = normalize_energy_data([atoms1, atoms2])
    want = [
        {
            'Mg': 0.5,
            'Cu': 0.5,
            'energy': 1.0
        },
        {
            'Mg': 0.5,
            'Cu': 0.5,
            'energy': 0.5
        }
    ]

    assert len(res) == len(want)
    for got, w in zip(res, want):
        for key in ['Mg', 'Cu', 'energy']:
            assert got[key] == pytest.approx(w[key])


def test_convex_hull():
    atoms1 = bulk('Al')
    atoms1.calc = SinglePointCalculator(atoms1, energy=0.1)

    atoms2 = bulk('Al')
    atoms2.calc = SinglePointCalculator(atoms2, energy=1.0)

    atoms3 = bulk('AlMg', crystalstructure='rocksalt', a=4.0)
    atoms3.calc = SinglePointCalculator(atoms3, energy=-0.5)

    atoms4 = bulk('AlMg', crystalstructure='rocksalt', a=4.0)
    atoms4.calc = SinglePointCalculator(atoms4, energy=-0.3)

    atoms5 = bulk('Mg')
    atoms5.calc = SinglePointCalculator(atoms5, energy=0.1)

    atoms6 = bulk('Mg')
    atoms6.calc = SinglePointCalculator(atoms6, energy=1.0)

    hull = convex_hull([atoms1, atoms2, atoms3, atoms4, atoms5, atoms6])

    # We expect all even indices to be on the convex hull
    for simplex in hull.simplices[hull.good]:
        for v in simplex:
            assert v % 2 == 0
