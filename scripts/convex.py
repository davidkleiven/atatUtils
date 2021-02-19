import argparse
from atatutils.utils import str2select_cond
from atatutils.convex import convex_hull
from ase.db import connect
from ase.calculators.singlepoint import SinglePointCalculator

def db_iterator(db, scond):
    if scond == "":
        yield from db.select()
    else:
        condition = str2select_cond(scond)
        yield from db.select(condition)

def points_on_convex(hull):
    idx = set()
    for s in hull.simplices[hull.good]:
        idx = idx.union(s)
    return idx


def main():
    parser = argparse.ArgumentParser(description="Program for analysing the convex hull from ASE DBs")
    parser.add_argument("db", type=str, help="Name of ASE database")
    parser.add_argument("--sel", default="", type=str, help="Select condition (comma separated "
                                                            "in the form key1=val1,key2=val2)")
    parser.add_argument("--obs", default=-1e3, type=float,
                        help="Observation point for evaluating the convex hull. See scipy convex hull"
                             "describing QG{N} section ")
    
    args = parser.parse_args()

    db = connect(args.db)
    atoms = []
    db_ids = []
    for row in db_iterator(db, args.sel):
        a = row.toatoms()
        E = row.get('energy', None)
        if E is not None:
            a.calc = SinglePointCalculator(a, energy=E)
            atoms.append(a)
            db_ids.append(row.id)
    hull = convex_hull(atoms)

    # Print information of structure on convex hull
    N = 58
    print('-'*N)
    print('DB ID | Chem. formula        | Energy (eV/atom)')
    print('-'*N)
    for idx in points_on_convex(hull):
        a = atoms[idx]
        chemical_formula = a.get_chemical_formula()
        E = a.get_potential_energy()/len(a)
        print(f'{db_ids[idx]:5d} | {chemical_formula:<20} | {E:27.4f}')
    print('-'*N)


if __name__ == '__main__':
    main()