import argparse
from atatutils.utils import atoms_with_singlepoint_calc
from atatutils.convex import convex_hull
from ase.db import connect

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
    atoms, db_ids = atoms_with_singlepoint_calc(db, args.sel)
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