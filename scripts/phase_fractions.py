import argparse
from atatutils.utils import atoms_with_singlepoint_calc
from atatutils.utils import conc_str2conc_dict, build_linear_problem, get_phase_fractions
from ase.db import connect
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Program for determining the phase fraction at zero kelvin from an ASE DB")
    parser.add_argument("db", type=str, help="Name of ASE database")
    parser.add_argument("--sel", default="", type=str, help="Select condition (comma separated "
                                                            "in the form key1=val1,key2=val2)")
    parser.add_argument("--conc", default="", type=str,
                        help="Comma separated list with concentrations. Must sum to one. Example: "
                             "Al=0.4,Mg=0.3,Cu=0.3")
    parser.add_argument("--minpf", type=float, default=1e-6,
                        help="Phases with lower phase fraction than the passed argument are "
                             "not listed in the output")

    args = parser.parse_args()
    db = connect(args.db)
    atoms, db_ids = atoms_with_singlepoint_calc(db, args.sel)

    conc = conc_str2conc_dict(args.conc)

    if abs(sum(conc.values()) - 1.0) > 1e-8:
        print("Concentrations do not sum to 1!")
        print(f"Extracted concenratations: {conc}")
        return

    print(f"Finding phase fractions for concentration {conc}")
    print("Building linear problem...")
    problem = build_linear_problem(atoms, conc)

    num_unknown = len(problem.cost_vector)
    num_constraints = problem.A_eq.shape[0]
    print(f"Build problem with: {num_unknown} unknowns and {num_constraints} equality constraints")

    print("Solving linear problem...")
    res = get_phase_fractions(problem)

    srt_idx = np.argsort(res)[::-1]
    N = 65
    print('-'*N)
    print('DB ID | Chem. formula        | Phase fraction | Energy (eV/atom)')
    print('-'*N)
    for idx in srt_idx:
        chemical_formula = atoms[idx].get_chemical_formula()
        pf = res[idx]
        E = atoms[idx].get_potential_energy()/len(atoms[idx])
        
        if pf < args.minpf:
            break
        print(f'{db_ids[idx]:5d} | {chemical_formula:<20} | {pf:14.4f} | {E:16.4f}')
    print('-'*N)


if __name__ == '__main__':
    main()
