from collections import Counter
from typing import List, Tuple, Union, NamedTuple, Dict, Sequence
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np
from scipy.optimize import linprog


ASE_SELECT_COND = Tuple[str, str, Union[str, int]]

COMPARATORS = ['>=', '<=', '>', '<', '=']

def find_comparator(query: str) -> str:
    for c in COMPARATORS:
        if c in query:
            return c
    raise ValueError("Did not find any of {COMPARATORS} in the query")

def concentration(atoms: Atoms) -> Dict[str, float]:
    cnt = Counter(atoms.symbols)
    return {k: v/len(atoms) for k, v in cnt.items()}


def str2select_cond(query: str) -> List[ASE_SELECT_COND]:
    select_conds = []
    for item in query.split(','):
        comparator = find_comparator(query)
        key, value = item.split(comparator)
        
        if value.lower() == 'true':
            value = 1
        elif value.lower() == 'false':
            value = 0
        
        try:
            value = int(value)
        except:
            pass

        select_conds.append((key, comparator, value))
    return select_conds


class PhaseFractionProblem(NamedTuple):
    cost_vector: np.array
    A_eq: np.array
    b_eq: np.array


def build_linear_problem(atoms: Sequence[Atoms], conc: Dict[str, float]) -> PhaseFractionProblem:
    energies = [a.get_potential_energy()/len(a) for a in atoms]
    concs = [concentration(a) for a in atoms]
    A_eq = np.zeros((len(conc)+1, len(energies)))
    b_eq = np.zeros(len(conc)+1)
    
    # Phase fractions must sum to one
    A_eq[0, :] = 1.0
    b_eq[0] = 1.0

    # Concentrations must match
    all_species = sorted(list(conc.keys()))  # Just make sure species always come in the same order
    for i, species in enumerate(all_species):
        for j, c in enumerate(concs):
            A_eq[i+1, j] = c.get(species, 0.0)
        b_eq[i+1] = conc.get(species, 0.0)

    assert A_eq.shape[0] == len(b_eq)
    assert A_eq.shape[1] == len(energies)

    return PhaseFractionProblem(cost_vector=energies, A_eq=A_eq, b_eq=b_eq)

def get_phase_fractions(problem: PhaseFractionProblem) -> np.ndarray:
    res = linprog(problem.cost_vector, A_eq=problem.A_eq, b_eq=problem.b_eq)
    return res.x

DB_IDs = List[int]

def db_iterator(db, scond):
    if scond == "":
        yield from db.select()
    else:
        condition = str2select_cond(scond)
        yield from db.select(condition)

def atoms_with_singlepoint_calc(db, query) -> Tuple[List[Atoms], DB_IDs]:
    atoms = []
    db_ids = []
    for row in db_iterator(db, query):
        a = row.toatoms()
        E = row.get('energy', None)
        if E is not None:
            a.calc = SinglePointCalculator(a, energy=E)
            atoms.append(a)
            db_ids.append(row.id)
    return atoms, db_ids


def conc_str2conc_dict(conc_str: str) -> Dict[str, float]:
    """
    Converts a concentration string into an equivalent dictionary

    Example:
    Al=0.4,Mg=0.3,Cu=0.3 is converted to{"Al": 0.4, "Mg": 0.3, "Cu": 0.3}

    :param conc_str: String with the concentration
    """
    res = {}
    for item in conc_str.split(','):
        key, value = item.split('=')
        res[key] = float(value)
    return res
