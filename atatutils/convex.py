from typing import Sequence, Set, Dict, List
from collections import Counter
from ase import Atoms
from scipy.spatial import ConvexHull
import numpy as np


def unique_elements(atoms: Sequence[Atoms]) -> Set[str]:
    """
    Return the unique elements in the list of atoms
    """
    unique_elems = set()
    for a in atoms:
        unique_elems = unique_elems.union(set(a.symbols))
    return unique_elems


def normalize_energy_data(atoms: Sequence[Atoms]) -> List[Dict[str, float]]:
    """
    Returns concentrations and normalized energies
    """
    res = []
    for a in atoms:
        conc = Counter(a.symbols)
        new_res = {k: v/len(a) for k, v in conc.items()}
        new_res['energy'] = a.get_potential_energy()/len(a)
        res.append(new_res)
    return res


def convex_hull(atoms: Sequence[Atoms], obs_point: float = -1e3) -> ConvexHull:
    """
    Return the (negative part) convex hull from a sequence of atoms

    :param atoms: Sequence of atoms
    :param obs_point: Observeration point where "good" facets must be seen. See
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
        for further information
    """
    # Skip one element since sum of concentration must be 1
    elems = sorted(list(unique_elements(atoms)))[:-1]

    norm_data = normalize_energy_data(atoms)
    data = np.zeros((len(atoms)+1, len(elems)+1))

    for i, item in enumerate(norm_data):
        data[i, :len(elems)] = [item.get(k, 0.0) for k in elems]
        data[i, -1] = item['energy']

    # Add obseveration point
    data[-1, :len(elems)] = 0.5
    data[-1, -1] = obs_point
    return ConvexHull(data, qhull_options=f'QG{data.shape[0]-1}')