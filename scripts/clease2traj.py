from ase.io.trajectory import TrajectoryWriter
from ase.db import connect
from ase.calculators.singlepoint import SinglePointCalculator
import argparse

def db2traj(db_name, traj_file):
    traj = TrajectoryWriter(traj_file)
    db = connect(db_name)

    counter = 0
    for row in db.select([('converged', '=', 1)]):
        fid = row.get('final_struct_id', None)
        if fid is None:
            continue
        
        energy = db.get(id=fid).get('energy', None)
        if energy is None:
            continue
            
        atoms = row.toatoms()
        atoms.calc = SinglePointCalculator(atoms, energy=energy)
        traj.write(atoms)
        counter += 1
    print(f"Wrote {counter} structures to {traj_file}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('db', help="CLEASE DB")
    parser.add_argument('traj', help="Trajectory file")
    args = parser.parse_args()
    db2traj(args.db, args.traj)

if __name__ == '__main__':
    main()