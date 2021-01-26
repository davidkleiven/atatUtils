#!/usr/bin/env python
from pymatgen.io.atat import Mcsqs
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.trajectory import TrajectoryReader
import argparse
import os

def convert(traj, outfolder, start=0):
    for i, atoms in enumerate(traj):
        folder = outfolder + f'/{start+i}'
        os.mkdir(folder)
        structure = AseAtomsAdaptor.get_structure(atoms)
        msc = Mcsqs(structure)
        with open(folder + '/str.in', 'w') as outfile:
            outfile.write(msc.to_string())
        with open(folder + '/energy', 'w') as outfile:
            outfile.write(f"{atoms.get_potential_energy()}\n")
        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('traj', help="Trajectory file")
    parser.add_argument('outfolder', help="Folder with outputs")
    parser.add_argument('--start', type=int, default=0, help="Number where ATAT enumeration starts")
    args = parser.parse_args()

    traj = TrajectoryReader(args.traj)
    convert(traj, args.outfolder, args.start)


if __name__ == '__main__':
    main()