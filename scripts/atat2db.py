import argparse
import os
from ase.db import connect
import subprocess
from ase.io import read

ATAT_GENERATED = 'atat_generated'
SCRIPT_NAME = '.str2cif_wrapper.sh'

def str2cif_sh_script():
    return "str2cif < $1 > $2\n"

def transfer_to_db(folder, db, structure_file):
    print("Transferring structures from ATAT folder system to DB")
    num_transferred = 0
    with open(SCRIPT_NAME, 'w') as out:
        out.write(str2cif_sh_script())

    for root, _, files in os.walk(folder):
        if structure_file in files:
            # Check if the folder already exists in the DB
            exists = True
            try:
                db.get(folder=root, struct_type=ATAT_GENERATED)
            except KeyError:
                exists = False

            if not exists:
                str_file = root + f"/{structure_file}"
                cif_file = root + "/structure.cif"

                # Convert to a CIF file
                subprocess.run(['sh', SCRIPT_NAME, str_file, cif_file])

                atoms = read(cif_file)
                db.write(atoms, folder=root, struct_type=ATAT_GENERATED)
                num_transferred += 1
    print(f"Transferred {num_transferred} structures to the database")
    os.remove(SCRIPT_NAME)


def parse_prop(prop):
    if ':' not in prop:
        return prop, prop
    splitted = prop.split(':')
    assert len(splitted) == 2
    return splitted

def transfer_prop_from_db(db, prop):
    num_prop_transferred = 0
    db_prop, atat_prop = parse_prop(prop)

    print(f"Writing DB prop {db_prop} to atat prop {atat_prop}")
    for row in db.select(struct_type=ATAT_GENERATED):
        folder= row.folder
        value = row.get(db_prop, None)
        if value is None:
            continue
        
        out_file = folder + f"/{atat_prop}"
        with open(out_file, 'w') as out:
            out.write(f"{value}\n")
        num_prop_transferred += 1
    print(f"Transferred {prop} for {num_prop_transferred} items")


def main():
    parser = argparse.ArgumentParser(description="Transfer an ATAT project to an ASE DB")
    parser.add_argument("folder", type=str, help="ATAT project folder")
    parser.add_argument("db", type=str, help="ASE database")
    parser.add_argument("--prop", type=str, default='energy',
                        help="Property to be exported. If the name expected by ATAT does not "
                             "match the key in the database, a mapping from database to ATAT "
                             "can be specified via <db item>:<atat name>. Thus, if the quantity "
                             "that ATAT want (say energy) is stored in  the DB as for example "
                             "relaxed_energy, specify relaxed_energy:energy")
    parser.add_argument("--str", type=str, default="str.out", help="Name of the structure files")

    args = parser.parse_args()
    db = connect(args.db)
    transfer_to_db(args.folder, db, args.str)
    transfer_prop_from_db(db, args.prop)
   

if __name__ == '__main__':
    main()
