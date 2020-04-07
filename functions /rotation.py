from ase import Atoms
from ase.io import Trajectory
from ase.visualize import view
from ase.io import read, write
from ase.io.trajectory import TrajectoryWriter
import numpy as np
from ase.calculators.emt import EMT


atoms = read('tin_acetate.xyz')

positions = atoms.get_positions()
np.array(positions).tolist()
dihedral = atoms.get_dihedral(9, 7, 2, 14)
print(dihedral)

with open('output.txt', 'w') as fin:
    for theta in range(int(dihedral), int(dihedral + 360), 45):
        atoms.set_dihedral(9, 7, 2, 14, theta, mask=[1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        atoms.set_calculator(EMT())
        energy = atoms.get_potential_energy()
        p = print('Energy is:', theta, energy)
        to_string1 = str(theta)
        to_string2 = str(energy)
        fin.write(to_string1 +'\n')
        fin.write(to_string2 +'\n')
        traj = Trajectory('output'+to_string1+'.traj', 'w', atoms)
        traj.write()
