from ase.io.trajectory import Trajectory
import numpy as np
from ase.io import read
from ase.geometry.analysis import Analysis
#old file see bonds.py for better version
file= read('tin_acetate.xyz')

traj = Trajectory('tin_acetate.traj', 'w')
traj.write(file)
file = read('tin_acetate.traj')

ana = Analysis(file)

SiOBonds = ana.get_bonds('Sn', 'O')
SiOSiAngles = ana.get_angles('O', 'Sn', 'O')

print("there are {} Si-O bonds in BETA".format(len(SiOBonds[0])))
print("there are {} Si-O-Si angles in BETA".format(len(SiOSiAngles[0])))

SiOBondsValues = ana.get_values(SiOBonds)
SiOSiAngleValues = ana.get_values(SiOSiAngles)

print("bond length data:")
print("the average Si-O bond length is {}.".format(np.average(SiOBondsValues)))
print("the minimum Si-O Distance is:", np.amin(SiOBondsValues))
print("the maximum Si-O Distance is:", np.amax(SiOBondsValues))

print("bond angle data:")
print("the average Si-O-Si angle is {}.".format(np.average(SiOSiAngleValues)))
print("the maximum Si-O-Si angle is:", np.amax(SiOSiAngleValues))
print("the minimum Si-O-Si angle is:", np.amin(SiOSiAngleValues))
