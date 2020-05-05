from neb import switch_indices as Switch
from neb import check_interpolation
from ase.io import read
from ase.visualize import view
index1 = 202
index2 = 202

initial = read("initial.traj")
final = read("final.traj")
changed = Switch('final.traj', index1, index2)

check_interpolation('initial.traj','final.traj', 11)
#check_interpolation('initial.traj', changed, 11, interpolation="idpp")


#interpolation = read("interpolation.traj")
#view(interpolation)
