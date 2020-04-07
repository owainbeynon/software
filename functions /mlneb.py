from ase.io import read
from ase import Atoms
from ase.neb import NEB
from ase.optimize import BFGS
import matplotlib.pyplot as plt
from catlearn.optimize.mlneb import MLNEB
from ase.neb import NEBTools
from catlearn.optimize.tools import plotneb
from ase.calculators.aims import Aims
from ase.build import molecule

def my_calc():
    return Aims(xc='pbesol',
                spin='none',
                relativistic=('atomic_zora', 'scalar'),
                vdw_correction_hirshfeld='true',
                k_grid=(2, 2, 2,),
                compute_forces=True, 
                final_forces_cleaned=True)

initial = read('initial.traj')
final = read('final.traj')

n = 7

calculator = my_calc()

neb_catlearn = MLNEB(start=initial,
                     end=final,
                     ase_calc=calculator,
                     n_images=n,
                     interpolation='idpp', restart=False)
neb_catlearn.run(fmax=0.01, trajectory='ML-NEB.traj', full_output=False)

