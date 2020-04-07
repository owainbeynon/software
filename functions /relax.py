from ase.io import read
from ase import Atoms
from ase.calculators.aims import Aims


BEA = read('x.cif')

calc = Aims(xc='pbesol',
            k_grid=(2, 2, 2),
            spin='none',
            relativistic=('atomic_zora', 'scalar'),
            vdw_correction_hirshfeld='true',
            relax_geometry='trm 1E-2',
            relax_unit_cell='full',
            sc_accuracy_forces='1E-3')

BEA.set_calculator(calc)

e_BEA = BEA.get_potential_energy()

print(e_BEA)
