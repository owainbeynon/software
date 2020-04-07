from ase import Atoms
from ase.io import read
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.calculators.aims import Aims
from ase.vibrations import Vibrations

​BEAH = read('BEAH.cif')

calc = Aims(xc='pbesol',
            relativistic='atomic_zora', 'scalar',
            k_grid=(2,2,2),
            vdw_correction='hirshfeld',
            sc_accuracy_etot=1e-4,
            sc_accuracy_eev=1e-4,
            sc_accuracy_rho=1e-4,
            sc_accuracy_forces=1e-4))


BEAH.set_calculator(calc)

from ase.constraints import FixAtoms

c = FixAtoms(indices=[i for i in range(0,192)])
BEAH.set_constraint(c)

e0 = BEAH.get_potential_energy()

QuasiNewton(BEAH).run(fmax=0.05)

​from ase.vibrations import Vibrations

vib = Vibrations(BEAH, indices=[193])
vib.run()
vib.summary()
vib.write_mode()
