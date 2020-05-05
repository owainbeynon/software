#!/usr/bin/env python3

import matplotlib.pyplot as plt
from mulliken import parse_mulliken_file
from graphs import get_graph_linetype, get_graph_colour, set_graph_axes_mulliken
from bonds import get_indices_of_elements

# Read in atoms information
output_file = "stage4.traj"
from ase.io import read
atoms = read(output_file)

# Read in Mulliken data from file
mulliken_file = "Mulliken.out"
mulliken_data = parse_mulliken_file(mulliken_file)

#### Assertion statements ####
#assert(mulliken_data.get_natoms() == 6)
#assert(mulliken_data.get_nspin() == 1)
#assert(mulliken_data.get_nkpts() == 184)
#assert(mulliken_data.get_nstates() == 61)
#assert(mulliken_data.get_data_integrity())
#####

# Collect all the density of states data to plot
x, all_data = mulliken_data.get_all_plot_data()

# Collect the indices for each element we are interested in
sn_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'sn')
o_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'o')
c_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'c')
si_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'si')
h_indices = get_indices_of_elements(atoms.get_chemical_symbols(), 'h')

# Collect the density of states data to plot as a function of atomic label
x, sn_p = mulliken_data.get_plot_data(sn_indices, range(mulliken_data.get_nspin()),
                                    range(mulliken_data.get_nkpts()), 'p')
x, o_p = mulliken_data.get_plot_data(o_indices, range(mulliken_data.get_nspin()),
                                    range(mulliken_data.get_nkpts()), 'p')
x, c_p = mulliken_data.get_plot_data(c_indices, range(mulliken_data.get_nspin()),
                                    range(mulliken_data.get_nkpts()), 'p')
x, h_p = mulliken_data.get_plot_data(h_indices, range(mulliken_data.get_nspin()),
                                    range(mulliken_data.get_nkpts()), 'p')
x, si_p = mulliken_data.get_plot_data(si_indices, range(mulliken_data.get_nspin()),
                                    range(mulliken_data.get_nkpts()), 'p')


# Put this at the end so it covers everything else and shows the outline of the DOS correctly
for sp in range(len(all_data)):
    if sp == 0:
        plt.plot(x, sn_p[sp], lw=2, label='Sn 5p', color=get_graph_colour(0), ls=get_graph_linetype())
        plt.plot(x, sn_p[sp]+o_p[sp], lw=2, label='O 2p', color=get_graph_colour(1), ls=get_graph_linetype())
        plt.plot(x, sn_p[sp]+o_p[sp]+c_p[sp], lw=2, label='C 2p', color=get_graph_colour(2), ls=get_graph_linetype())
        plt.plot(x, sn_p[sp]+o_p[sp]+c_p[sp]+h_p[sp],lw=2, label='Si 3p', color=get_graph_colour(3), ls=get_graph_linetype())
        plt.plot(x, sn_p[sp]+o_p[sp]+c_p[sp]+h_p[sp]+si_p[sp], lw=2, label='H 1s', color=get_graph_colour(4), ls=get_graph_linetype())
        plt.plot(x, all_data[sp], lw=2, label='All', color='black', ls=get_graph_linetype())
    else: # (sp == 1)
        plt.plot(x, -(sn_p[sp]), lw=2, color=get_graph_colour(0), ls=get_graph_linetype())
        plt.plot(x, -(sn_p[sp]+o_p[sp]), lw=2, color=get_graph_colour(1), ls=get_graph_linetype())
        plt.plot(x, -(sn_p[sp]+o_p[sp]+c_p[sp]), lw=2, color=get_graph_colour(2), ls=get_graph_linetype())
        plt.plot(x, -(sn_p[sp]+o_p[sp]+c_p[sp]+h_p[sp]), lw=2, color=get_graph_colour(3), ls=get_graph_linetype())
        plt.plot(x, -(sn_p[sp]+o_p[sp]+c_p[sp]+h_p[sp]+si_p[sp]), lw=2, color=get_graph_colour(4), ls=get_graph_linetype())
        plt.plot(x, -(all_data[sp]), lw=2, label='All', color='black,', ls=get_graph_linetype())

# Work to rescale axes. Extracts the maximum y-value
set_graph_axes_mulliken(plt, x, all_data, mulliken_data.get_homo(), mulliken_data.get_graph_xlabel())

# Add a legend
plt.legend()

# Display the graphs
plt.show()