Monte Carlo DPD
===============

This C program uses Monte Carlo simulations and dissipative particle dynamics (DPD) to investigate quantities like energy and pressure of a system of particles. Future plans are to add a polymer, wall and nanopore.

Installation
------------

To get latest version:

	git clone https://github.com/tdunn19/dpd

Usage
-----

To compile and run:

	cd dpd/src
	make opt3
	./dpd.run

Edit src/dpd.inp to change parameters.

Changelog
---------
Version 2.0 (May 24, 2013)
*	the polymer has been added to the system
*   note: brute force method currently does not work with polymer
* 	new variables: n_mon and a_ij broken up into a_mm, a_ms and a_ss
*   calc_energy_list has been split into calc_energy_dpd and calc_energy_mon
*   random_move has been split into random_move_dpd and random_move_mon
*   energy_ij has been renamed to energy_c for conservative
*   new function energy_fene: takes in two bonded monomers, outputs FENE bond potential
    *   r_max: max distance between bonded monomers before the bond breaks
    *   r_eq: equilibrium bond length
    *   r_0: defined as r_max - r_eq
    *   bond_break: system parameter assigned the value of 1 (true) if a monomer-monomer FENE bond breaks (r_ij > r_max)

Version 1.1 (May 23, 2013)
*	some efficiency changes:
*	a new cell list will only be generated when a particles moves into a new cell
*	new energies will only be calculated for particles in neighboring (2 neighbors deep) cells of the randomly chosen particle

Version 1.0 (May 20, 2013)
*	employed cell list method of calculation

