Monte Carlo DPD
===============

This C program uses the Metropolis Monte Carlo method and dissipative particle dynamics (DPD) to investigate quantities like energy and pressure of a system of particles. Future plans are to add a polymer, wall and nanopore.

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

Edit src/dpd.inp to change parameters. The following is a description of all the input parameters as well as typical values in brackets:

    n_mon - the number of monomers
    density_s - the solvent density (3)
    calc_list - 1 to use cell list method, 0 to use brute force method of calculation

    length_x - system box length in the x (10)
    length_y - in the y (10)
    length_z - in the z (long enough to fit the polymer)

    density_w - density of wall particles (3)
    n_layers - number of wall layers (2)

    r_c - cutoff radius of particle interaction (1.0)
    dr_max_dpd - max displacement of a dpd monte carlo move (0.2)
    dr_max_mon - max displacement of a monomer monte carlo move (0.1)

    a_mm - monomer-monomer interaction strength (25.0)
    a_ms - monomer-solvent interaction strength (15.0)
    a_ss - solvent-solvent interaction strength (25.0)
    a_sw - solvent-wall interaction strength (9.01)

    pol_init_z - initial z coordinate of the first monomer (>length_z/2)
    pol_init_bl - initial bond length of the polymer (0.5-1.0)

    n_steps - number of monte carlo moves to attempt (1000000)
    mc_ratio - percentage of mc moves that will choose a monomer over a solvent particle (0.5)
    temp - system temperature (1.0)
    freq_sample - sample the system every x steps (100)
    freq_monitor - monitor the system every x steps (100)

    iseed - random seed (change this before every run)

The run.sh and submit.sh scripts were written specifically for ACEnet clusters.

Changelog
---------

Version 3.0 (June 7, 2013)
* the wall has been added
    * new input parameter density_w: set density of wall particles
    * new input parameter n_layers: set number of wall layers
    * new parameter n_wall: number of wall particles
    * new parameter a_sw: interaction strength between solvent and wall particles
    * new parameter r_wall: distance between wall particles
* new function init_wall: calculates n_wall, allocates memory for both the wall and the solvent particles in the array part_dpd and initializes the wall particle positions
* renamed function setup_coords to init_part: initializes monomers and solvent particles
* new function check_wall: takes in a particle position vector r and checks if it is inside the wall
    * adjusts new parameter wall_overlap accordingly
    * the random move functions will now check for this overlap and immediately reject the move if necessary
* new function calc_cm: calculates center of mass of the polymer
* new function calc_rg: calculates radius of gyration of the polymer
* new input parameter mc_ratio: allows the user to adjust the frequency that monomers or solvent particles are chosen
    * values in the range 0 to 1
    * 0 = solvent particles moved only
    * 1 = monomers moved only
* new input parameter pol_init_z: sets the initial position of the first monomer of the chain
* new input parameter pol_init_bl: sets the initial bond length of the polymer
* split up sys.length into length.x, length.y and length.z to allow variable box length in each dimension
    * split up n_cell and r_cell as well to accomodate this change
* system volume calculation now subtracts the volume occupied by the wall
* modified the sample function to find absolute monomer positions (in the absence of periodic boundary conditions)
    * calc_re, calc_bond_length and calc_cm should now work properly if a periodic boundary is encountered
* renamed function periodic_bc to periodic_bc_dr: takes a displacement vector dr and adjusts according to periodic boundary conditions
* new function periodic_bc_r: takes a position vector r and adjusts according to periodic boundary conditions
* calc_pressure now accounts for wall-solvent and wall-monomer ineractions
* changed functions in init.c to store initial positions in part.r as well as part.ro
* removed a redundant call of new_list in main.c, as it is called in init.c already

Version 2.3 (May 31, 2013)
*   bug fixes
*   when old neighbour cells were considered, the head of chain was set to -1 to indicate that cell has already been counted. This caused calc_energy_dpd and calc_energy_mon to also not consider those cells, which was not intended
    * made a new 3D array to deal with this: sys.hoc_copy
    * can freely alter sys.hoc_copy under random_move functions without affecting calc_energy functions
    * in the future, may use memcpy function for efficiency, but getting segfaults with it now
*   if a move is accepted, particles will now overwrite the old energies with the new energies
    * that is, part.Eo = part.E
    * this ensures that if a move gets rejected and all the particle energies are reset, only those affected by the move will be reset to part.Eo
*   typo on lines mc.c:124-126
    * should be part_dpd, not part_mon
*   rewrote periodic_bc to use pointers
*   typo in calc_bond_length: i < sys.n_mon-2 should be i < sys.n_mon-1
*   new input parameters:
    * dr_max_dpd and dr_max_mon: set maximum displacement of a MC move
    * n_attempt_dpd and n_attempt_mon: counters for number of MC moves
    * n_accept_dpd and n_attempt_mon: counters for number of accepted MC moves
*   print_stats now prints selection and acceptance statistics of the monte carlo moves

Version 2.2 (May 27, 2013)
*   added monitor structure:
    * new variables: mon.energy, mon.re2, mon.rex, mon.rey, mon.rez, mon.bond_length, monitor_step
    * new function monitor_mem: to allocate memory for monitored quantities
    * new input parameter freq_monitor: note that this must be equal to freq_sample (may be changed in the future)
*   added some new sample quantities: energy, re2, re2x, re2y, re2z
    * new functions: calc_pressure, calc_re, calc_bond_length
    * calc_pressure was modified to account for polymer
*   new function periodic_bc: takes in dr vector and adjusts it based on periodic boundary conditions
*   new module energy.c: functions from calc.c and cell.c which related to energy were moved here
*   new function calc_energy_brute: brute force calculation of all monomer and dpd particle energies (deleted old calc_energy)
*   new function check_bond(i): checks for breakage in the bonds on either side of monomer i and adjust sys.bond_break accordingly
    * removed this same functionality from energy_fene, as it is now redundant
*   fixed a bug: conservative soft pair potentials were missing a factor of 1/2
*   renamed output to print_stats, a more descriptive name
*   adjusted write_mon to include new monitored quantities
*   fixed a bug: under monte_carlo, if a move was rejected, only one type of particle energy was reset (i.e. just dpd or just mon)
*   still to do: change calc_re to take into account periodic boundary conditions (might be tricky)

Version 2.1 (May 25, 2013)
*   some bug fixes and new functionality
*   rewrote the criteria for calculating new energies after a particle moves
    * check for new cell
    * find dix, diy, diz (difference between cell i and j)
    * using these, calculate energies for old neighboring cells
    * old method: 5*5*5 = 125 cells, every time
    * new method: 3*3*3 = 27 cells, if particle stays in the same cell
    *   = 27 + 9 cells, if one cell dimension changes (e.g. ix != ixo)
    *   = 27 + 9 + 6 cells, if two cell dimensions change (e.g. ix, iy != ixo, iyo)
    *   = 27 + 9 + 6 + 4 cells, if three cell dimensions change (e.g. ix, iy, iz != ixo, iyo, izo)
*   new function accept_move: returns 1 (true) if MC move is accepted, or 0 (false) if it is rejected
*   new function check_cell: returns 1 (true) if a particle has entered a new cell, 0 (false) otherwise
*   periodic boundary conditions for vector dr were incorrectly incrementing/decrementing by sys.length/2 instead of sys.length
*   monomers were not being initialized in the middle of the x-y plane

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

To-do list
----------

*   add the nanopore
*   add functionality to completely remove the wall and/or polymer using just input parameters
*   instead of using a global 3d array hoc_copy, use the memcpy function
