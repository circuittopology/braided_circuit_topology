# braided_circuit_topology
Braided Circuit Topology framework for multiple polymer molecules analysis 

Algorithmic flow of the different parts of the simulation:
First compile the C programs by running the Makefile

<ol>
  <li> Initialise polymers by running 'initial_condition.exe #polymers #monomers per polymer box_y_size box_z_size' </li>
  This produces a file called 'chains_ic.data'
  <li> LAMMPS simulation by running 'main_batch.topo' in LAMMPS with the bending stiffness as argument. </li>
  This produces trajectory information, as well as radius of gyration, etc.
  <li> Perform Circuit Topology analysis by running 'CT_analysis_batch.exe cutoff_start cutoff_end nr_of_runs' </li>
  This produces a file 'CT_content_K%0.1f.txt' containing the CT motif fractions for the chosen stiffness and number of runs
  <li> Reshape files for PPA analysis by running 'reshape.exe nr_of_runs' </li>
  Reshapes the trajectories for using them in the PPA analysis.
  <li> LAMPPS PPA Analysis by running 'main_batch.ppa' in LAMMPS without arguments. </li>
  Produces new trajectory data that is subsequently used for the braid analysis.
  <li> Do a Braid Analysis by running the Matlab script 'entanglement.m'. WARNING: this requires that the BraidLab package is installed and that the function 'data_to_braid.m' is in the Matlab PATH. </li>
</ol>

To facilitate use and to allow for the systematic scanning of a stiffness range $[\kappa_-, \kappa_+]$, the BASH script can be run. Make sure to change the number of runs/iterations per stiffness value in the script, as well as in main.topo, main.ppa, reshape.c and CT_analysis_batch.c, and then recompile using the Makefile.
