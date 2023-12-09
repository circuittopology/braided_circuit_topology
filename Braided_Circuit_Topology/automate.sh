#!/bin/bash

start=`date +%s`
# WATCH OUT!!! Change nr. of iterations in main.topo, main.ppa, here in reshape.exe and CT_analysis_batch.exe


# Uncomment to remove all trailing initial data files (clean up)
#rm braid_measures_*.txt
#rm RG_*.txt
#rm CT_content_*.txt

# Change start and stop stiffness values, as well as increment in seq loop

for k in $(seq 11.2 0.2 12.0)
do
	echo $k > ktmp.txt
	echo "Running LAMMPS simulation for k = $k"
	mpirun ./lmp_mpi -in main_batch.topo -var bending_stiff $k
	echo "Performing Circuit Topology analysis..."
	./CT_analysis_batch.exe 1.4 1.4 1000 $k
	echo "Reshaping files..."
	./reshape.exe 1000
	echo "Running PPA algorithm for k = $k"
	mpirun ./lmp_mpi -in main_batch.ppa
	echo "Starting Matlab entanglement..."
	matlab -batch "entanglement"
	echo "All done, enjoy!"
done

rm RG_*.txt

echo "Full loop completed."

end=`date +%s`

runtime=$((end-start))

echo $runtime
