# Script to automate gmx

for dir in */; do
    cd $dir
    echo "$dir"
    mkdir MMPBSA_results
    for run_dir in 2_run 2_run_second 2_run_third; do
	cd "$run_dir"
	echo "$run_dir"
    	cp ../../mmpbsa.in .
    	mpirun -np 120 gmx_MMPBSA MPI -i mmpbsa.in -cs md.tpr -ci index_mmpbsa_adjusted.ndx -cg 28 29 -ct md_mmpbsa.xtc -cp topol.top -nogui &
    	wait $!
    	dirname=$(basename `pwd`)
    	cp FINAL_RESULTS_MMPBSA.dat MMPBSA_$dirname.dat
    	cp MMPBSA_*.dat ../MMPBSA_results
    	cd ..
    done
    cd ..
done

echo "The automated calculation has finished"
