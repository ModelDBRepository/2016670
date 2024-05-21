all: Knclmp

Knclmp: Kncl.cpp current_classes.h 
	g++ Kncl.cpp -O3 -lm -fopenmp -o Knclmp

run_baseline: Knclmp
	./Knclmp inpNet -clusters sc_dmn10_rat -cluster_scale 0.35 -randseed 10 

run_ach: Knclmp
	./Knclmp inpNet -clusters sc_dmn10_rat -cluster_scale 0.35 -p_gKld 0.04324 -network_ampa_file dmn_ampa -randseed 10 
	 
clean:
	-rm Knclmp micK avg_V connectivity conn rconn time_* spikes py_rates py_in logfile in_py in_rates E_py_py *.eps *.fig *.mat STDIN.* 2>/dev/null


