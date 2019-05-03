# ClonOr_cpp

Clonal Origin program with multiple annealed importance sampling (AIS) steps for faster inference of homologous recombination in bacteria. This code is based on previous work: https://github.com/xavierdidelot/ClonalOrigin. Introducing AIS steps into a reversible jump MCMC (rjMCMC) leads to a reduction in asymptotic variance for the chain.

HOW TO RUN IN XCODE USING CPP FILES ONLY

	Download cpp files;

	File -> New project;

	Add cpp files;

	Download gsl;

	Configure in building settings:
		OTHER_LDFLAGS = -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas
		HEADER_SEARCH_PATHS = /usr/local/include/
  
	Build project;
	
Running the program with argument -h shows help file with arguments needed. E.g. the program will run introducing the following arguments:
	
	-a 1,1,1,2,2,1,1,1,0,0,0 -x 0 -y 500000 -z 10 -D50 -T 10 -R 5 -A 300 true_tree.nwk simulatedData.xmfa test.xml
	where
	-A denotes the desired number of annealing steps.


IN THE CLUSTER: compile cpp files using
	g++ *.cpp -o "NAME" -lgsl -lgslcblas -lm
then run
	./"NAME" -a 1,1,1,2,2,1,1,1,0,0,0 -x 0 -y 500000 -z 10 -D50 -T 10 -R 5 -A 300 true_tree.nwk simulatedData.xmfa test.xml
