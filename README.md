# ClonOr_cpp

Clonal Origin program with multiple annealed importance sampling (AIS) steps for faster inference of homologous recombination in bacteria. This code is based on previous work: https://github.com/xavierdidelot/ClonalOrigin. Introducing AIS steps into a reversible jump MCMC (rjMCMC) leads to a reduction in asymptotic variance for the chain.

HOW TO RUN IN XCODE USING CPP FILES ONLY

	Download cpp files;

	File -> New project;

	Add cpp files;

	Download gsl;

	Configure GSL in building settings:
		OTHER_LDFLAGS = -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas
		HEADER_SEARCH_PATHS = /usr/local/include/
  
	Build project;

EXTRA SETTINGS FOR PARALLEL IMPLEMENTAION IN XCODE
	
	1. On terminal run "brew install llvm"
	2. Specify a user-defined setting CC with the following value 
		"/usr/local/opt/llvm/bin/clang".
	3. Add to the Library Search Paths the value 
		"/usr/local/opt/llvm/lib",
	to allow for the precompiled libraries (for example, libiomp5) to be dynamically linked to your executable.
	4. Add to the Header Search Paths the value 
		"/usr/local/opt/llvm/lib/clang/8.0.0/include",
	   this allows inclusion of headers provided by llvm installation (for example, omp.h). 
	   NB1 The installed version of Clang compiler might be different, so the value 8.0.0 of the Header Search Path should be changed accordingly. 
	   NB2 We already added the GSL configuration to Header Search Paths, separate values should be separated by a semicolo, i.e. the field should look like: 
		"/usr/local/include/; /usr/local/opt/llvm/lib/clang/8.0.0/include"
	5. Add the "-fopenmp" compilation flag to Other C Flags
	6. Link llvm OpenMP library dynamically to your executable.
	
	The code should compile now. 
	The error
		"clang: error: cannot specify -o when generating multiple output files" 
	might pop up. To fix this go to build settings > build options > Enable Index-While-Building Functionality and set the value to "No".
	Sources and more details:
		http://antonmenshov.com/2017/09/09/clang-openmp-setup-in-xcode/
		https://stackoverflow.com/questions/46527662/xcode-9-clang-error-cannot-specify-o-when-generating-multiple-output-files
	


Running the program with argument -h shows help file with arguments needed. E.g. the program will run introducing the following arguments:
	
	-a 1,1,1,2,2,1,1,1,0,0,0 -x 0 -y 500000 -z 10 -D50 -T 10 -R 5 -m 0.25,100,10 true_tree.nwk simulatedData.xmfa test.xml
	where
	-m requires 3 values (gamma, T, N) of types (double, int, int). 
	   The value:
		gamma is the power used in the annealing computation, gamma>1 means the steps adding a recombination concentrate near the start, gamma=1 means that the steps are equidistant;
        	T corresponds to the number of steps used in the AIS procedure for adding or deleting a recombination.
		N denotes the number of replications used in the MAIS procedure for adding or deleting a recombination. 
	   NB The parameter N is the one introducing the parallelisation, if you want run the code with N=1 (which corresponds to a pure AIS-RJMCMC) you should use N=0, which will run the code without a parallel implementation.
			
	


IN THE CLUSTER: 
	
	1. Compile cpp files using
		g++ *.cpp -o <<NAME>> -lgsl -lgslcblas -lm -fopenmp
	2. Then run
		./<<NAME>> -a 1,1,1,2,2,1,1,1,0,0,0 -x 0 -y 500000 -z 10 -D50 -T 10 -R 5 -m 0.25,100,1 true_tree.nwk simulatedData.xmfa test.xml
