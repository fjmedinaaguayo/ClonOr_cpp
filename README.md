# ClonOr_cpp

Clonal Origin program with multiple AIS for faster inference of homologous recombination in bacteria.

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
