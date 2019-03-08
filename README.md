# ClonOr_cpp

Clonal Origin program with multiple AIS for faster inference of homologous recombination in bacteria.

HOW TO COMPILE USING CPP FILES ONLY

Download cpp files;
File -> New project;
Add cpp files;
Download gsl;
Configure in building settings:
	OTHER_LDFLAGS = -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas
	HEADER_SEARCH_PATHS = /usr/local/include/
  
