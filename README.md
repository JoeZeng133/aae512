The code demonstrates mesh generation and transient heat conduction problem in a single connected domain.

1. COMPILING
	A cmake of version higher than 3.9 is required. Make sure the compiler supports c++11 standard.
	
2. DEMO
	In order to succesuffly run the two demo programs (solver_demo.m, mesh_demo.m) in the matlab folder, make sure the executable is placed in the matlab folder as well.

3. STANDALONE RUNNING
	command 1:main.exe boundaryfilelocation meshfilelocation
		reads boundary information from boundaryfilelocation and generates mesh to be stored in meshfilelocation
	command 2:main.exe meshfilelocation configurationlocation outputlocation
		reads mesh from meshfilelocation, configuration from configurationlocation and calculates temperature distribution to be stored in outputlocation