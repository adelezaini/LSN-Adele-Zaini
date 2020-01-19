LSN_exercises_delivery

––––––––––––––––––– BASIC FILE ORGANISATION OF EACH FOLDER –––––––––––––––––––––––

- Jupyter Notebook containing an explanation of the exercises and representations of the results (eventually some –commented–Python scripts to manage the input/output files and run the executable quickier): Exercises*.ipynb
- source code: main.cpp (or analogous names, such as MolDyn_NVE.cpp, Monte_Carlo_ISING_1D.cpp...)
- files for the Random Number Generator: random.h, random.cpp, Primes, seed.in
- output files with simulation results: \*.out or output.\* (eventually contained in proper folders)
- Makefile
- executable: main.exe (or MolDyn_NVE.exe, Monte_Carlo_ISING_1D.exe...)

In some exercises you will also find:
- README when the file organisation needs a further explanation (or how to equilibrate the system, etc...)
- header files: main.h, MolDyn_NVE.h, Monte_Carlo_ISING_1D.h, TSP.h...
- input files with simulation parameters: input.dat (or input.*)
- output files with configurations: config.*
