1) Prior to running the "supp_simulation.py" script, ensure that all necessary dependencies are installed.

2) For optimal computational efficiency, it is recommended to utilize the Intel Fortran compiler. In cases where this is not feasible, the GNU Fortran compiler can be employed as a last resort.

3) Modify the 'n_proc' and 'compiler' variables located at the beginning of the specified script to specify the number of processors and the compiler to be used, respectively.

4) The script generates a set of directories, each identified by the PDE ID and resolution mentioned in their respective names. These directories house Fortran codes, data files, and visualizations associated with the simulation.