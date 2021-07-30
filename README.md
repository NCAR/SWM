SWM is a simplified kernel representing the nonlinear PDE's governing geophysical fluid flow.
Please contact the Point of Contact below before distributing it further, or if any questions or issues arise.

COPY vs SWAP:

There is an optimization that avoids copies in loop 300 by swapping pointers. COPY can be turned on and off by defining _COPY_ or not, as desired.

Building:

Use __./buildRun*.sh__ to build and run your desired version of SWM. Use __./buildRun*.sh__ __-h__ to see additional optinos for building and runs.

Validation: 

Compare results.txt with the provided file results_ref.txt, particularly the last values printed, namely those below-
 diagonal elements of p
 diagonal elements of u
 diagonal elements of v

Performance: 

Performance reported in results_ref.txt is that obtained on a single core of a 2.6 GHz Intel Core i7 Macbook Pro with 2400 MHz DDR4 memory.

Point of Contact:

Dr. Richard Loft
National Center for Atmospheric Research
loft@ucar.edu
