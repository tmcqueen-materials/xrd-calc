# xrd-calc
Anomalous X-ray Scattering Quantitative Analysis Toolkit
Journal Reference: J. Am. Chem. Soc. 132, 16185-16190 (2010)

This page serves as a repository for all of the necessary computer code (written in C) to extract quantitative information regarding elemental occupancy of different crystallographic sites from a series of single crystal x-ray datasets collected near the absorption edges of the elements.

To use the code, download the source. Modify Makefile to match the local compiler settings (RHEL 5 / CentOS 5 should require no changes, and see below for OS X 10.6; other operating systems may require changes).

Build the executables using "make all", and then run search_f on each of the x-ray datasets (see paper for analysis details). More complex problems may require you to modify the routine in search_f, a la the example code in search_hbs_f. At some point in the future, if there is interest, a general-purpose program and GUI will be created based on this code.

Important Notes

Beyond atomic positions, the most important ins parameters are the x-ray wavelength on the CELL line and DISP commands for every element (see next item).

This program does not currently include functionality to automatically determine appropriate initial anomalous scattering factors. Thus you must provide the appropriate DISP lines giving (at least) f' and f'' for all elements present, at the measured wavelength. If you have portable code that would let me include this feature, let me know!

In general, the best quantitation is achieved using an interative process: assume no site mixing and run search_f at all wavelengths. Following the paper method, compute the level of mixing on each site. Then construct new ins files that include this level of site mixing, and use as new inputs for search_f. Repeat until the values are converged. For elements that have very similar electron counts (e.g. Zn and Cu), no iteration is, in practice, needed.

Mac OS X 10.6 Installation Instructions

1. Install XCode (use the OS X installer disc, or register and download from http://developer.apple.com/technologies/xcode.html).
2. Download the source files into a folder called xrd-calc.
3. Open a new terminal window. This is best done by opening "X11" (Applications->Utilities->X11), which opens an xterm.
4. Goto the xrd-calc directory: cd ~/Desktop/xrd-calc
5. Edit the file Makefile:
   emacs Makefile
   (you can use another editor; a brief emacs tutorial is here: http://www2.lib.uchicago.edu/keith/tcl-course/emacs/tutorial.html): Delete "/usr/lib64/libm.so" from the COMMON line. Save and exit.
6. Build the programs. Type: make
   You should see a bunch of lines like:
   gcc -O3 -c -o xrd_calc.o xrd_calc.c
   gcc -O3 -c -o compute.o compute.c
   gcc -O3 -c -o read.o read.c
   ...
   And no error messages (example error message: make: [xrd_calc] Error 1)
7. Assuming no error messages, you will be back at the command prompt. Run the tests to make sure everything is working:
   ./test_main | more
   All the tests should pass (the vertical line in a pipe symbol, above the \ key)
8. At this point you should be able to run search_f. A typical run would be:
   ./search_hbs_f hbs.ins hbs.hkl
