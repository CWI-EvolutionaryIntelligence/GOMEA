# CONTRIBUTE 

This page is now incomplete, but will include information on how to contribute to the GOMEA library.

## Notes on including sub-module
- Default constructor should be available.
- Destructor should not free objects that are not initialized in the default constructor.
- No exit statements should be used throughout the code.
- No output to stdout and no files.
- Options should be passed to the C++ class through a C++ config class.
- Default parameters should be uniquely defined in the Cython .pyx file wrapping the sub-module. These default parameters are then used to initialize the C++ config class.
