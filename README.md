# Code of paper "Noisy Dual Principal Component Pursuit", ICML 2019

## Synthetic Experiments

- Requirements

	- Matlab
	- C++
	- [Armadillo](http://arma.sourceforge.net/download.html) (C++ scientific computing library)
	- Python3
	- numpy, matplotlib
	- [optional] [OpenBLAS](http://www.openblas.net)

	Basically, we use MATLAB for simple tasks, C++ for heavy jobs and Python to do some plotting after obtaining data through MATLAB/C++ programs.
	
- Usage (tested under Mac OS)

	- Matlab program can be executed directly
	- Generally, each `.cpp` file is paired with a `driver.py` file (just run the driver file is enough)
	- An installation of OpenBLAS will further accelerate the C++ programs but the compilation in `driver.py` needs to be changed accordingly

	

	