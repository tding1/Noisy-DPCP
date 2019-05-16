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

	
## 3D Roadplane Estimation

- `demo.m` is a toy example that runs single subspace learning algorithms on real 3D road plane detection data. Once we have selected the frame and click run button, it instantly runs the algorithms and returns the clustering metrics, geometric metrics and algorithmic metrics as mentioned in the paper. Also, a poster showing the projections of the separated point clouds onto the image is generated after the execution of the program.

- `/data` is a folder containing annotations for point clouds and corresponding images.

- `/algorithms` is a folder containing various single subspace learning algorithms.