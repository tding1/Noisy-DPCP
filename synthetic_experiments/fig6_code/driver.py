import subprocess
import os

subprocess.call(['/bin/bash', '-i', '-c', 'g++ -std=c++11 -larmadillo generate_data.cpp -O2 -o generate_data -DNDEBUG -framework Accelerate'])

D = 30
d = 29
N = 1000

os.system('./generate_data -D %d -d %d -N %d' % (D, d, N))
