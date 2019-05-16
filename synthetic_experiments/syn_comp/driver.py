import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt

D = 30
d = 29
N = 1000
r = 0.1
M = int(r * N / (1 - r))
num_seg = 100
sigma_limit = 0.2


if not os.path.isdir("./files"):
    os.system('mkdir files')
subprocess.call(['/bin/bash', '-i', '-c', 'g++ -std=c++11 -larmadillo compare_PSGD_REAPER_GGD.cpp -O2 -o comp -DNDEBUG -framework Accelerate'])
os.system('./comp -D %d -d %d -N %d -r %f -num %d -sig %f' % (D, d, N, r, num_seg, sigma_limit))

fileOpt1 = ['cos_phi_PSGD.ty', 'cos_phi_REAPER.ty', 'cos_phi_GGD.ty']

filename = fileOpt1

data = dict()
for file in filename:
    with open('./files/' + file, 'r') as f:
        vals = f.readlines()
    data[file] = []
    for each in vals:
        data[file].append(float(each.strip()))
    data[file] = np.array(data[file])


sigma = np.linspace(0, sigma_limit, num_seg)

if filename == fileOpt1:
    fig = plt.figure(figsize=(10, 8))
    plt.plot(sigma, data['cos_phi_PSGD.ty'], 'r', linewidth=3)
    plt.plot(sigma, data['cos_phi_REAPER.ty'], 'g', linewidth=3)
    plt.plot(sigma, data['cos_phi_GGD.ty'], 'k', linewidth=3)
    plt.xlabel(r'$\sigma$', fontsize=35)
    plt.ylabel(r'$\cos(\phi_*)$', fontsize=35)
    plt.legend(['PSGD', 'REAPER', 'GGD'], fontsize=35)

    plt.xticks([0, 0.05, 0.1, 0.15, 0.2], fontsize=35)
    plt.yticks([0, 0.02, 0.04, 0.06], fontsize=35)

    fig.savefig("foo.pdf", bbox_inches='tight')

plt.show()
