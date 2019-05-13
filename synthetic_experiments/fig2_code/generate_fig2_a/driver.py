import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

D = 30
d = 29
N = 1500
r = 0.7
M = int(r * N / (1 - r))
num_seg = 100
sigma_limit = 0.1

if not os.path.isdir("./files"):
    os.system('mkdir files')
subprocess.call(['/bin/bash', '-i', '-c', 'g++ -larmadillo generate_fig2_a.cpp -O2 -o generate_fig2_a -DNDEBUG -framework Accelerate'])
os.system('./generate_fig2_a -D %d -d %d -N %d -r %f -num %d -sig %f' % (D, d, N, r, num_seg, sigma_limit))

filename = ['chatemax.ty', 'chatxmin.ty', 'etaO.ty', 'alpha.ty', 'beta.ty']

data = dict()
for file in filename:
    with open('./files/' + file, 'r') as f:
        vals = f.readlines()
    data[file] = []
    for each in vals:
        data[file].append(float(each.strip()))
    data[file] = np.array(data[file])


sigma = np.linspace(0, sigma_limit, num_seg)

b = N * sigma
ff = plt.figure(figsize=(10, 8))
plt.plot(sigma, (data['etaO.ty'] + D * 1) / data['chatxmin.ty'], 'g', linewidth=6)
plt.plot(sigma, data['beta.ty'], 'r--', linewidth=6)
plt.legend([r'$R_{\mathcal{O}/\widehat{\mathcal{X}}}$',
            r'$R_{\widehat{\mathcal{E}}/\widehat{\mathcal{X}}}$'], fontsize=45, loc='best')
plt.xlabel(r'$\sigma$', fontsize=40)

plt.xticks([0, 0.05, 0.1], fontsize=40)
plt.yticks([0.1, 0.2, 0.3, 0.4], fontsize=40)
plt.axis([0, 0.1, 0, 0.41])

plt.show()
ff.savefig("foo.pdf", bbox_inches='tight')
