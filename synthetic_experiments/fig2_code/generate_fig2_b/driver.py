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
r_limit = 0.7
num_seg = 100
sigma = 0.05

if not os.path.isdir("./files"):
    os.system('mkdir files')
subprocess.call(['/bin/bash', '-i', '-c', 'g++ -larmadillo generate_fig2_b.cpp -O2 -o generate_fig2_b -DNDEBUG -framework Accelerate'])
os.system('./generate_fig2_b -D %d -d %d -N %d -r %f -num %d -sig %f' % (D, d, N, r_limit, num_seg, sigma))


filename = ['chatemax.ty', 'chatxmin.ty', 'etaO.ty', 'alpha.ty', 'beta.ty']

data = dict()
for file in filename:
    with open('./files/' + file, 'r') as f:
        vals = f.readlines()
    data[file] = []
    for each in vals:
        data[file].append(float(each.strip()))
    data[file] = np.array(data[file])


rr = np.linspace(0, r_limit, num_seg)

ff = plt.figure(figsize=(10, 8))
print(rr[0])

plt.plot(rr, (data['etaO.ty'] + D * 1) / data['chatxmin.ty'], 'g', linewidth=6)
plt.plot(rr, data['beta.ty'], 'r--', linewidth=6)
plt.legend([r'$R_{\mathcal{O}/\widehat{\mathcal{X}}}$',
            r'$R_{\widehat{\mathcal{E}}/\widehat{\mathcal{X}}}$'], fontsize=43, loc='best')
plt.xlabel(r'$M/(M+N)$', fontsize=40)

plt.xticks([0.1, 0.3, 0.5, 0.7], fontsize=40)
plt.yticks([0, 0.1, 0.2, 0.3, 0.4], fontsize=40)

plt.axis([0, 0.7, 0, 0.45])

plt.show()
ff.savefig("foo.pdf", bbox_inches='tight')
