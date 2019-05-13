import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

D = 30
N = 1500
num_seg = 200
sigma_limit = 0.1

if not os.path.isdir("./files"):
    os.system('mkdir files')

subprocess.call(['/bin/bash', '-i', '-c', 'g++ -std=c++11 -larmadillo generate_fig7_a.cpp -O2 -o generate_fig7_a -DNDEBUG -framework Accelerate'])

d = 29
r = 0.1
os.system('./generate_fig7_a -D %d -d %d -N %d -r %f -num %d -sig %f' % (D, d, N, r, num_seg, sigma_limit))
os.system('mv ./files/t1.ty ./files/t1_d29_o0.1.ty')

d = 29
r = 0.7
os.system('./generate_fig7_a -D %d -d %d -N %d -r %f -num %d -sig %f' % (D, d, N, r, num_seg, sigma_limit))
os.system('mv ./files/t1.ty ./files/t1_d29_o0.7.ty')

d = 1
r = 0.1
os.system('./generate_fig7_a -D %d -d %d -N %d -r %f -num %d -sig %f' % (D, d, N, r, num_seg, sigma_limit))
os.system('mv ./files/t1.ty ./files/t1_d1_o0.1.ty')

d = 1
r = 0.7
os.system('./generate_fig7_a -D %d -d %d -N %d -r %f -num %d -sig %f' % (D, d, N, r, num_seg, sigma_limit))
os.system('mv ./files/t1.ty ./files/t1_d1_o0.7.ty')


filename = ['t1_d1_o0.1.ty', 't1_d1_o0.7.ty', 't1_d29_o0.1.ty', 't1_d29_o0.7.ty']

data = dict()
for file in filename:
    with open("./files/" + file, 'r') as f:
        vals = f.readlines()
    data[file] = []
    for each in vals:
        data[file].append(float(each.strip()))
    data[file] = np.array(data[file])

num_seg = 200
sigma_limit = 0.1
sigma = np.linspace(0, sigma_limit, num_seg)
ind = [2 * i for i in range(100)]

LINEWIDTH = 6
FIG_SIZE = (10, 8)
FONT_SIZE = 40
LEGEND_SIZE = 28
ff = plt.figure(figsize=FIG_SIZE)
plt.plot(sigma[ind], data['t1_d1_o0.1.ty'][ind], '-', linewidth=LINEWIDTH)
plt.plot(sigma[ind], data['t1_d1_o0.7.ty'][ind], ':', linewidth=LINEWIDTH)
plt.plot(sigma[ind], data['t1_d29_o0.1.ty'][ind], '--', linewidth=LINEWIDTH)
plt.plot(sigma[ind], data['t1_d29_o0.7.ty'][ind], '-.', linewidth=LINEWIDTH)

plt.xlabel(r'$\sigma$', fontsize=FONT_SIZE)
plt.ylabel(r'$t_1$', fontsize=FONT_SIZE, rotation=0)
plt.legend([r'$d=1, \frac{M}{N+M}=0.1$', r'$d=1, \frac{M}{N+M}=0.7$',
            r'$d=29, \frac{M}{N+M}=0.1$', r'$d=29, \frac{M}{N+M}=0.7$'], fontsize=LEGEND_SIZE, loc='best')

plt.xticks([0, 0.05, 0.1], fontsize=FONT_SIZE)
plt.yticks([0.1, 0.2, 0.3, 0.4], fontsize=FONT_SIZE)
plt.axis([0, 0.1, 0, 0.45])

plt.show()

ff.savefig("foo.pdf", bbox_inches='tight')
