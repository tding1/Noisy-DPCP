import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

M = np.loadtxt('M.txt')
N = np.loadtxt('N.txt')
cosphi = np.loadtxt('cos_phi.txt')

fig = plt.figure(figsize=(10, 8))

cs = plt.scatter(N, M, s=100, marker='s', c=cosphi, cmap=cm.Greys)
x = np.linspace(0, 500, 250)
y = 0.0083 * x**2 + 0.8333 * x
plt.plot(x, y, 'r', linewidth=6)
plt.text(330, 750, r'$M=O(N^2)$', fontsize=30)
plt.arrow(330, 820, -25, 80, head_width=15, head_length=50, fc='k', ec='k')
plt.axis([0, 500, 0, 3000])
plt.xlabel(r'$N$', fontsize=40)
plt.ylabel(r'$M$', rotation=0, fontsize=40)
plt.xticks(fontsize=40)
plt.yticks([0, 1000, 2000, 3000], fontsize=40)
t = plt.colorbar(cs)
t.ax.tick_params(labelsize=40)

plt.show()

fig.savefig("foo.pdf", bbox_inches='tight')
