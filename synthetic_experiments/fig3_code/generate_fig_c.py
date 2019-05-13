import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

x = np.loadtxt('x_val.txt')
y = np.loadtxt('y_val.txt')
t1 = np.loadtxt('t1_val.txt')
t2 = np.loadtxt('t2_val.txt')


fig = plt.figure(figsize=(15, 10))
ax = fig.add_subplot(111)

# choose c = t1 or t2 to plot fig a or b
cs = plt.scatter(x, y, s=5, marker='o', c=t1, cmap=cm.Greys)
plt.axis([0, 1, 0, 0.25])
plt.xlabel(r'$R_{\mathcal{O}/\widehat{\mathcal{X}}}$', fontsize=60)
plt.ylabel(r'$R_{\widehat{\mathcal{E}}/\widehat{\mathcal{X}}}$', rotation=0, fontsize=60)
plt.xticks(fontsize=60)
plt.yticks(fontsize=60)
t = plt.colorbar(cs)
t.ax.tick_params(labelsize=40)
ax.yaxis.set_label_coords(-0.15, 0.45)
plt.xticks([0, 0.5, 1])
plt.yticks([0.25])

plt.show()

fig.savefig("foo.pdf", bbox_inches='tight')
