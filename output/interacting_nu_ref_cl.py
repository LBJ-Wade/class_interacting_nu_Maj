import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Dropbox/class_interacting_nu_Maj/output/interacting_nu_ref_cl.dat', '/Users/poulin/Dropbox/class_interacting_nu_Maj/output/interacting_nu_test_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['interacting_nu_ref_cl', 'interacting_nu_test_cl']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 1, data[1]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()