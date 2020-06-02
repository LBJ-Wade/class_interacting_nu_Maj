import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Dropbox/class_interacting_nu_Maj/output/interacting_nu_test_perturbations_k0_s.dat', '/Users/poulin/Dropbox/class_interacting_nu_Maj/output/interacting_nu_ref_perturbations_k0_s.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['interacting_nu_test_perturbations_k0_s', 'interacting_nu_ref_perturbations_k0_s']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'delta_ncdm[0]', u'delta_ncdm[1]']
tex_names = ['delta_ncdm[0]', 'delta_ncdm[1]']
x_axis = 'tau [Mpc]'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 17]))
ax.loglog(curve[:, 0], abs(curve[:, 21]))

index, curve = 1, data[1]
y_axis = [u'delta_ncdm[0]', u'delta_ncdm[1]']
tex_names = ['delta_ncdm[0]', 'delta_ncdm[1]']
x_axis = 'tau [Mpc]'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 17]))
ax.loglog(curve[:, 0], abs(curve[:, 21]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('tau [Mpc]', fontsize=16)
plt.show()