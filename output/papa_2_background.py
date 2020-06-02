import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/class_public/output/papa_2_background.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['papa_2_background']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = ['conf.time[Mpc]', 'propertime[Gyr]']
tex_names = ['proper time [Gyr]', 'conf. time [Mpc]']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 2]))
ax.loglog(curve[:, 0], abs(curve[:, 1]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()