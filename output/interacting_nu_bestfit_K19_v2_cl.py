import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/class_public/output/interacting_nu_bestfit_K19_v2_cl.dat', '/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/class_interacting_nu/output/interacting_nu_bestfit_K19_v2_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['interacting_nu_bestfit_K19_v2_cl', 'interacting_nu_bestfit_K19_v2_cl']

fig, ax = plt.subplots()
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()