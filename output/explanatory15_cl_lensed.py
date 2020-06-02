import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/class_public/output/explanatory15_cl_lensed.dat', '/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/AxiCLASS_merging/output/explanatory11_cl_lensed.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['explanatory15_cl_lensed', 'explanatory11_cl_lensed']

fig, ax = plt.subplots()
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()