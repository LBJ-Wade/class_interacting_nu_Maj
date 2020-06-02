import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/class_public/output/explanatory19_pk_nl.dat', '/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/class_public/output/explanatory20_pk_nl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['explanatory19_pk_nl', 'explanatory20_pk_nl']

fig, ax = plt.subplots()
y_axis = [u'P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ax.set_xlabel('k (h/Mpc)', fontsize=16)
plt.show()