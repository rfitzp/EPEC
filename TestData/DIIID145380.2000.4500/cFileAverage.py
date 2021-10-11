'''
Python script to average all cFiles
'''

import glob
import matplotlib.pyplot as plt

cFiles = glob.glob("cFiles/c*")
nfiles = len(cFiles)

nrows = 100

column1 = [0.]*nrows
column2 = [0.]*nrows
column3 = [0.]*nrows
column4 = [0.]*nrows
column5 = [0.]*nrows

for f in cFiles:
    
    f = open('cFiles/c145380.4000','r')

    lines=f.readlines()[1:]

    data1 = []
    data2 = []
    data3 = []
    data4 = []
    data5 = []
    for x in lines:
        data1.append(float(x.split()[0]))
        data2.append(float(x.split()[1]))
        data3.append(float(x.split()[2]))
        data4.append(float(x.split()[3]))
        data5.append(float(x.split()[4]))

    f.close()

    for i in range(0,nrows):
        column1[i] = column1[i] + data1[i]
        column2[i] = column2[i] + data2[i]
        column3[i] = column3[i] + data3[i]
        column4[i] = column4[i] + data4[i]
        column5[i] = column5[i] + data5[i]

for i in range(0,nrows):
    column1[i] = column1[i] /float(nfiles)
    column2[i] = column2[i] /float(nfiles)
    column3[i] = column3[i] /float(nfiles)
    column4[i] = column4[i] /float(nfiles)
    column5[i] = column5[i] /float(nfiles)

lim = 10.;

for i in range(0,nrows):
    if column1[i] > lim:
        column1[i] = lim
    if column2[i] > lim:
        column2[i] = lim
    if column3[i] > lim:
        column3[i] = lim
    if column4[i] > lim:
        column4[i] = lim
    if column5[i] > lim:
        column5[i] = lim

f = open('cFile','w')
f.write('%4d\n' % nrows)
for i in range(0, nrows):
    f.write('%11.4e %11.4e %11.4e %11.4e %11.4e\n' % (column1[i], column2[i], column3[i], column4[i], column5[i]))
f.close()    

fig = plt.figure(figsize=(16.0, 8.0))
fig.canvas.manager.set_window_title("Average cFile for DIII-D Shot 145380")

plt.subplot(2, 2, 1)
plt.xlim(0.0, 1.0)
plt.plot(column1, column2, 'r--', linewidth=0.4)
plt.plot(column1, column2, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel('$\\chi_\phi(m^2/s)$')

plt.subplot(2, 2, 2)
plt.xlim(0.0, 1.0)
plt.plot(column1, column3, 'r--', linewidth=0.4)
plt.plot(column1, column3, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\chi_e(m^2/s)$")

plt.subplot(2, 2, 3)
plt.xlim(0.0, 1.0)
plt.plot(column1, column4, 'r--', linewidth=0.4)
plt.plot(column1, column4, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$D_\\perp(m^2/s)$")

plt.subplot(2, 2, 4)
plt.xlim(0.0, 1.0)
plt.plot(column1, column5, 'r--', linewidth=0.4)
plt.plot(column1, column5, 'bo', markersize=5)
plt.xlabel('$\\Psi_N$')
plt.ylabel("$\\chi_i(m^2/s)$")

plt.tight_layout(pad=0.5)

plt.show()
