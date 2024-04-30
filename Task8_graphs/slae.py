import matplotlib.pyplot as pyplot
import numpy as np

n = 25
input1 = np.loadtxt('build/Conjugate_gradient.txt')
input2 = np.loadtxt('build/GMRES.txt')
ys1 = np.log(input1[0:input1.size//2])
ys2 = np.log(np.fabs(input2[0:n]))
xs1 = np.arange(0, ys1.size, 1)
xs2 = np.arange(0, ys2.size, 1)
print(ys2)

fig, ax = pyplot.subplots(figsize=(12, 8), dpi=400)

ax.minorticks_on()
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
ax.grid()
ax.grid(which='minor', linestyle=':')

ax.plot(xs1, ys1, zorder=10, linewidth=2, label='Метод сопряжённых градиентов (CG)')
ax.plot(xs2, ys2, zorder=10, linewidth=2, label='GMRES')
ax.set_xlabel('Итерация', fontsize=14, fontname='Serif')
ax.set_ylabel('Логарифм невязки', fontsize=14, fontname='Serif')
ax.legend(loc='upper right')
pyplot.savefig("Graph 1.jpg", bbox_inches='tight')

