import matplotlib.pyplot as pyplot
import numpy as np

n = 25
input1 = np.loadtxt('build/Conjugate_gradient.txt')
input2 = np.loadtxt('build/GMRES.txt')
ys1 = np.log(input1[0:input1.size//2])
ys2 = np.log(np.fabs(input2[0:n]))
xs1 = input1[input1.size//2:]
xs2 = input2[n:-2]
xs2 = np.append(xs2, (input2[-1]))
ys2 = np.append(ys2, (np.log(input2[-2])))
xs1 = np.cumsum(xs1)
xs2 = np.cumsum(xs2)
print(ys2)

fig, ax = pyplot.subplots(figsize=(12, 8), dpi=400)

ax.minorticks_on()
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
ax.grid()
ax.grid(which='minor', linestyle=':')

ax.plot(xs1, ys1, zorder=10, linewidth=2, label='Метод сопряжённых градиентов (CG)')
ax.plot(xs2, ys2, zorder=10, linewidth=2, label='GMRES')
ax.set_xlabel('Время, мкс', fontsize=14, fontname='Serif')
ax.set_ylabel('Логарифм невязки', fontsize=14, fontname='Serif')
ax.legend(loc='upper right')
ax.text(500, 0, "Последний прямой участок GMRES - время решения СЛАУ гауссом")
pyplot.savefig("Graph 2.jpg", bbox_inches='tight')

