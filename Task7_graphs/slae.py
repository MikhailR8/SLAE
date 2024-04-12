import matplotlib.pyplot as pyplot
import numpy as np

ys1 = np.log(np.loadtxt('build/Conjugate_gradient.txt'))
ys5 = np.log(np.loadtxt('build/Simple_iter.txt'))
ys6 = np.log(np.loadtxt('build/Simple_iter_boost.txt'))
ys7 = np.log(np.loadtxt('build/Fastest_descent.txt'))
ys1 = ys1[0:ys1.size//2]
ys5 = ys5[0:ys5.size//2]
ys6 = ys6[0:ys6.size//2]
ys7 = ys7[0:ys7.size//2]
xs1 = np.arange(0, ys1.size, 1)
xs5 = np.arange(0, ys5.size, 1)
xs6 = np.arange(0, ys6.size, 1)
xs7 = np.arange(0, ys7.size, 1)

fig, ax = pyplot.subplots(figsize=(12, 8), dpi=400)

ax.minorticks_on()
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
ax.grid()
ax.grid(which='minor', linestyle=':')

ax.plot(xs1, ys1, zorder=10, linewidth=2, label='Метод Сопряжённых градиентов')
ax.plot(xs5, ys5, zorder=10, linewidth=2, label='Метод простой итерации')
ax.plot(xs6, ys6, zorder=10, linewidth=2, label='Метод простой итерации с ускорением')
ax.plot(xs7, ys7, zorder=10, linewidth=2, label='Наискорейший градиентный спуск')
ax.set_xlabel('Итерация', fontsize=14, fontname='Serif')
ax.set_ylabel('Логарифм невязки', fontsize=14, fontname='Serif')
ax.legend(loc='upper right')
pyplot.savefig("Graph 1.jpg", bbox_inches='tight')

