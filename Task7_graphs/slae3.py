import matplotlib.pyplot as pyplot
import numpy as np

ys1 = np.loadtxt('build/Conjugate_gradient.txt')
ys5 = np.loadtxt('build/Simple_iter.txt')
ys6 = np.loadtxt('build/Simple_iter_boost.txt')
ys7 = np.loadtxt('build/Fastest_descent.txt')
xs1 = np.cumsum(ys1[ys1.size//2+1:])
xs5 = np.cumsum(ys5[ys5.size//2:])
xs6 = np.cumsum(ys6[ys6.size//2+1:])
xs7 = np.cumsum(ys7[ys7.size//2:])

ys1 = np.log(np.loadtxt('build/Conjugate_gradient.txt'))
ys5 = np.log(np.loadtxt('build/Simple_iter.txt'))
ys6 = np.log(np.loadtxt('build/Simple_iter_boost.txt'))
ys7 = np.log(np.loadtxt('build/Fastest_descent.txt'))
ys1 = ys1[0:ys1.size//2]
ys5 = ys5[0:ys5.size//2]
ys6 = ys6[0:ys6.size//2]
ys7 = ys7[0:ys7.size//2]

fig, ax = pyplot.subplots(figsize=(12, 8), dpi=400)

ax.minorticks_on()
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
ax.grid()
ax.grid(which='minor', linestyle=':')

ax.plot(xs1, ys1, zorder=10, linewidth=2, label='Метод Сопряжённых градиентов')
ax.plot(xs5, ys5, zorder=10, linewidth=2, label='Метод простой итерации')
ax.plot(xs6, ys6, zorder=10, linewidth=2, label='Метод простой итерации с ускорением')
ax.plot(xs7, ys7, zorder=10, linewidth=2, label='Наискорейший градиентный спуск')
ax.set_xlabel('Время, мкс', fontsize=14, fontname='Serif')
ax.set_ylabel('Логарифм невязки', fontsize=14, fontname='Serif')
ax.legend(loc='upper left')
pyplot.savefig("Graph 3.jpg", bbox_inches='tight')