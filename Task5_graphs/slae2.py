import matplotlib.pyplot as pyplot
import numpy as np

ys1 = np.cumsum(np.loadtxt('result5.txt'))
ys2 = np.cumsum(np.loadtxt('result6.txt'))
ys3 = np.cumsum(np.loadtxt('result7.txt'))
ys4 = np.cumsum(np.loadtxt('result8.txt'))[1:]
xs = np.arange(0, 10, 1)

fig, ax = pyplot.subplots(figsize=(12, 8), dpi=400)

ax.minorticks_on()
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
ax.grid()
ax.grid(which='minor', linestyle=':')

ax.plot(xs, ys1, zorder=10, linewidth=2, label='Метод Якоби')
ax.plot(xs, ys2, zorder=10, linewidth=2, label='Метод Гаусса-Зейделя')
ax.plot(xs, ys3, zorder=10, linewidth=2, label='Метод простой итерации')
ax.plot(xs, ys4, zorder=10, linewidth=2, label='Метод простой итерации с ускорением')
ax.set_xlabel('Итерация', fontsize=14, fontname='Georgia')
ax.set_ylabel('Время', fontsize=14, fontname='Georgia')
ax.legend(loc='upper left')
pyplot.savefig("Graph 2.jpg", bbox_inches='tight')

