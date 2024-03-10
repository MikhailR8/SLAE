import matplotlib.pyplot as pyplot
import numpy as np

ys1 = np.log(np.loadtxt('result1.txt'))[0:10]
ys2 = np.log(np.loadtxt('result2.txt'))[0:10]
ys3 = np.log(np.loadtxt('result3.txt'))[0:10]
ys4 = np.log(np.loadtxt('result4.txt'))[0:10]
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
ax.set_ylabel('Логарифм неувязки', fontsize=14, fontname='Georgia')
ax.legend(loc='upper left')
pyplot.savefig("Graph 1.jpg", bbox_inches='tight')

