import matplotlib.pyplot as pyplot
import numpy as np

ys1 = np.loadtxt('build/result4.txt')
ys2 = np.loadtxt('build/result5.txt')
ys3 = np.loadtxt('build/result6.txt')
xs = np.arange(1000, 2000, 20)

fig, ax = pyplot.subplots(figsize=(12, 8), dpi=400)

ax.minorticks_on()
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
ax.grid()
ax.grid(which='minor', linestyle=':')

ax.plot(xs, ys1, zorder=10, linewidth=2, label='Метод Якоби')
ax.plot(xs, ys2, zorder=10, linewidth=2, label='Метод Гаусса-Зейделя')
ax.plot(xs, ys3, zorder=10, linewidth=2, label='Матрица простой итерации')
ax.set_xlabel('Размер одной стороны квадратной матрицы со строгим диагональным преобладанием', fontsize=14, fontname='Georgia')
ax.set_ylabel('Число итераций (max 100)', fontsize=14, fontname='Georgia')
ax.legend(loc='upper left')
pyplot.savefig("Graph 2.jpg", bbox_inches='tight')

