import matplotlib.pyplot as pyplot
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker

name = 'ass2.txt'
ys1 = np.loadtxt('result1.txt')
ys2 = np.loadtxt('result2.txt')
ys3 = np.loadtxt('result3.txt')
ys4 = np.loadtxt('result4.txt')
xs = np.arange(1, 1505, 5)

fig, ax = pyplot.subplots(figsize=(12, 8), dpi=400)

ax.minorticks_on()
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
ax.grid()
ax.grid(which='minor', linestyle=':')

ax.plot(xs, ys1, zorder=10, linewidth=2, label='Плотная матрица в обычном виде')
ax.plot(xs, ys2, zorder=10, linewidth=2, label='Матрица с нулями в обычном виде')
ax.plot(xs, ys3, zorder=10, linewidth=2, label='Плотная матрица CSR')
ax.plot(xs, ys4, zorder=10, linewidth=2, label='Матрица с нулями CSR')
ax.set_xlabel('Размер одной стороны квадратной матрицы', fontsize=14, fontname='Georgia')
ax.set_ylabel('Время в микросекундах', fontsize=14, fontname='Georgia')
ax.legend(loc='lower right')
pyplot.savefig("Graph 1.jpg", bbox_inches='tight')

pyplot.show()
