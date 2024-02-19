import matplotlib.pyplot as pyplot
import numpy as np

ys1 = np.loadtxt('build/result1.txt')
ys2 = np.loadtxt('build/result2.txt')
ys3 = np.loadtxt('build/result3.txt')
ys4 = np.loadtxt('build/result4.txt')
xs = np.arange(100, 2000, 10)

fig, ax = pyplot.subplots(figsize=(12, 8), dpi=400)

ax.minorticks_on()
# ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
ax.grid()
ax.grid(which='minor', linestyle=':')

ax.plot(xs, ys1, zorder=10, linewidth=2, label='Плотная матрица, хранимая в обычном виде')
ax.plot(xs, ys2, zorder=10, linewidth=2, label='Плотная матрица, хранимая в виде CSR')
ax.plot(xs, ys3, zorder=10, linewidth=2, label='Матрица с нулями, хранимая в обычном виде')
ax.plot(xs, ys4, zorder=10, linewidth=2, label='Матрица с нулями, хранимая в виде CSR')
ax.set_xlabel('Размер одной стороны квадратной матрицы', fontsize=14, fontname='Georgia')
ax.set_ylabel('Время в микросекундах', fontsize=14, fontname='Georgia')
ax.legend(loc='upper left')
pyplot.savefig("Graph 1.jpg", bbox_inches='tight')

