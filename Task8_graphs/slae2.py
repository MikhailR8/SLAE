import matplotlib.pyplot as pyplot
import numpy as np

n = 5
input1 = np.loadtxt('build/Conjugate_gradient.txt')
input2 = np.loadtxt('build/GMRES.txt')
ys1 = np.log(input1[0:input1.size//2])
xs1 = input1[input1.size//2:]
xs1 = np.cumsum(xs1)
xs2 = []
ys2 = []
s = input2.size // (n+1) // 2
print(s)
for i in range(s):
    # print(input2[(2*i*(n+1)) : (2*i*(n+1) + n)])
    ys2.append(np.log(np.fabs(input2[(2*i*(n+1)) : (2*i*(n+1) + n)])))
    # print(input2[(2*i*(n+1) + n) : (2*i*(n+1) + 2 * (n))])
    xs2.append(input2[(2*i*(n+1) + n) : (2*i*(n+1) + 2 * (n))])
    ys2.append(np.array([np.log(np.fabs(input2[(2*i*(n+1) + 2 * (n))]))]))
    xs2.append(np.array([input2[(2*i*(n+1) + 2 * (n))+1]]))

xs2 = np.concatenate(xs2, axis=0)
ys2 = np.concatenate(ys2, axis=0)
xs2 = np.cumsum(xs2)
print(xs2)
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
ax.text(300, 0, "Прямые участки GMRES - время решения СЛАУ гауссом")
pyplot.savefig("Graph 3.jpg", bbox_inches='tight')

