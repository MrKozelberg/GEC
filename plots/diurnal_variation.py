import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

start = 0
final = 49
point1 = 18
point2 = 42

a = np.genfromtxt('diurnal_variation.txt', delimiter='\t')

s = 0
for i in range(start,final):
    s += a[i,1]
s = s / np.size(a[start:final,1])

fig, ax = plt.subplots()

ax.set_title("Diurnal Variation (2015-12-31)")
ax.set_xlabel(r'$t$  [hours]')
ax.set_ylabel(r'')
ax.set_xlim([0,48])
ax.set_ylim([0,1.5])

ax.plot(a[start:point1+1, 0], a[start:point1+1,1] / s,'--r',
        a[point1:point2+1, 0], a[point1:point2+1, 1] / s,'-r',
        a[point2:final, 0], a[point2:final,1] / s,'--r', linewidth = 2)


#  Устанавливаем интервал основных и
#  вспомогательных делений:
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))


#  Добавляем линии основной сетки:
ax.grid(which='major',
        color = 'gray')

#  Включаем видимость вспомогательных делений:
ax.minorticks_on()
#  Теперь можем отдельно задавать внешний вид
#  вспомогательной сетки:
ax.grid(which='minor',
        color = 'gray')

fig.set_figwidth(12)
fig.set_figheight(4)

plt.show()

