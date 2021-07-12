import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

start = 18
final = 42
point1 = 18
point2 = 42

a = np.genfromtxt('diurnal_variation.txt', delimiter='\t')
p = np.polyfit(a[:,0], a[:,1],1)
b = np.polyval(p, a[:,0])

fig, ax = plt.subplots()

ax.set_title("Diurnal Variation (2015-12-31)")
ax.set_xlabel(r'$t$  [hours]')
ax.set_ylabel(r'IP as a fraction of the mean')
ax.set_xlim([18,42])
ax.set_ylim([0.95,1.05])

ax.plot(a[:,0], (a[:,1] - b + 240) / 240 , '-r', linewidth = 2)

#  Устанавливаем интервал основных и
#  вспомогательных делений:
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
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

fig.savefig('Diurnal Variation (2015-12-31) 18-42.png', dpi = 300)

plt.show()

