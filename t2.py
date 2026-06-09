import matplotlib.pyplot as plt
import numpy as np

a = np.arange(64) * 8
b = np.arange(40) * 6 + 136
c = np.arange(20) * 4 + 216

arr = np.concatenate([a, b, b, b, c, c, c, c, c])

arrs = [arr, arr + 50, arr + 100, arr + 150]

parts = plt.violinplot(arrs, np.arange(1, 5), showmeans=True)
for pc in parts['bodies']:
    pc.set_facecolor('#D43F3A')
    pc.set_edgecolor('black')
    pc.set_alpha(1)
plt.show()
plt.violinplot(arrs, np.arange(1, 5), showmedians=True)
plt.show()