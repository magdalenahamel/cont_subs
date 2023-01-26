import numpy as np
import ppxf
import matplotlib.pyplot as plt
plt.switch_backend('agg')
a = np.linspace(0,10, 10)
b = np.random.random(10)
print('a',a)
print('b', b)
plt.savefig('example')
