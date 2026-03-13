import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 100)
y = np.sin(x) + 0.5 * np.random.randn(100)
plt.figure(figsize=(8, 4))
plt.plot(x, y, "o", label="Data")

plt.show()
