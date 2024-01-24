from matplotlib import pyplot as plt
import math

I_file = open("../data/I_data.txt")
N_file = open("../data/N_data.txt")
U_file = open("../data/U_data.txt")


N_vals = []
N_err  = []
I_err  = []
U_err = []

for line in N_file:
    N_vals.append(float(line.split()[0]))
    N_err.append(float(line.split()[1]))

for line in I_file:
    I_err.append(float(line.split()[1]))

for line in U_file:
    U_err.append(float(line.split()[1]))

plt.xlabel("number of grid points $\\sqrt{N}$")
plt.ylabel("integral value $I_N$")
plt.scatter(N_vals, N_err, s=2, alpha=0.7, label="weighted inportance sampling $X = \\Phi^{-1}(W)$", color="red")
plt.scatter(N_vals, I_err, s=2, alpha=0.7, label="uniform sampling $W$", color="black")
plt.scatter(N_vals, U_err, s=2, alpha=0.7, label="unweighted inportance sampling $X = \\Phi^{-1}(W)$", color="green")
plt.legend()
plt.show()
