from matplotlib import pyplot as plt
import math

I_file = open("../data_backup/wall/#2/I_data.txt")
N_file = open("../data_backup/wall/#2/N_data.txt")
U_file = open("../data_backup/wall/#2/U_data.txt")


N_vals = []
N_err  = []
I_err  = []
U_err = []

for line in N_file:
    N_vals.append(float(line.split()[0]))
    N_err.append(float(line.split()[3]))

for line in I_file:
    I_err.append(float(line.split()[3]))

for line in U_file:
    U_err.append(float(line.split()[3]))

plt.yscale("log")
plt.xlabel("number of grid points $\\sqrt{N}$")
plt.ylabel("real error $|1 - \\frac{I_N}{I}|$")
plt.scatter(N_vals, N_err, s=2, alpha=0.7, label="weighted inportance sampling $X = \\Phi^{-1}(W)$", color="red")
plt.scatter(N_vals, I_err, s=2, alpha=0.7, label="uniform sampling $W$", color="black")
plt.scatter(N_vals, U_err, s=2, alpha=0.7, label="unweighted inportance sampling $X = \\Phi^{-1}(W)$", color="green")
plt.legend()
plt.show()
