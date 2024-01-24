from matplotlib import pyplot as plt
import math

N_file = open("../data/q_factor.txt")

cycles = []
q_factor  = []

i = 1
for line in N_file:
    cycles.append(i)
    q_factor.append(float(line.split()[0]))
    i += 1


plt.xlabel("number of cycles$")
plt.ylabel("q-factor $q$")
plt.scatter(cycles, q_factor, s=2, alpha=0.5, color="black")
plt.show()
