from matplotlib import pyplot as plt

data_file = open("../data_backup/consecutive_runs/peak/consecutive_runs.txt", "r")

counter    = []
quality    = []
opt_vals   = []
int_vals   = []
n_int_vals = []
n_est_err  = []
n_real_err = []
u_int_vals = []
u_est_err  = []
u_real_err = []
i_int_vals = []
i_est_err  = []
i_real_err = []

center1 = []
center2 = []

i = 1
j = 1
for line in data_file:
    if (i % 2 == 0):
        int_vals.append(float(line.split()[0]))
        n_int_vals.append(float(line.split()[1]))
        n_est_err.append(float(line.split()[2]))
        n_real_err.append(float(line.split()[3]))
        u_int_vals.append(float(line.split()[4]))
        u_est_err.append(float(line.split()[5]))
        u_real_err.append(float(line.split()[6]))
        i_int_vals.append(float(line.split()[7]))
        i_est_err.append(float(line.split()[8]))
        i_real_err.append(float(line.split()[9]))
    else:
        counter.append(j)
        opt_vals.append(float(line.split()[3]))
        quality.append(float(line.split()[4]))
        center1.append(float(line.split()[1]))
        center2.append(float(line.split()[2]))
        j += 1
    i += 1
    

# plt.xlabel("number of runs")
# plt.ylabel("integral value $I_N$")
# plt.scatter(counter, n_int_vals, s=2, alpha=0.7, label="weighted inportance sampling $X = \\Phi^{-1}(W)$", color="red")
# plt.scatter(counter, u_int_vals, s=2, alpha=0.7, label="unweighted inportance sampling $X = \\Phi^{-1}(W)$", color="green")
# plt.scatter(counter, int_vals, s=2, alpha=0.7, label="real integral value", color="black")
# plt.legend()
# plt.show()

plt.yscale("log")
plt.xlabel("number of runs")
plt.ylabel("estimated error $\delta I_N$")
plt.scatter(counter, n_est_err, s=2, alpha=0.7, label="weighted inportance sampling $X = \\Phi^{-1}(W)$", color="red")
plt.scatter(counter, u_est_err, s=2, alpha=0.7, label="unweighted inportance sampling $X = \\Phi^{-1}(W)$", color="green")
plt.scatter(counter, i_est_err, s=2, alpha=0.7, label="uniform sampling $W$", color="black")
plt.legend()
plt.show()

# plt.yscale("log")
# plt.xlabel("number of runs")
# plt.ylabel("real error $|1 - \\frac{I_N}{I}|$")
# plt.scatter(counter, n_real_err, s=2, alpha=0.7, label="weighted inportance sampling $X = \\Phi^{-1}(W)$", color="red")
# plt.scatter(counter, u_real_err, s=2, alpha=0.7, label="unweighted inportance sampling $X = \\Phi^{-1}(W)$", color="green")
# plt.scatter(counter, i_real_err, s=2, alpha=0.7, label="uniform sampling $W$", color="black")
# plt.legend()
# plt.show()
