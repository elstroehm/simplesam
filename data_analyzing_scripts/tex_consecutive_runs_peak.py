#from matplotlib import pyplot as plt
#import numpy as np
#import pandas as pd

# data_file = open("../data/consecutive_runs.txt")
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
# height  = []

i = 1
j = 1
for line in data_file:
    if (i % 2 == 0):
        int_vals.append(round(float(line.split()[0]), 2))
        n_int_vals.append(round(float(line.split()[1]), 2))
        n_est_err.append(round(float(line.split()[2]), 2))
        n_real_err.append(round(float(line.split()[3]), 5))
        u_int_vals.append(round(float(line.split()[4]), 2))
        u_est_err.append(round(float(line.split()[5]), 2))
        u_real_err.append(round(float(line.split()[6]), 5))
        i_int_vals.append(round(float(line.split()[7]), 2))
        i_est_err.append(round(float(line.split()[8]), 2))
        i_real_err.append(round(float(line.split()[9]), 5))
    else:
        counter.append(j)
        opt_vals.append(round(float(line.split()[3]), 3))
        quality.append(round(float(line.split()[4]), 3))
        center1.append(round(float(line.split()[1]), 2))
        center2.append(round(float(line.split()[2]), 2))
        # height.append(float(line.split()[1]))
        j += 1
    i += 1

i_arr = []
n_arr = []
u_arr = []
center_arr = []

j = 0
for l in counter:
    i_arr.append(str(i_int_vals[j]) + " $\\pm$ " + str(i_est_err[j]))
    n_arr.append(str(n_int_vals[j]) + " $\\pm$ " + str(n_est_err[j]))
    u_arr.append(str(u_int_vals[j]) + " $\\pm$ " + str(u_est_err[j]))
    center_arr.append("$(" + str(center1[j]) + " , " + str(center2[j]) + ")$")
    j += 1
    
latex_array_1 = []
latex_array_2 = []

for k in range(0, len(counter), 1):
    latex_array_1.append( str(counter[k]) + " & " + center_arr[k] + " & " + str(int_vals[k]) + " & "
    + n_arr[k] + " & " + u_arr[k] + " & " + i_arr[k] + "\\\\ \n \\hline \n" )
    latex_array_2.append( str(counter[k]) + " & " + center_arr[k] + " & " + str(quality[k]) + " & " + str(opt_vals[k]) + " & " + str(n_real_err[k]) + " & " + str(u_real_err[k]) + " & " + str(i_real_err[k])
    + "\\\\ \n \\hline \n" )

tex_file_1 = open("../data_backup/consecutive_runs/peak/tex_data_1.txt", "w")
tex_file_2 = open("../data_backup/consecutive_runs/peak/tex_data_2.txt", "w")

for k in range(0, len(latex_array_1), 1):
    tex_file_1.write(latex_array_1[k])
    tex_file_2.write(latex_array_2[k])
    

# plt.xlabel("number of runs")
# plt.ylabel("integral value $I_N$")
# plt.scatter(counter, n_int_vals, s=2, alpha=0.7, label="weighted inportance sampling $X = \\Phi^{-1}(W)$", color="red")
# plt.scatter(counter, u_int_vals, s=2, alpha=0.7, label="unweighted inportance sampling $X = \\Phi^{-1}(W)$", color="green")
# plt.scatter(counter, int_vals, s=2, alpha=0.7, label="real integral value", color="black")
# plt.legend()
# plt.show()

# plt.yscale("log")
# plt.xlabel("number of runs")
# plt.ylabel("estimated error $\delta I_N$")
# plt.scatter(counter, n_est_err, s=2, alpha=0.7, label="weighted inportance sampling $X = \\Phi^{-1}(W)$", color="red")
# plt.scatter(counter, u_est_err, s=2, alpha=0.7, label="unweighted inportance sampling $X = \\Phi^{-1}(W)$", color="green")
# plt.legend()
# plt.show()

# plt.yscale("log")
# plt.xlabel("number of runs")
# plt.ylabel("real error $|1 - \\frac{I_N}{I}|$")
# plt.scatter(counter, n_real_err, s=2, alpha=0.7, label="weighted inportance sampling $X = \\Phi^{-1}(W)$", color="red")
# plt.scatter(counter, u_real_err, s=2, alpha=0.7, label="unweighted inportance sampling $X = \\Phi^{-1}(W)$", color="green")
# plt.legend()
# plt.show()
