import numpy as np
from matplotlib import pyplot as plt
import os
import os.path

plt.ion();

net_name = "email.edgelist";

lambd = 1.0;
alpha = 1.0;
delta_1 = 0.1;
delta_2 = 0.1;
delta_3 = 0.1;
beta = 1.0;
eta = 0.01;
Tf = 100;

cmd = "./dtmc " +  net_name + " " + str(lambd) + " " + str(alpha) + " " + str(delta_1) + " " + str(delta_2) + " " + str(delta_3) + " " + str(eta) + " " + str(beta) + " " + str(Tf) + " cp both";
out_mc = "%s_MicroFinal_RumorMC_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat" % (net_name, lambd, alpha, delta_1, delta_2, delta_3, eta, beta, Tf);
out_mk = "%s_Final_Nodal_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat" % (net_name, lambd, alpha, delta_1, delta_2, delta_3, eta, beta, Tf);


if os.path.isfile(out_mk) == False:
	os.system(cmd);

DataMC = np.loadtxt(out_mc);
DataMK = np.loadtxt(out_mk);

fig = plt.figure();
#plt.add_subplot(1,1,1, rasterized=True)

plt.scatter(DataMK[:,0], DataMC[:,0], s=50, alpha=.35, marker='o', c ='#c92f2f', zorder=3, label="Ignorant");
plt.scatter(DataMK[:,1], DataMC[:,1], s=50, alpha=.35, marker='s', c ='#15c8f9', zorder=2, label="Spreader");
plt.scatter(DataMK[:,2], DataMC[:,2], s=50, alpha=.35, marker='d', c ='#42269b', zorder=1, label="Stifler");


plt.ylim([0, 1])
plt.xlim([0, 1])
plt.legend(loc='lower right')

plt.rcParams.update({'font.size': 18})
plt.xlabel('Analytical', fontsize=18)
plt.ylabel('MC', fontsize=18)

plt.tight_layout();
plt.savefig("email_MK_CP_Micro.pdf")

plt.show()
