import numpy as np
from matplotlib import pyplot as plt
import os
import os.path

import matplotlib.patches as patches

plt.ion();

fig = plt.figure();
ax = fig.add_subplot(111);


Org = np.loadtxt("Data_time_Higgs.txt");

ax.plot(Org[:,0], Org[:,1], 'k', label="Real spreaders", lw=2, ls='--')
ax.plot(Org[:,0], Org[:,2], 'k', label="Real stiflers", lw=2, ls='-.')

# Part 1

lambd = 0.00025;
alpha = 0.0002;
delta_1 = 0.0;
delta_2 = 0.0001;
delta_3 = 0.0;
beta = 0.0;
eta = 0.85;
Tf = 78;

#out_nodal = "Nodes_InitialCondition.txt"
#cmd = "./dtmc higgs-social_network_undirected_full_0.edgelist " + str(lambd) + " " + str(alpha) + " " + str(delta_1) + " " + str(delta_2) + " " + str(delta_3) + " " + str(eta) + " " + str(beta) + " " + str(Tf) + " " + out_nodal + " rp markov";

cmd = "./dtmc higgs-social_network_undirected_full_0.edgelist " + str(lambd) + " " + str(alpha) + " " + str(delta_1) + " " + str(delta_2) + " " + str(delta_3) + " " + str(eta) + " " + str(beta) + " " + str(Tf) + " " + "5.543315915571e-02 " + "2.422791718314e-02" + " rp markov";
out = "higgs-social_network_undirected_full_0.edgelist_Time_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat" % (lambd, alpha, delta_1, delta_2, delta_3, eta, beta, Tf);
out_nodal = "higgs-social_network_undirected_full_0.edgelist_Final_Nodal_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat" % (lambd, alpha, delta_1, delta_2, delta_3, eta, beta, Tf);


Ex0Bool = False;

if os.path.isfile(out) == False:
	os.system(cmd);
	Ex0Bool = True;

Exp0 = np.loadtxt(out);

ax.plot(Exp0[:,0], '-r')
ax.plot(Exp0[:,1], '-b')
ax.plot(Exp0[:,2], '-g')
#ax.plot(Exp0[:,2] + Exp0[:,1], 'darkred')

# Part 2
lambd = 0.021;
alpha = 0.00075;
delta_1 = 0.0;
delta_2 = 0.0015;
delta_3 = 0.0;
beta = 0.00;
eta = 0.17;
Tf2 = 30;

cmd = "./dtmc higgs-social_network_undirected_full_0.edgelist " + str(lambd) + " " + str(alpha) + " " + str(delta_1) + " " + str(delta_2) + " " + str(delta_3) + " " + str(eta) + " " + str(beta) + " " + str(Tf2) + " " + out_nodal + " rp markov";
out = "higgs-social_network_undirected_full_0.edgelist_Time_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat" % (lambd, alpha, delta_1, delta_2, delta_3, eta, beta, Tf2);
out_nodal = "higgs-social_network_undirected_full_0.edgelist_Final_Nodal_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat" % (lambd, alpha, delta_1, delta_2, delta_3, eta, beta, Tf2);


if os.path.isfile(out) == False or Ex0Bool == True:
	os.system(cmd);

Exp1 = np.loadtxt(out);

ax.plot(np.array(range(Tf2))+(Tf-1), Exp1[:,0], '-r')
ax.plot(np.array(range(Tf2))+(Tf-1), Exp1[:,1], '-b')
ax.plot(np.array(range(Tf2))+(Tf-1), Exp1[:,2], '-g')
#ax.plot(np.array(range(Tf2))+(Tf-1), Exp1[:,2] + Exp1[:,1], 'darkred', label="Spreader or Stifler")


# Part 3

lambd = 0.065;
alpha = 0.002;
delta_1 = 0.0;
delta_2 = 0.002;
delta_3 = 0.0;
beta = 0.00;
eta = 0.01;
Tf3 = 75;

cmd = "./dtmc higgs-social_network_undirected_full_0.edgelist " + str(lambd) + " " + str(alpha) + " " + str(delta_1) + " " + str(delta_2) + " " + str(delta_3) + " " + str(eta) + " " + str(beta) + " " + str(Tf3) + " " + out_nodal + " rp markov";
out = "higgs-social_network_undirected_full_0.edgelist_Time_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat" % (lambd, alpha, delta_1, delta_2, delta_3, eta, beta, Tf3);
out_nodal = "higgs-social_network_undirected_full_0.edgelist_Final_Nodal_RumorMarkov_lambda_%.6f_alpha_%.6f_delta1_%.6f_delta2_%.6f_gamma_%.6f_eta_%.6f_beta_%.6f_tmax_%ld.dat" % (lambd, alpha, delta_1, delta_2, delta_3, eta, beta, Tf3);

if os.path.isfile(out) == False or Ex0Bool == True:
	os.system(cmd);

Exp1 = np.loadtxt(out);

ax.plot(np.array(range(Tf3))+(Tf+Tf2-2), Exp1[:,0], '-r', label="Ignorant")
ax.plot(np.array(range(Tf3))+(Tf+Tf2-2), Exp1[:,1], '-b', label="Spreader")
ax.plot(np.array(range(Tf3))+(Tf+Tf2-2), Exp1[:,2], '-g', label="Stifler")
#ax.plot(np.array(range(Tf2))+(Tf-1), Exp1[:,2] + Exp1[:,1], 'darkred', label="Spreader or Stifler")


#plt.xlim([0, Tf+Tf2])
ax.set_xlim([0, 168])
ax.set_ylim([0, np.max(np.max(Exp0))])


p1 = patches.Rectangle((0, 0), Tf-1, np.max(np.max(Exp0)), alpha=0.1, edgecolor="none");
ax.add_patch(p1)
p2 = patches.Rectangle((Tf-1, 0), Tf2-1, np.max(np.max(Exp0)), alpha=0.1, facecolor="red", edgecolor="none");
ax.add_patch(p2)
p3 = patches.Rectangle((Tf2+Tf-2, 0), Tf3, np.max(np.max(Exp0)), alpha=0.1, facecolor="blue", edgecolor="none");
ax.add_patch(p3)

ax.legend(loc='center left')

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.rcParams.update({'font.size': 12})
plt.xlabel('Time (hours)', fontsize=18)
plt.ylabel('Number of individuals', fontsize=18)

plt.tight_layout();

plt.savefig("DTMC_Higgs.pdf")

plt.show()
