import numpy as np
import csv
import os
import matplotlib.pyplot as plt


N_directories = ['osc000']
x_dir_ibd = np.array([])
y_dir_ibd = np.array([])
z_dir_ibd = np.array([])
x_dir_es = np.array([])
y_dir_es = np.array([])
z_dir_es = np.array([])

def ratio_err(a,b,siga,sigb):
  return np.sqrt(float(siga*siga)/float(b*b) + float(sigb*sigb)*(float(a*a)/float(b**4)))

for num,d in enumerate(N_directories):

  for subdir, dirs, files in os.walk('/Users/noahsteinberg/Physics/James_Research/NSSI-Supernova/nu_events/'+ d):
    for file in files:
      filepath = subdir + os.sep + file
      if filepath.endswith("1.txt"):# and ('es' in filepath):
        with open(filepath) as inp:
          reader = csv.reader(inp, delimiter=" ")
          for row in reader:
            if len(row) == 6:

              if(filepath.endswith('ibd1.txt')):
                x_dir_ibd = np.append(x_dir_ibd,[float(row[2])])
                y_dir_ibd = np.append(y_dir_ibd,[float(row[3])])
                z_dir_ibd = np.append(z_dir_ibd,[float(row[4])])

              if(filepath.endswith('es1.txt')):
                x_dir_es = np.append(x_dir_es,[float(row[2])])
                y_dir_es = np.append(y_dir_es,[float(row[3])])
                z_dir_es = np.append(z_dir_es,[float(row[4])])


r_ibd =np.sqrt(x_dir_ibd**2 + y_dir_ibd**2 + z_dir_ibd**2)
phi_ibd = np.arctan2(y_dir_ibd,x_dir_ibd)
theta_ibd = np.arccos(np.divide(y_dir_ibd,r_ibd))


r_es =np.sqrt(x_dir_es**2 + y_dir_es**2 + z_dir_es**2)
phi_es = np.arctan(y_dir_es,x_dir_es)
theta_es = np.arccos(np.divide(y_dir_es,r_es))


theta_ibd,theta_bin_edges_ibd = np.histogram(theta_ibd,range = (0,np.pi),bins=20)
phi_ibd,phi_bin_edges_ibd = np.histogram(phi_ibd,range = (-np.pi,np.pi),bins=20)


theta_es,theta_bin_edges_es = np.histogram(theta_es,range = (0,np.pi),bins=20)
phi_es,phi_bin_edges_es = np.histogram(phi_es,range = (-np.pi,np.pi),bins=20)

fig = plt.figure()

ax2 = fig.add_subplot(121)
ax2.scatter(theta_bin_edges_es[:-1],theta_es,label = 'ES events')
ax2.scatter(theta_bin_edges_ibd[:-1],theta_ibd, label='IBD events')

ax1 = fig.add_subplot(122)
ax1.scatter(phi_bin_edges_es[:-1],phi_es, label='ES events')
ax1.scatter(phi_bin_edges_ibd[:-1],phi_ibd, label='IBD events')

ax1.set_ylabel('Events')
ax1.set_xlabel('phi (radians)')
ax1.set_title('Phi')
ax1.set_yscale('log')
ax1.set_ylim(1)
ax1.legend(loc='upper right');

ax2.set_ylabel('Events')
ax2.set_xlabel('theta (radians)')
ax2.set_title('Theta')
ax2.set_yscale('log')
ax2.set_ylim(1)
ax2.legend(loc='upper right');

plt.show()
