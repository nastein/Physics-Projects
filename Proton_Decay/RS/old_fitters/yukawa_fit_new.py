import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.optimize
from scipy.optimize import fsolve
import sys
import csv
from numpy import sqrt as sqrt
from numpy import log as log
pi = np.pi

#Kr as a function of Lambda_IR
def kr(lam):
	return np.log((2.435e18)/lam)/pi

#quark yukawa coupling as a function of cl,cr, and Lambda_IR
def qyuk(cl,cr,lam):
	return sqrt((np.exp(2*(1 - cl - cr)*pi*kr(lam))*(.5 - cl)*(.5 - cr))/((np.exp((1 - 2*cl)*pi*kr(lam)) - 1)*(np.exp((1 - 2*cr)*pi*kr(lam)) - 1)))

mu = []
md = []
ms = []
mc = []
mb = []
mt = []
me = []
mmu = []
mtau = []

lambdas = np.logspace(3,17,60)
lambdas = [lambdas[0],lambdas[1],lambdas[2]]

i=0
with open('quark_masses.csv','r') as f:
	reader = csv.reader(f)
	for row in reader:
		if i == 0:
			mu = np.asfarray(np.array(row),float)
		if i == 1:
			md = np.asfarray(np.array(row),float)
		if i == 2:
			ms = np.asfarray(np.array(row),float)
		if i == 3:
			mc = np.asfarray(np.array(row),float)
		if i == 4:
			mb = np.asfarray(np.array(row),float)
		if i == 5:
			mt = np.asfarray(np.array(row),float)
		i = i + 1




j=0
with open('lepton_masses.csv','r') as h:
	reader = csv.reader(h)
	for row in reader:
		if j == 0:
			me = np.asfarray(np.array(row),float)
		if j == 1:
			mmu = np.asfarray(np.array(row),float)
		if j == 2:
			mtau = np.asfarray(np.array(row),float)
		j = j+1

yu = []
yd = []
ye = []

for i in range(0,len(lambdas)):
	yu.append([[(sqrt(2)*mu[i])/246, 0., 0.],[0., (sqrt(2)*mc[i])/246, 0.],[0., 0., (sqrt(2)*mt[i])/246]])
	yd.append([[(sqrt(2)*md[i])/246, 0., 0.],[0., (sqrt(2)*ms[i])/246, 0.],[0., 0., (sqrt(2)*mb[i])/246]])
	ye.append([[me[i], 0., 0.],[0., mmu[i], 0.],[0., 0., mtau[i]]])

Vckm = [[.97427, .22534, .00351],[.22520, .97344, .0413],[.00867, .0404, .999146]]



def cal_c(lam,yu,yd,ye):

	def f(c):
		cq1,cq2,cq3,cu1,cu2,cu3,cd1,cd2,cd3,ce1,ce2,ce3,lu1,lu2,lu3,lu4,lu5,lu6,lu7,lu8,lu9,du1,du2,du3,du4,du5,du6,du7,du8,du9 = c


		Gu = [[qyuk(cq1,cu1,lam),qyuk(cq1,cu2,lam),qyuk(cq1,cu3,lam)],
	  		  [qyuk(cq2,cu1,lam),qyuk(cq2,cu2,lam),qyuk(cq2,cu3,lam)],
	  		  [qyuk(cq3,cu1,lam),qyuk(cq3,cu2,lam),qyuk(cq3,cu3,lam)]]

		Gd = [[qyuk(cq1,cd1,lam),qyuk(cq1,cd2,lam),qyuk(cq1,cd3,lam)],
	  		  [qyuk(cq2,cd1,lam),qyuk(cq2,cd2,lam),qyuk(cq2,cd3,lam)],
	  		  [qyuk(cq3,cd1,lam),qyuk(cq3,cd2,lam),qyuk(cq3,cd3,lam)]]

		Gu_squared = np.matmul(Gu,np.transpose(Gu))
		Gd_squared = np.matmul(Gd,np.transpose(Gd))

		Ul = [[lu1,lu2,lu3],[lu4,lu5,lu6],[lu7,lu8,lu9]]
		Dl = [[du1,du2,du3],[du4,du5,du6],[du7,du8,du9]]

		rul = np.matmul(Ul,np.matmul(yu,np.matmul(yu,np.transpose(Ul))))
		rdl = np.matmul(Dl,np.matmul(yd,np.matmul(yd,np.transpose(Dl))))

		rvckm = np.matmul(np.transpose(Ul),Dl)


		return (
			Gu_squared[0][0] - rul[0][0],Gu_squared[0][1] - rul[0][1],Gu_squared[0][2] - rul[0][2],
			Gu_squared[1][0] - rul[1][0],Gu_squared[1][1] - rul[1][1],Gu_squared[1][2] - rul[1][2],
			Gu_squared[2][0] - rul[2][0],Gu_squared[2][1] - rul[2][1],Gu_squared[2][2] - rul[2][2],
			Gd_squared[0][0] - rdl[0][0],Gd_squared[0][1] - rdl[0][1],Gd_squared[0][2] - rdl[0][2],
			Gd_squared[1][0] - rdl[1][0],Gd_squared[1][1] - rdl[1][1],Gd_squared[1][2] - rdl[1][2],
			Gd_squared[2][0] - rdl[2][0],Gd_squared[2][1] - rdl[2][1],Gd_squared[2][2] - rdl[2][2],
			qyuk(ce1,ce1,lam) - ye[0][0],qyuk(ce2,ce2,lam) - ye[1][1],qyuk(ce3,ce3,lam) - ye[2][2],
			Vckm[0][0] - rvckm[0][0],Vckm[0][1] - rvckm[0][1],Vckm[0][2] - rvckm[0][2],
			Vckm[1][0] - rvckm[1][0],Vckm[1][1] - rvckm[1][1],Vckm[1][2] - rvckm[1][2],
			Vckm[2][0] - rvckm[2][0],Vckm[2][1] - rvckm[2][1],Vckm[2][2] - rvckm[2][2],
			)
	
	#Construct an initial guess
	c0 = np.empty(30)
	c0.fill(.7)

	#Get the solution for the system
	c = scipy.optimize.least_squares(f,c0,method='lm').x

	cq = [c[0],c[1],c[2]]
	cu = [c[3],c[4],c[5]]
	cd = [c[6],c[7],c[8]]
	ce = [c[9],c[10],c[11]]
	
	#Gather predictions for yukawa matrices and CKM matrix
	ul_predict = [[c[12],c[13],c[14]],[c[15],c[16],c[17]],[c[18],c[19],c[20]]]
	dl_predict = [[c[21],c[22],c[23]],[c[24],c[25],c[26]],[c[27],c[28],c[29]]]

	#Predicted CKM matrix
	vckm_predict = np.matmul(np.transpose(ul_predict),dl_predict)

	#Predicted undiagonalized yukawa couplings
	Gu = [[qyuk(cq[0],cu[0],lam),qyuk(cq[0],cu[1],lam),qyuk(c[0],cu[2],lam)],
	  	  [qyuk(cq[1],cu[0],lam),qyuk(cq[1],cu[1],lam),qyuk(c[1],cu[2],lam)],
	  	  [qyuk(cq[2],cu[0],lam),qyuk(cq[2],cu[1],lam),qyuk(c[2],cu[2],lam)]]

	Gd = [[qyuk(cq[0],cu[0],lam),qyuk(cq[0],cu[1],lam),qyuk(cq[0],cu[2],lam)],
	      [qyuk(cq[1],cu[0],lam),qyuk(cq[1],cu[1],lam),qyuk(cq[1],cu[2],lam)],
	  	  [qyuk(cq[2],cu[0],lam),qyuk(cq[2],cu[1],lam),qyuk(cq[2],cu[2],lam)]]

	#Predicted undiagonalized yukawa couplings squared
	yusq_predict = np.matmul(np.transpose(ul_predict),np.matmul(Gu,np.matmul(np.transpose(Gu),ul_predict)))
	ydsq_predict = np.matmul(np.transpose(dl_predict),np.matmul(Gd,np.matmul(np.transpose(Gd),dl_predict)))

	return cq,cu,cd,ce,vckm_predict,yusq_predict,ydsq_predict

elec_c = []
mu_c = []
tau_c = []
up_c = []
down_c = []
strange_c = []
charm_c = []
bottom_c = []
top_c = []
left_up_c = []
left_charm_c = []
left_top_c = []

ckm_predict_ = []
yusq_ = []
ydsq_ = []
yusq_predict_ = []
ydsq_predict_ = []

pred_up_yukawa = []
pred_strange_yukawa = []
pred_top_yukawa = []

up_yukawa = []
strange_yukawa = []
top_yukawa = []

for i in range(0,3):
	c = cal_c(lambdas[i],yu[i],yd[i],ye[i])
	elec_c.append(c[3][0])
	mu_c.append(c[3][1])
	tau_c.append(c[3][2])
	up_c.append(c[1][0])
	down_c.append(c[2][0])
	strange_c.append(c[2][1])
	charm_c.append(c[1][1])
	bottom_c.append(c[2][2])
	top_c.append(c[1][2])
	left_up_c.append(c[0][0])
	left_charm_c.append(c[0][1])
	left_top_c.append(c[0][2])
	ckm_predict_.append(c[4])
	yusq_predict_.append(c[5])
	ydsq_predict_.append(c[6])
	yusq_.append(np.matmul(yu[i],yu[i]))
	ydsq_.append(np.matmul(yd[i],yd[i]))
	pred_up_yukawa.append(sqrt(yusq_predict_[i][0][0]))
	pred_strange_yukawa.append(sqrt(ydsq_predict_[i][1][1]))
	pred_top_yukawa.append(sqrt(yusq_predict_[i][2][2]))


fig = plt.figure()

ax1 = fig.add_subplot(311)
ax1.set_xscale('log')
#ax1.scatter(lambdas,elec_c,c='red',label='electron c')
#ax1.scatter(lambdas,mu_c,c='brown',label='muon c')
#ax1.scatter(lambdas,tau_c,c='black',label='tau c')
ax1.scatter(lambdas,left_up_c,c='green', label='left up c')
ax1.scatter(lambdas,left_charm_c,c='purple',label='left charm c')
ax1.scatter(lambdas,left_top_c,c='blue',label='left top c')
ax1.legend(loc=3)
ax1.set_title('5d mass parameters for left handed quarks and leptons')
#ax1.set_ylim([-2,5])

ax2 = fig.add_subplot(312)
ax2.set_xscale('log')
ax2.scatter(lambdas,up_c,c='red', label='right up c')
ax2.scatter(lambdas,down_c,c='brown',label='right down c')
ax2.scatter(lambdas,strange_c,c='black',label='right strange c')
ax2.scatter(lambdas,charm_c,c='green',label='right charm c')
ax2.scatter(lambdas,bottom_c,c='purple',label='right bottom c')
ax2.scatter(lambdas,top_c,c='blue',label='right top c')
ax2.legend(loc=3)
ax2.set_title('5d mass parameters for right handed quarks')
#ax2.set_ylim([-10,5])

ax3 = fig.add_subplot(313)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.scatter(lambdas, pred_up_yukawa,c='red',label='up quark yukawa')
ax3.scatter(lambdas, pred_strange_yukawa,c='black',label='strange quark yukawa')
ax3.scatter(lambdas, pred_top_yukawa,c='blue',label='top quark yukawa')
ax3.legend(loc=3)
ax3.set_title('Predicted quark yukawa couplings')


plt.show()






