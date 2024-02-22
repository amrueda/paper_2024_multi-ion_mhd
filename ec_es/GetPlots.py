#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 15:21:40 2020

@author: andresrueda
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches



CFL=np.array([0.025,0.05,0.1,0.2,0.4])

folderNames=['CFL_0.025','CFL_0.050','CFL_0.100','CFL_0.200','CFL_0.400']
folderN=['','']

# READ DATA

dS_EC=np.zeros([2,5])
dS_ES=np.zeros([2,5])
dS_llf=np.zeros([2,5])

for i in range(5):
    print(i)
    data=np.genfromtxt("./out_ec_cfl"+str(CFL[i])+"/analysis.dat", names=True,invalid_raise=False)
    dS_EC[0][i]=(data['entropy_density'][1] - data['entropy_density'][-1])
    dS_EC[1][i]=max(data['dsdu_ut'][1:-1])
    print(len(data['entropy_density']))
        
    data=np.genfromtxt("./out_es_cfl"+str(CFL[i])+"/analysis.dat", names=True,invalid_raise=False)
    dS_ES[0][i]=(data['entropy_density'][1] - data['entropy_density'][-1])
    dS_ES[1][i]=max(data['dsdu_ut'][1:-1])
    print(len(data['entropy_density']))
    
    data=np.genfromtxt("./out_llf_cfl"+str(CFL[i])+"/analysis.dat", names=True,invalid_raise=False)
    dS_llf[0][i]=(data['entropy_density'][1] - data['entropy_density'][-1])
    dS_llf[1][i]=max(data['dsdu_ut'][1:-1])
    print(len(data['entropy_density']))
    
#
fig = plt.figure()
plt.loglog(CFL,dS_EC[0],'.-',fillstyle='none',linewidth=1, label="EC")
plt.loglog(CFL,dS_ES[0],'s--',fillstyle='none',linewidth=1, label="ES")
plt.loglog(CFL,dS_llf[0],'^:',fillstyle='none',linewidth=1, label="EC+LLF")

x1 = [0.1, 0.2]
y1 = [1e-12, 1e-12]
plt.plot(x1, y1,'k-')
x1 = [0.1, 0.2]
y1 = [1e-12, 1.6e-11]
plt.plot(x1, y1,'k-')
x1 = [0.2, 0.2]
y1 = [1e-12, 1.6e-11]
plt.plot(x1, y1,'k-')
plt.text(0.21,3e-12,'4', fontsize=12)
plt.text(0.14,2e-13,'1', fontsize=12)

plt.legend()
#plt.legend(["$N=4$, EC","$N=4$, ES","$N=6$, EC","$N=6$, ES"],loc='center left')
plt.xticks(CFL,labels=CFL)

plt.xlabel('CFL', fontsize=12)
plt.ylabel('$S_{\Omega}(t=0) - S_{\Omega}(t=0.4)$', fontsize=12) #frac{\eta}{\eta_0}
#

# ######################
# # Show a detail
# ax = plt.gca()
# x1 =0.024
# x2 =0.401
# y1 =9.5e-5
# y2 =9.59e-5
# rect = patches.Rectangle((x1,y1), x2-x1,y2-y1, linewidth=1, edgecolor='k', facecolor='none')
# ax.add_patch(rect)

# plt.plot([0.025,x1],[1e-8,y1],'k--', linewidth=1)
# plt.plot([0.050,x2],[1e-12,y2],'k--', linewidth=1)
# ax_new = fig.add_axes([0.15, 0.4, 0.3, 0.25]) 
# plt.loglog(CFL,dS_ES[0],'s--',fillstyle='none',linewidth=1, label="ES")
# plt.loglog(CFL,dS_llf[0],'^:',fillstyle='none',linewidth=1, label="EC+LLF")
# plt.xlim(x1,x2)
# plt.ylim(y1,y2)
# plt.xticks([]),plt.yticks([])
# plt.tick_params(left = False) 
# ########################

plt.show()

plt.figure()
plt.loglog(CFL,abs(dS_EC[1]),'.-',fillstyle='none',linewidth=1, label="EC ($|\min(-\dot{S}_{\Omega})|$)")
plt.loglog(CFL,-dS_ES[1],'s--',fillstyle='none',linewidth=1, label="ES")
plt.loglog(CFL,-dS_llf[1],'^:',fillstyle='none',linewidth=1, label="EC+LLF")

plt.xlabel('CFL', fontsize=12)
plt.ylabel('$\min (-\dot{S}_{\Omega})$', fontsize=12) #frac{\eta}{\eta_0}



plt.legend()
plt.xticks(CFL,labels=CFL)

