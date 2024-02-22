#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:25:46 2023

@author: andresrueda
"""
%matplotlib qt
import numpy as np
import matplotlib.pyplot as plt

khi_glm_H_H2 = np.genfromtxt("khi_glm_H_H2_ecllf/out/analysis.dat", names=True, invalid_raise=False) # , delimiter="," #,skip_header=2
khi_glm_H_H2_es = np.genfromtxt("khi_glm_H_H2_es/out/analysis.dat", names=True, invalid_raise=False) # , delimiter="," #,skip_header=2
khi_noglm_H_H2 = np.genfromtxt("khi_noglm_H_H2_ecllf/out/analysis.dat", names=True, invalid_raise=False) # , delimiter="," #,skip_header=2
khi_noglm_H_H2_es = np.genfromtxt("khi_noglm_H_H2_es/out/analysis.dat", names=True, invalid_raise=False) # , delimiter="," #,skip_header=2
khi_glm_H_H2_central = np.genfromtxt("khi_glm_H_H2_std/out/analysis.dat", names=True, invalid_raise=False) # , delimiter="," #,skip_header=2

def plot_poloidal(data, linestyl = "-", label = "", **kwargs):
    plt.plot(data['time'], data['poloidal_field'] / data['poloidal_field'][0], linestyl, label = label, **kwargs)

def plot_divb_l2(data, linestyl = "-", label = "", **kwargs):
    plt.plot(data['time'], data['l2_divb'], linestyl, label = label, **kwargs)

#########################

# Plot poloidal field for multi-species case
plot_poloidal(khi_glm_H_H2_es, 'b-', "split DG ES")
plot_poloidal(khi_glm_H_H2, 'k:', "split DG EC+LLF")
plot_poloidal(khi_noglm_H_H2_es, 'm', "split DG ES (no GLM)", linestyle=(0, (3, 1, 1, 1, 1, 1)))
plot_poloidal(khi_noglm_H_H2, 'g--', "split DG EC+LLF (no GLM)")

plot_poloidal(khi_glm_H_H2_central, 'r-.', "std DG ")

# Plot x at the end of crashing simus
plt.plot(khi_glm_H_H2_central['time'][-1], khi_glm_H_H2_central['poloidal_field'][-1] / khi_glm_H_H2_central['poloidal_field'][0], 'rx')
plt.plot(khi_noglm_H_H2['time'][-1], khi_noglm_H_H2['poloidal_field'][-1] / khi_noglm_H_H2['poloidal_field'][0], 'gx')
plt.plot(khi_noglm_H_H2_es['time'][-1], khi_noglm_H_H2_es['poloidal_field'][-1] / khi_noglm_H_H2_es['poloidal_field'][0], 'mx')
plt.xlabel("time")
plt.ylabel("$<B_p^2>$")
plt.legend()

plt.show()

#########################
# Plot divB for multi-species case
plt.figure()

plot_divb_l2(khi_glm_H_H2_es, 'b-', "split DG ES")
plot_divb_l2(khi_glm_H_H2, 'k:', "split DG EC+LLF")
plot_divb_l2(khi_noglm_H_H2_es, 'm', "split DG ES (no GLM)", linestyle=(0, (3, 1, 1, 1, 1, 1)))
plot_divb_l2(khi_noglm_H_H2, 'g--', "split DG EC+LLF (no GLM)")

plot_divb_l2(khi_glm_H_H2_central, 'r-.', "std DG ")

# Plot x at the end of crashed
plt.plot(khi_glm_H_H2_central['time'][-1], khi_glm_H_H2_central['l2_divb'][-1], 'rx')
plt.plot(khi_noglm_H_H2['time'][-1], khi_noglm_H_H2['l2_divb'][-1], 'gx')
plt.plot(khi_noglm_H_H2_es['time'][-1], khi_noglm_H_H2_es['l2_divb'][-1], 'mx')
plt.xlabel("time")
plt.ylabel("$\\|\\nabla \\cdot \\vec{B}\\|$")
plt.legend()

plt.show()

