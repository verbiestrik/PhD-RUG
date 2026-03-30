import simulationRSWME
import pdeRSWME
import meshRSWME
import spatialDiscretizationRSWME
import timeIntegrationRSWME
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import configparser
import timeit
from matplotlib import rcParams
from matplotlib import rc
import csv
rc('text', usetex=False)
rcParams['text.usetex'] = False
from Functions_for_mainRSWME import *
import configparser
from pathlib import Path
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import csv

def read_csv_file(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        metadata = [next(reader) for _ in range(5)]
        headers = next(reader)
        next(reader)
        data = []
        for row in reader:
            numeric_row = [float(value) for value in row]
            data.append(numeric_row)
    return np.array(data)


NUM_CONFIGS = 9
filenames = [f"Numerical_Results/TestIb/data_set{i}.csv" for i in range(1, NUM_CONFIGS + 1)]
vareps = [0.01,0.01,0.01,0.1,0.1,0.1,1.0,1.0,1.0] 
data_sets = [
            (vareps[i],
            read_csv_file(filenames[i]), read_csv_file(filenames[i+1]), read_csv_file(filenames[i+2]))
            for i in range(0, NUM_CONFIGS, 3)
        ]

fig, axes = plt.subplots(2, 3, figsize=(15, 8))
       
h_min, h_max, u_min, u_max, a1_min, a1_max = [], [], [], [], [], []
for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
    h_min.append(min(np.min(swme[:,1]),np.min(swe[:,1]),np.min(r_swme[:,1])))
    h_max.append(max(np.max(swme[:,1]),np.max(swe[:,1]),np.max(r_swme[:,1])))
    u_min.append(min(np.min(swme[:,2]),np.min(swe[:,2]),np.min(r_swme[:,2])))
    u_max.append(max(np.max(swme[:,2]),np.max(swe[:,2]),np.max(r_swme[:,2])))
    a1_min.append(min(np.min(swme[:,3]),np.min(r_swme[:,3])))
    a1_max.append(max(np.max(swme[:,3]),np.max(r_swme[:,3])))
hmin, hmax = min(h_min), max(h_max)
umin, umax = min(u_min), max(u_max)
a1min, a1max = min(a1_min), max(a1_max)
margeb = abs((a1max+umax)-(a1min+umin))
margea = abs(a1max-a1min)
margeu = abs((umax)-(umin))
margeh = abs(hmax-hmin)

for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
    axes[0, i].plot(
        swe[:,0], swe[:, 1],
        label="SWE", linewidth=2.5, color='red'
    )
           
    axes[0, i].plot(
        r_swme[:,0], r_swme[:, 1], linestyle='--',
        label="RSWME1", linewidth=2.5, color='limegreen'
        )
    
    axes[0, i].plot(
        swme[:,0], swme[:, 1],
        label="SWME1", linewidth=2.5, color='blue', linestyle=':'
    )
           
    axes[0, i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axes[0, i].yaxis.set_major_locator(MaxNLocator(nbins=6))
    axes[0, i].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
    offset_text = axes[0,i].yaxis.get_offset_text()
    offset_text.set_fontsize(14)  
    axes[0, i].set_xticks([-1, -0.5, 0, 0.5, 1])
    axes[0, i].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
    axes[0, i].set_title(f"$\epsilon = {eps:.2f}$", fontsize=18, fontweight='bold')
    axes[0, i].set_ylabel(r'$h$', fontsize=17, fontweight='bold')
    axes[0, i].set_ylim(hmin-0.05*margeh , hmax+0.05*margeh)
    axes[0, i].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
    axes[0, i].tick_params(axis='both', which='major', labelsize=14)
    axes[0, i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
    
    
    axes[1, i].plot(
        swe[:,0], swe[:, 2],
        label="SWE", linewidth=2.5, color='red'
    )
    axes[1, i].plot(
        r_swme[:,0], r_swme[:, 2], linestyle='--',
        label="RSWME1", linewidth=2.5, color='limegreen'
    )
    
    axes[1, i].plot(
        swme[:,0], swme[:, 2],
        label="SWME1", linewidth=2.5, color='blue', linestyle=':'
    )

    axes[1, i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axes[1, i].yaxis.set_major_locator(MaxNLocator(nbins=6))
    axes[1, i].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
    offset_text = axes[1,i].yaxis.get_offset_text()
    offset_text.set_fontsize(14)  
    axes[1, i].set_xticks([-1, -0.5, 0, 0.5, 1])
    axes[1, i].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
    axes[1, i].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
    axes[1, i].set_ylabel(r"$u_m$", fontsize=17, fontweight='bold')
    axes[1, i].set_ylim(umin-0.05*margeu,umax+0.05*margeu)
    axes[1, i].tick_params(axis='both', which='major', labelsize=14)
    axes[1, i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)   
lines, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(lines, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08),
    fontsize=16, framealpha=1, edgecolor='black', shadow=True)

plt.tight_layout()
plt.savefig("OverviewFirstPlotIb.pdf",format='pdf', bbox_inches='tight')
plt.close()        

fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))

for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
    
    axes2[i].plot(
        swme[:,0], 0.0*swme[:, 3],
        label="SWE", linewidth=2.5, color='red'
    )
    axes2[i].plot(
        r_swme[:,0], r_swme[:, 3],
        label="RSWME1", linewidth=2.5, color='limegreen', linestyle="--"
    )
    axes2[i].plot(
        swme[:,0], swme[:, 3],
        label="SWME1", linewidth=2.5, color='blue', linestyle=':'
    )
    
    axes2[i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axes2[i].yaxis.set_major_locator(MaxNLocator(nbins=6))
    axes2[i].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
    offset_text = axes2[i].yaxis.get_offset_text()
    offset_text.set_fontsize(14)  
    axes2[i].set_xticks([-1, -0.5, 0, 0.5, 1])
    axes2[i].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
    axes2[i].set_title(f"$\epsilon = {eps:.2f}$", fontsize=18, fontweight='bold')
    axes2[i].set_ylabel(r'$\alpha_1$', fontsize=17, fontweight='bold')
    axes2[i].set_ylim(-0.03, 0.01)
    print("sharpgradienttest")
    print(f"max-val for alpha1: {a1max:.5f}")
    print(f"min-val for alpha1: {a1min:.5f}")
    axes2[i].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
    axes2[i].tick_params(axis='both', which='major', labelsize=14)
    axes2[i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

    lines, labels = axes2[0].get_legend_handles_labels()
    fig2.legend(lines, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.10),
        fontsize=16, framealpha=1, edgecolor='black', shadow=True)

plt.tight_layout()
plt.savefig("OverviewFirstPlotMomentIb.pdf", format='pdf', bbox_inches='tight', dpi=500)
plt.close()

fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))

for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
    
    axes3[i].plot(
        swe[:,0], swe[:,2],
        label="SWE", linewidth=2.5, color='red'
    )

    axes3[i].plot(
        r_swme[:,0], r_swme[:,2]+r_swme[:, 3],
        label="RSWME1", linewidth=2.5, color='limegreen', linestyle='--'
    )
    axes3[i].plot(
        swme[:,0], swme[:,2]+swme[:, 3],
        label="SWME1", linewidth=2.5, color='blue', linestyle=':'
    )
    axes3[i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axes3[i].yaxis.set_major_locator(MaxNLocator(nbins=6))
    axes3[i].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
    offset_text = axes3[i].yaxis.get_offset_text()
    offset_text.set_fontsize(14)  
    axes3[i].set_xticks([-1, -0.5, 0, 0.5, 1])
    axes3[i].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
    axes3[i].set_title(f"$\epsilon = {eps:.2f}$", fontsize=18, fontweight='bold')
    axes3[i].set_ylabel(r'$u_b$', fontsize=17, fontweight='bold')
    axes3[i].set_ylim(a1min+umin-0.05*margeb, a1max+umax+0.05*margeb)
    axes3[i].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
    axes3[i].tick_params(axis='both', which='major', labelsize=14)
    axes3[i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
    
     
    lines, labels = axes3[0].get_legend_handles_labels()
    fig3.legend(lines, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.10),
        fontsize=16, framealpha=1, edgecolor='black', shadow=True)

plt.tight_layout()
plt.savefig("OverviewFirstPlotBottomIb.pdf",format='pdf', bbox_inches='tight')
plt.close()

dh_min, dh_max = [], []
dx = 2.0/1000.0
for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
    dh_min.append(min(np.min(Godunovmin(swe[:,1],dx)), np.min(Godunovmin(r_swme[:,1],dx)), np.min(Godunovmin(swme[:,1],dx)))*(1.05))
    dh_max.append(max(np.max(Godunovmin(swe[:,1],dx)), np.max(Godunovmin(r_swme[:,1],dx)), np.max(Godunovmin(swme[:,1],dx)))*(1.05))           
dhmin, dhmax = min(dh_min), max(dh_max)
margedh = abs(hmax-hmin)
dh_lim_min, dh_lim_max = dhmin-0.05*margedh, dhmax+0.05*margedh


fig4, axes4 = plt.subplots(1, 3, figsize=(15, 5))

for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
    
    axes4[i].plot(
        swe[:,0], Godunovmin(swe[:,1],dx),
        label="SWE", linewidth=2.5, color='red'
    )    
    axes4[i].plot(
        r_swme[:,0], Godunovmin(r_swme[:,1],dx),
        label="RSWME1", linewidth=2.5, color='limegreen', linestyle='--'
    )
    axes4[i].plot(
        r_swme[:,0], Godunovmin(swme[:,1],dx),
        label="SWME1", linewidth=2.5, color='blue', linestyle=':'
    )

    axes4[i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axes4[i].yaxis.set_major_locator(MaxNLocator(nbins=6))
    axes4[i].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
    offset_text = axes4[i].yaxis.get_offset_text()
    offset_text.set_fontsize(14)       
    axes4[i].set_xticks([-1, -0.5, 0, 0.5, 1])
    axes4[i].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
    axes4[i].set_title(f"$\epsilon = {eps:.2f}$", fontsize=18, fontweight='bold')
    axes4[i].set_ylabel(r'$D(-h^4)$', fontsize=17, fontweight='bold')
    axes4[i].set_ylim(dh_lim_min, dh_lim_max)
    axes4[i].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
    axes4[i].tick_params(axis='both', which='major', labelsize=14)
    axes4[i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
    


    lines, labels = axes4[0].get_legend_handles_labels()
    fig4.legend(lines, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.10),
        fontsize=16, framealpha=1, edgecolor='black', shadow=True)

plt.tight_layout()
plt.savefig("OverviewFirstPlotDhIb.pdf",format='pdf', bbox_inches='tight')
plt.close()


