import simulationRSWME
import pdeRSWME
import meshRSWME
import spatialDiscretizationRSWME
import timeIntegrationRSwME
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


NUM_CONFIGS = 7
filenames = [f"Numerical_Results/TestIII/data_set{i}.csv" for i in range(1, NUM_CONFIGS + 1)]
Nlist = [2,2,4,4,6,6] 
data_sets = [read_csv_file(filenames[i])
            for i in range(0, NUM_CONFIGS-1)
        ]
SWE =  read_csv_file(filenames[NUM_CONFIGS-1])

fig, axes = plt.subplots(1, 2, figsize=(15, 8))
       
axes[0].plot(
        SWE[:,0], SWE[:,1],
        label="SWE", linewidth=2.5, color='red'
) 
    
axes[0].plot(
        data_sets[1][:,0], data_sets[1][:,1],
        label="RSWME2", linewidth=2.5, color='limegreen'
    )
axes[0].plot(
        data_sets[3][:,0], data_sets[3][:,1],
        label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
    )

axes[0].plot(
        data_sets[5][:,0], data_sets[5][:,1],
        label="RSWME6", linewidth=2.5, color='darkslategray', linestyle=':'
    )  

axes[0].plot(
        data_sets[0][:,0], data_sets[0][:,1],
        label="SWME2", linewidth=2.5,  color='blue'
    )

axes[0].plot(
        data_sets[2][:,0], data_sets[2][:,1],
        label="SWME4", linewidth=2.5, color='cyan', linestyle='-.'
    ) 

axes[0].plot(
        data_sets[4][:,0], data_sets[4][:,1],
        label="SWME6", linewidth=2.5,  color='magenta', linestyle=':'
    )

    
axes[0].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes[0].yaxis.set_major_locator(MaxNLocator(nbins=6))
axes[0].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
offset_text = axes[0].yaxis.get_offset_text()
offset_text.set_fontsize(14)  
axes[0].set_xticks([-1, -0.5, 0, 0.5, 1])
axes[0].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
axes[0].set_ylabel("$h$", fontsize=17, fontweight='bold')
axes[0].set_xlabel("$x$", fontsize=17, fontweight='bold')
axes[0].tick_params(axis='both', which='major', labelsize=14)
axes[0].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

axes[1].plot(
        SWE[:,0], SWE[:,2],
        label="SWE", linewidth=2.5, color='red'
) 
    
axes[1].plot(
        data_sets[1][:,0], data_sets[1][:,2],
        label="RSWME2", linewidth=2.5, color='limegreen'
    )
axes[1].plot(
        data_sets[3][:,0], data_sets[3][:,2],
        label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
    )

axes[1].plot(
        data_sets[5][:,0], data_sets[5][:,2],
        label="RSWME6", linewidth=2.5, color='darkslategray', linestyle=':'
    )  

axes[1].plot(
        data_sets[0][:,0], data_sets[0][:,2],
        label="SWME2", linewidth=2.5,  color='blue'
    )

axes[1].plot(
        data_sets[2][:,0], data_sets[2][:,2],
        label="SWME4", linewidth=2.5, color='cyan', linestyle='-.'
    ) 

axes[1].plot(
        data_sets[4][:,0], data_sets[4][:,2],
        label="SWME6", linewidth=2.5,  color='magenta', linestyle=':'
    )

    
axes[1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes[1].yaxis.set_major_locator(MaxNLocator(nbins=6))
axes[1].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
offset_text = axes[1].yaxis.get_offset_text()
offset_text.set_fontsize(14)  
axes[1].set_xticks([-1, -0.5, 0, 0.5, 1])
axes[1].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
axes[1].set_ylabel("$u_m$", fontsize=17, fontweight='bold')
axes[1].set_xlabel("$x$", fontsize=17, fontweight='bold')
axes[1].tick_params(axis='both', which='major', labelsize=14)
axes[1].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)


lines, labels = axes[0].get_legend_handles_labels()
rswme_handles = [lines[1], lines[2], lines[3]]
rswme_labels = [labels[1], labels[2], labels[3]]

swme_handles = [lines[4], lines[5], lines[6]]
swme_labels = [labels[4], labels[5], labels[6]]

swe_handles = [lines[0]]
swe_labels = [labels[0]]
#all_handles = rswme_handles + swme_handles + swe_handles
#all_labels = rswme_labels + swme_labels + swe_labels
#fig.legend(all_handles, all_labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15),
#   fontsize=16, framealpha=1, edgecolor='black', shadow=True)

legend1 = fig.legend(
    rswme_handles, rswme_labels,
    loc='upper left',
    bbox_to_anchor=(0.20, 1.15),
    fontsize=16,
    framealpha=1,
    edgecolor='black',
    shadow=True
)

legend2 = fig.legend(
    swme_handles, swme_labels,
    loc='upper center',
    bbox_to_anchor=(0.5, 1.15),
    fontsize=16,
    framealpha=1,
    edgecolor='black',
    shadow=True
)

legend3 = fig.legend(
    swe_handles, swe_labels,
    loc='upper right',
    bbox_to_anchor=(0.80, 1.15),
    fontsize=16,
    framealpha=1,
    edgecolor='black',
    shadow=True
)

# Add all legends to the axes
axes[0].add_artist(legend1)
axes[0].add_artist(legend2)
axes[0].add_artist(legend3)
plt.tight_layout()
plt.savefig("OverviewplotsCaseIIIb.pdf", bbox_inches='tight')
plt.close()
    



fig2, axes2 = plt.subplots(3, 2, figsize=(15, 14))
         
axes2[0,0].plot(
        SWE[:,0], 0.0*SWE[:,0],
        label="SWE", linewidth=2.5, color='red'
) 
    
axes2[0,0].plot(
        data_sets[1][:,0], data_sets[1][:,3],
        label="RSWME2", linewidth=2.5, color='limegreen'
    )
axes2[0,0].plot(
        data_sets[3][:,0], data_sets[3][:,3],
        label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
    )

axes2[0,0].plot(
        data_sets[5][:,0], data_sets[5][:,3],
        label="RSWME6", linewidth=2.5, color='darkslategray' , linestyle=':'
    )

axes2[0,0].plot(
        data_sets[0][:,0], data_sets[0][:,3],
        label="SWME2", linewidth=2.5,  color='blue'
    )

axes2[0,0].plot(
        data_sets[2][:,0], data_sets[2][:,3],
        label="SWME4", linewidth=2.5, color='cyan' , linestyle='-.'
    )
    
axes2[0,0].plot(
        data_sets[0][:,0], data_sets[4][:,3],
        label="SWME6", linewidth=2.5,  color='magenta', linestyle=':'
    )
axes2[0, 0].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes2[0, 0].yaxis.set_major_locator(MaxNLocator(nbins=6))
axes2[0, 0].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
offset_text = axes2[0, 0].yaxis.get_offset_text()
offset_text.set_fontsize(14)  
axes2[0, 0].set_xticks([-1, -0.5, 0, 0.5, 1])
axes2[0, 0].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
axes2[0, 0].set_ylabel(r'$\alpha_1$', fontsize=17, fontweight='bold')
axes2[0, 0].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
axes2[0, 0].tick_params(axis='both', which='major', labelsize=14)
axes2[0, 0].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

axes2[0, 1].plot(
        SWE[:,0], 0.0*SWE[:,0],
        label="SWE", linewidth=2.5, color='red'
) 
    
axes2[0, 1].plot(
        data_sets[1][:,0], data_sets[1][:,4],
        label="RSWME2", linewidth=2.5, color='limegreen'
    )
axes2[0,1].plot(
        data_sets[3][:,0], data_sets[3][:,4],
        label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
    )
axes2[0,1].plot(
        data_sets[5][:,0], data_sets[5][:,4],
        label="RSWME6", linewidth=2.5, color='darkslategray' , linestyle=':'
    )
    

axes2[0,1].plot(
        data_sets[0][:,0], data_sets[0][:,4],
        label="SWME2", linewidth=2.5,  color='blue'
    )

axes2[0,1].plot(
        data_sets[2][:,0], data_sets[2][:,4],
        label="SWME4", linewidth=2.5, color='cyan' , linestyle='-.'
    )

axes2[0,1].plot(
        data_sets[4][:,0], data_sets[4][:,4],
        label="SWME6", linewidth=2.5,  color='magenta', linestyle=':'
    )

axes2[0, 1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes2[0, 1].yaxis.set_major_locator(MaxNLocator(nbins=6))
axes2[0, 1].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
offset_text = axes2[0, 1].yaxis.get_offset_text()
offset_text.set_fontsize(14)  
axes2[0, 1].set_xticks([-1, -0.5, 0, 0.5, 1])
axes2[0, 1].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
axes2[0, 1].set_ylabel(r'$\alpha_2$', fontsize=17, fontweight='bold')
axes2[0, 1].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
axes2[0, 1].tick_params(axis='both', which='major', labelsize=14)
axes2[0, 1].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

axes2[1,0].plot(
        SWE[:,0], 0.0*SWE[:,0],
        label="SWE", linewidth=2.5, color='red'
) 
    
axes2[1,0].plot(
        data_sets[1][:,0], 0.0*data_sets[1][:,0],
        label="RSWME2", linewidth=2.5, color='limegreen', 
    )
axes2[1,0].plot(
        data_sets[3][:,0], data_sets[3][:,5],
        label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
    )

axes2[1,0].plot(
        data_sets[5][:,0], data_sets[5][:,5],
        label="RSWME6", linewidth=2.5, color='darkslategray' , linestyle=':'
    )

axes2[1,0].plot(
        data_sets[0][:,0], 0.0*data_sets[0][:,0],
        label="SWME2", linewidth=2.5,  color='blue'
)

axes2[1,0].plot(
        data_sets[2][:,0], data_sets[2][:,5],
        label="SWME4", linewidth=2.5, color='cyan' , linestyle='-.'
    )

axes2[1,0].plot(
        data_sets[4][:,0], data_sets[4][:,5],
        label="SWME6", linewidth=2.5,  color='magenta', linestyle=':'
    )


axes2[1, 0].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes2[1, 0].yaxis.set_major_locator(MaxNLocator(nbins=6))
axes2[1, 0].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
offset_text = axes2[1, 0].yaxis.get_offset_text()
offset_text.set_fontsize(14)  
axes2[1, 0].set_xticks([-1, -0.5, 0, 0.5, 1])
axes2[1, 0].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
axes2[1, 0].set_ylabel(r'$\alpha_3$', fontsize=17, fontweight='bold')
axes2[1, 0].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
axes2[1, 0].tick_params(axis='both', which='major', labelsize=14)
axes2[1, 0].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

axes2[1,1].plot(
        SWE[:,0], 0.0*SWE[:,0],
        label="SWE", linewidth=2.5, color='red'
    )
axes2[1,1].plot(
        data_sets[1][:,0], 0.0*data_sets[1][:,0],
        label="RSWME2", linewidth=2.5, color='limegreen'
    )
axes2[1,1].plot(
        data_sets[3][:,0], data_sets[3][:,6],
        label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
    )

axes2[1,1].plot(
        data_sets[5][:,0], data_sets[5][:,6],
        label="RSWME6", linewidth=2.5, color='darkslategray' , linestyle=':'
    )
    

axes2[1,1].plot(
        data_sets[0][:,0], 0.0*data_sets[0][:,0],
        label="SWME2", linewidth=2.5,  color='blue'
    )

axes2[1,1].plot(
        data_sets[2][:,0], data_sets[2][:,6],
        label="SWME4", linewidth=2.5, color='cyan' , linestyle='-.'
    )

axes2[1,1].plot(
        data_sets[4][:,0], data_sets[4][:,6],
        label="SWME6", linewidth=2.5,  color='magenta', linestyle=':'
    )

axes2[1, 1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes2[1, 1].yaxis.set_major_locator(MaxNLocator(nbins=6))
axes2[1, 1].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
offset_text = axes2[1, 1].yaxis.get_offset_text()
offset_text.set_fontsize(14)  
axes2[1, 1].set_xticks([-1, -0.5, 0, 0.5, 1])
axes2[1, 1].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
axes2[1, 1].set_ylabel(r'$\alpha_4$', fontsize=17, fontweight='bold')
axes2[1, 1].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
axes2[1, 1].tick_params(axis='both', which='major', labelsize=14)
axes2[1, 1].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

axes2[2,0].plot(
        SWE[:,0], 0.0*SWE[:,0],
        label="SWE", linewidth=2.5, color='red'
) 
    
axes2[2,0].plot(
        data_sets[1][:,0],0.0*data_sets[1][:,0],
        label="RSWME2", linewidth=2.5, color='limegreen'
    )
axes2[2,0].plot(
        data_sets[3][:,0],0.0*data_sets[3][:,0],
        label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
    )

axes2[2,0].plot(
        data_sets[5][:,0], data_sets[5][:,7],
        label="RSWME6", linewidth=2.5, color='darkslategray' , linestyle=':'
    )
    

axes2[2,0].plot(
        data_sets[0][:,0],0.0*data_sets[0][:,0],
        label="SWME2", linewidth=2.5,  color='blue'
    )

axes2[2,0].plot(
        data_sets[2][:,0],0.0*data_sets[2][:,0],
        label="SWME4", linewidth=2.5, color='cyan', linestyle='-.'
    )

axes2[2,0].plot(
        data_sets[4][:,0],data_sets[4][:,7],
        label="SWME6", linewidth=2.5,  color='magenta', linestyle=':'
    )

axes2[2, 0].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes2[2, 0].yaxis.set_major_locator(MaxNLocator(nbins=6))
axes2[2, 0].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
offset_text = axes2[2, 0].yaxis.get_offset_text()
offset_text.set_fontsize(14)  
axes2[2, 0].set_xticks([-1, -0.5, 0, 0.5, 1])
axes2[2, 0].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])
axes2[2, 0].set_ylabel(r'$\alpha_5$', fontsize=17, fontweight='bold')
axes2[2, 0].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
axes2[2, 0].tick_params(axis='both', which='major', labelsize=14)
axes2[2, 0].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

axes2[2,1].plot(
        SWE[:,0], 0.0*SWE[:,0],
        label="SWE", linewidth=2.5, color='red'
    )
axes2[2,1].plot(
        data_sets[1][:,0], 0.0*data_sets[1][:,0],
        label="RSWME2", linewidth=2.5, color='limegreen', linestyle="--"
    )
axes2[2,1].plot(
        data_sets[3][:,0], 0.0*data_sets[3][:,0],
        label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
    )
axes2[2,1].plot(
        data_sets[5][:,0], data_sets[5][:,8],
        label="RSWME6", linewidth=2.5, color='darkslategray' , linestyle=':'
    )

axes2[2,1].plot(
        data_sets[0][:,0], 0.0*data_sets[0][:,0],
        label="SWME2", linewidth=2.5,  color='blue'
    )

axes2[2,1].plot(
        data_sets[2][:,0], 0.0*data_sets[2][:,0],
        label="SWME4", linewidth=2.5, color='cyan' , linestyle='-.'
    )

axes2[2,1].plot(
        data_sets[4][:,0], data_sets[4][:,8],
        label="SWME6", linewidth=2.5,  color='magenta', linestyle=':'
    )


axes2[2, 1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes2[2, 1].yaxis.set_major_locator(MaxNLocator(nbins=6))
axes2[2, 1].ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)  
offset_text = axes2[2, 1].yaxis.get_offset_text()
offset_text.set_fontsize(14)  
axes2[2, 1].set_xticks([-1, -0.5, 0, 0.5, 1])
axes2[2, 1].set_xticklabels(['-1', '-0.5', '0', '0.5', '1'])    
axes2[2, 1].set_ylabel(r'$\alpha_6$', fontsize=17, fontweight='bold')
axes2[2, 1].set_xlabel(r"$x$", fontsize=17, fontweight='bold')
axes2[2, 1].tick_params(axis='both', which='major', labelsize=14)
axes2[2, 1].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)


lines, labels = axes2[0, 0].get_legend_handles_labels()
rswme_handles = [lines[1], lines[2], lines[3]]
rswme_labels = [labels[1], labels[2], labels[3]]

swme_handles = [lines[4], lines[5], lines[6]]
swme_labels = [labels[4], labels[5], labels[6]]

swe_handles = [lines[0]]
swe_labels = [labels[0]]
#all_handles = rswme_handles + swme_handles + swe_handles
#all_labels = rswme_labels + swme_labels + swe_labels
#fig2.legend(all_handles, all_labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08),
#    fontsize=16, framealpha=1, edgecolor='black', shadow=True)

legend1 = fig2.legend(
    rswme_handles, rswme_labels,
    loc='upper left',
    bbox_to_anchor=(0.2, 1.08),
    fontsize=16,
    framealpha=1,
    edgecolor='black',
    shadow=True
)

legend2 = fig2.legend(
    swme_handles, swme_labels,
    loc='upper center',
    bbox_to_anchor=(0.5, 1.08),
    fontsize=16,
    framealpha=1,
    edgecolor='black',
    shadow=True
)

legend3 = fig2.legend(
    swe_handles, swe_labels,
    loc='upper right',
    bbox_to_anchor=(0.80, 1.08),
    fontsize=16,
    framealpha=1,
    edgecolor='black',
    shadow=True
)


axes2[0, 0].add_artist(legend1)
axes2[0, 0].add_artist(legend2)
axes2[0, 0].add_artist(legend3)

plt.tight_layout()
plt.savefig("OverviewThirdPlotMomentsIIIb.pdf",format='pdf', bbox_inches='tight')
plt.close()

