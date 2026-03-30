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
from pathlib import Path

def main():

    NUM_CONFIGS = 9
    IS_1D_KEY = '1D'
    IS_CLASSICAL_KEY = 'classical'
    CONFIG_DIR = Path("Config-filesII")

    config_files = [str(CONFIG_DIR / f"config{chr(96 + i)}.txt") for i in range(1, NUM_CONFIGS + 1)]

    configs = load_configs(config_files, NUM_CONFIGS)
    pdes = create_pdes(configs, NUM_CONFIGS)
    spatialDiscretizations = create_spatial_discretizations(configs, NUM_CONFIGS)
    timeintegrations = create_time_integrations(configs, NUM_CONFIGS)
    is_1d_list = extract_is_1d(configs, NUM_CONFIGS, IS_1D_KEY)
    is_classical_list = extract_is_classical(configs, NUM_CONFIGS, IS_CLASSICAL_KEY)
    
    if all(is_1d_list): 
        meshes  = create_meshes(configs, NUM_CONFIGS)
        if  all(is_classical_list):
            simulations = create_simulations(configs, pdes, meshes, spatialDiscretizations, timeintegrations, NUM_CONFIGS)
            results, timings = run_simulations(configs, simulations, NUM_CONFIGS)
        else:
            "Error no classical methods"
        
        data_sets = [
            (1 / configs[f"C{i}"]['pde_information'].getfloat('slipLength'),
            results[f"S{i}"], results[f"S{i+1}"], results[f"S{i+2}"])
            for i in range(1, NUM_CONFIGS+1, 3)
        ]


        for i in range(1,NUM_CONFIGS+1):
            headers = ['x', 'h', 'um', 'a1', 'a2'] 
            filename = f"Numerical_Results/TestII/data_set{i}.csv"
            eps_val =1/(float(configs[f"C{i}"]['pde_information']['slipLength']))
            order = configs[f"C{i}"]["numerical_method_information"].getfloat('order')
            time = timings[f"S{i}"]
            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([f"N={order}"])
                writer.writerow([f"cpu-time: $t={time}$"])
                writer.writerow([f"$\\varepsilon = {eps_val}$"])
                writer.writerow([configs[f"C{i}"]['pde_information']['pde_type']]) 
                writer.writerow([])
                writer.writerow(headers)
                writer.writerows(results[f"S{i}"])
       
        
        nx = configs['C1']['grid_information'].getint('resolutionX')
        for eps, ref, swe, rswe in data_sets:
            filename = f"Numerical_Results/TestII/Errors_eps_{eps}.csv"
            swe_h = np.linalg.norm(ref[:, 1] - swe[:, 1], 1)/np.linalg.norm(ref[:, 1] , 1)
            swe_u = np.linalg.norm(ref[:, 2] - swe[:, 2], 1)/np.linalg.norm(ref[:, 2] , 1)
            rswe_h = np.linalg.norm(ref[:, 1] - rswe[:, 1], 1)/np.linalg.norm(ref[:, 1] , 1)
            rswe_u = np.linalg.norm(ref[:, 2] - rswe[:, 2], 1)/np.linalg.norm(ref[:, 2] , 1)

            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([f"for eps = {eps}"])
                writer.writerow([f"Error for SWE in h: {swe_h}"])
                writer.writerow([f"Error for SWE in u: {swe_u}"])
                writer.writerow([f"Error for RSWE in h: {rswe_h}"])
                writer.writerow([f"Error for RSWE in u: {rswe_u}\n"])


        h_min, h_max, u_min, u_max, a1_min, a1_max, a2_min,a2_max = [], [], [], [], [], [],[],[]
        for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
            h_min.append(min(np.min(swme[:,1]),np.min(swe[:,1]),np.min(r_swme[:,1])))
            h_max.append(max(np.max(swme[:,1]),np.max(swe[:,1]),np.max(r_swme[:,1])))
            u_min.append(min(np.min(swme[:,2]),np.min(swe[:,2]),np.min(r_swme[:,2])))
            u_max.append(max(np.max(swme[:,2]),np.max(swe[:,2]),np.max(r_swme[:,2])))
            a1_min.append(min(np.min(swme[:,3]),np.min(r_swme[:,3])))
            a1_max.append(max(np.max(swme[:,3]),np.max(r_swme[:,3])))
            a2_min.append(min(np.min(swme[:,4]),np.min(r_swme[:,4])))
            a2_max.append(max(np.max(swme[:,4]),np.max(r_swme[:,4])))
        hmin, hmax = min(h_min), max(h_max)
        umin, umax = min(u_min), max(u_max)
        a1min, a1max =min(a1_min), max(a1_max)
        a2min, a2max =min(a2_min), max(a2_max)
        margeb = abs((a1max+a2max+umax)-(a1min+a2min+umin))
        margea1 = abs(a1max-a1min)
        margea2 = abs(a2max-a2min)
        margeu = abs((umax)-(umin))
        margeh = abs(hmax-hmin)
      
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
       
        for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
            axes[0, i].plot(
                meshes['M2'].cell_center_positions, swe[:, 1],
                label="SWE", linewidth=2.5, color='red'
            )

            axes[0, i].plot(
                meshes['M3'].cell_center_positions, r_swme[:, 1], linestyle='--',
                label="RSWME", linewidth=2.5, color='limegreen'
            )
           
            axes[0, i].plot(
                meshes['M1'].cell_center_positions, swme[:, 1],
                label="SWME", linewidth=2.5, color='blue', linestyle=':'
            )
 
            axes[0, i].set_title(f"$\epsilon = {eps:.2f}$", fontsize=18, fontweight='bold')
            axes[0, i].set_ylabel(r'$h$', fontsize=16, fontweight='bold')
            axes[0, i].set_ylim(hmin-0.05*margeh , hmax+0.05*margeh)
            axes[0, i].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
            axes[0, i].tick_params(axis='both', which='major', labelsize=14)
            axes[0, i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

            axes[1, i].plot(
                meshes['M2'].cell_center_positions, swe[:, 2],
                label="SWE", linewidth=2.5, color='red'
            )

            axes[1, i].plot(
                meshes['M3'].cell_center_positions, r_swme[:, 2], linestyle='--',
                label="RSWME", linewidth=2.5, color='limegreen'
            )
            
            axes[1, i].plot(
                meshes['M1'].cell_center_positions, swme[:, 2],
                label="SWME", linewidth=2.5, color='blue', linestyle=':'
            )

            axes[1, i].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
            axes[1, i].set_ylabel(r"$u_m$", fontsize=16, fontweight='bold')
            axes[1, i].set_ylim(umin-0.05*margeu,umax+0.05*margeu)
            axes[1, i].tick_params(axis='both', which='major', labelsize=14)
            axes[1, i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

        
        lines, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(lines, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08),
           fontsize=16, framealpha=1, edgecolor='black', shadow=True)
        
        plt.tight_layout()
        plt.savefig("OverviewSecondPlot.pdf",format='pdf', bbox_inches='tight')
        plt.close()
 
        fig2, axes2 = plt.subplots(2, 3, figsize=(15, 8))
       
        for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
            axes2[0, i].plot(
                meshes['M2'].cell_center_positions, 0.0*swme[:, 3],
                label="SWE", linewidth=2.5, color='red'
            ) 
            
            axes2[0, i].plot(
                meshes['M3'].cell_center_positions, r_swme[:, 3],
                label="RSWME", linewidth=2.5, color='limegreen', linestyle="--"
            )

            axes2[0, i].plot(
                meshes['M1'].cell_center_positions, swme[:, 3],
                label="SWME", linewidth=2.5, color='blue', linestyle=':'
            )

           
            axes2[0, i].set_title(f"$\epsilon = {eps:.2f}$", fontsize=18, fontweight='bold')
            axes2[0, i].set_ylabel(r'$\alpha_1$', fontsize=16, fontweight='bold')
            axes2[0, i].set_ylim(a1min-0.05*margea1, a1max+0.05*margea1)
            axes2[0, i].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
            axes2[0, i].tick_params(axis='both', which='major', labelsize=14)
            axes2[0, i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

            axes2[1, i].plot(
                meshes['M3'].cell_center_positions, r_swme[:, 4]*0,
                label="SWE", linewidth=2.5, color='red'
            )

            axes2[1, i].plot(
                meshes['M3'].cell_center_positions, r_swme[:, 4],
                label="RSWME", linewidth=2.5, color='limegreen', linestyle='--'
            )
           

            axes2[1, i].plot(
                meshes['M1'].cell_center_positions, swme[:, 4],
                label="SWME", linewidth=2.5, color='blue', linestyle=':'
            )


            axes2[1, i].set_title(f"$\epsilon = {eps:.2f}$", fontsize=18, fontweight='bold')
            axes2[1, i].set_ylabel(r'$\alpha_2$', fontsize=16, fontweight='bold')
            axes2[1, i].set_ylim(a2min-0.05*margea2, a2max+0.05*margea2)
            axes2[1, i].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
            axes2[1, i].tick_params(axis='both', which='major', labelsize=14)
            axes2[1, i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

        lines, labels = axes2[0, 0].get_legend_handles_labels()
        fig2.legend(lines, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08),
            fontsize=16, framealpha=1, edgecolor='black', shadow=True)
        
        plt.tight_layout()
        plt.savefig("OverviewSecondPlotMoments.pdf",format='pdf', bbox_inches='tight')
        #plt.savefig("OverviewSecondPlotMoments.png",format='png', bbox_inches='tight')
        plt.close()

        fig3, axes3 = plt.subplots(1, 3, figsize=(15, 8))
       
        for i, (eps, swme, swe, r_swme) in enumerate(data_sets):
            
            axes3[i].plot(
                meshes['M2'].cell_center_positions, swe[:,2],
                label="SWE", linewidth=2.5, color='red'
            )

            axes3[i].plot(
                meshes['M3'].cell_center_positions, r_swme[:,2]+r_swme[:, 3]+r_swme[:,4],
                label="RSWME", linewidth=2.5, color='limegreen', linestyle='--'
            )

            axes3[i].plot(
                meshes['M1'].cell_center_positions, swme[:,2]+swme[:, 3]+swme[:, 4],
                label="SWME", linewidth=2.5, color='blue', linestyle=':'
            )
            
           
            axes3[i].set_title(f"$\epsilon = {eps:.2f}$", fontsize=18, fontweight='bold')
            axes3[i].set_ylabel(r'$u_b$', fontsize=16, fontweight='bold')
            axes3[i].set_ylim(a1min+umin+a2min-0.05*margeb, a1max+umax+a2max+0.05*margeb)
            axes3[i].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
            axes3[i].tick_params(axis='both', which='major', labelsize=14)
            axes3[i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

            lines, labels = axes3[0].get_legend_handles_labels()
            fig3.legend(lines, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08),
                fontsize=16, framealpha=1, edgecolor='black', shadow=True)
        
        plt.tight_layout()
        plt.savefig("OverviewSecondPlotBottom.pdf",format='pdf', bbox_inches='tight')
        #plt.savefig("OverviewSecondPlotBottom.png",format='png', bbox_inches='tight')
        plt.close()

        
        
    else:
         print('2D not implemented yet')

  

if __name__ == '__main__':
    main()