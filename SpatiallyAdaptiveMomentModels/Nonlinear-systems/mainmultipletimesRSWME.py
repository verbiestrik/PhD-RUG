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
from Functions_for_mainRSWME import *
import configparser
from pathlib import Path

def main():
    
    NUM_CONFIGS = 9
    IS_1D_KEY = '1D'
    IS_CLASSICAL_KEY = 'classical'
    eps2 = 2
    CONFIG_DIR = Path("Config-filesMT")
    
    config_files = [str(CONFIG_DIR / f"config{chr(96 + i)}.txt") for i in range(1, NUM_CONFIGS + 1)]

    configs = load_configs(config_files, NUM_CONFIGS)
    pdes = create_pdes(configs, eps2, NUM_CONFIGS)
    spatialDiscretizations = create_spatial_discretizations(configs, NUM_CONFIGS)
    timeintegrations = create_time_integrations(configs, NUM_CONFIGS)
    is_1d_list = extract_is_1d(configs, NUM_CONFIGS, IS_1D_KEY)
    is_classical_list = extract_is_classical(configs, NUM_CONFIGS, IS_CLASSICAL_KEY)
    
    if all(is_1d_list): 
        meshes  = create_meshes(configs, NUM_CONFIGS)
        if  all(is_classical_list):
            simulations = create_simulations(configs, pdes, meshes, spatialDiscretizations, timeintegrations, NUM_CONFIGS, eps2)
            results,timings = run_simulations(configs, simulations, NUM_CONFIGS)
        else:
            "Error no classical methods"
        
        data_sets = [
            (configs[f"C{i}"]["numerical_method_information"].getfloat('t_end'),
            results[f"S{i}"], results[f"S{i+1}"], results[f"S{i+2}"])
            for i in range(1, NUM_CONFIGS+1, 3)
        ]

        nx = configs['C1']['grid_information'].getint('resolutionX')
        results_list = []
     
        for t_end, ref, swe, rswe in data_sets:
           
            errors = calculate_errors(ref, swe, rswe, nx=nx)
            results_list.append({
                "t": t_end,
                "SWE_h": errors[0],
                "SWE_u": errors[1],
                "RSWE_h": errors[2],
                "RSWE_u": errors[3]})
            
        for res in results_list:
            print(f"for t = {res['t']}")
            print(f"Error for SWE in h: {res['SWE_h']}")
            print(f"Error for SWE in u: {res['SWE_u']}")
            print(f"Error for RSWME in h: {res['RSWE_h']}")
            print(f"Error for RSWME in u: {res['RSWE_u']}\n")

        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
       
        for i, (eps, swme, swe, r_swme1) in enumerate(data_sets):
            
            axes[0, i].plot(
                meshes['M1'].cell_center_positions, swme[:, 1],
                label="SWME", linewidth=2.5, color='blue', linestyle=':'
            )

            axes[0, i].plot(
                meshes['M3'].cell_center_positions, r_swme1[:, 1],
                label="RSWME", linewidth=2.5, color='limegreen'
            )
           
            axes[0, i].plot(
                meshes['M2'].cell_center_positions, swe[:, 1],
                label="SWE", linewidth=2.5, color='red', linestyle='--'
            )
            axes[0, i].set_title(f"$t = {eps:.2f}$", fontsize=18, fontweight='bold')
            axes[0, i].set_ylabel("h", fontsize=16, fontweight='bold')
            #axes[0, i].set_ylim(1.03 , 1.17)
            axes[0, i].set_xlabel("x", fontsize=16, fontweight='bold')
            axes[0, i].tick_params(axis='both', which='major', labelsize=14)
            axes[0, i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

           
            axes[1, i].plot(
                meshes['M1'].cell_center_positions, swme[:, 2],
                label="SWME", linewidth=2.5, color='blue', linestyle=':'
            )

            axes[1, i].plot(
                meshes['M2'].cell_center_positions, r_swme1[:, 2],
                label="RSWME", linewidth=2.5, color='limegreen'
            )
            axes[1, i].plot(
                meshes['M3'].cell_center_positions, swe[:, 2],
                label="SWE", linewidth=2.5, color='red', linestyle='--'
            )
            axes[1, i].set_xlabel("x", fontsize=16, fontweight='bold')
            axes[1, i].set_ylabel("u_m", fontsize=16, fontweight='bold')
            #axes[1, i].set_ylim(-0.005,0.17)
            axes[1, i].tick_params(axis='both', which='major', labelsize=14)
            axes[1, i].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

        
        lines, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(lines, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08),
           fontsize=16, framealpha=1, edgecolor='black', shadow=True)
        
        plt.tight_layout()
        plt.savefig("OverviewTimeDiffIII.png", bbox_inches='tight', dpi=300)
        plt.close()

        
        
    else:
        print('2D not implemented yet')

  

if __name__ == '__main__':
    main()