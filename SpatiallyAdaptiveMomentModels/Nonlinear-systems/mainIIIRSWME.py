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
    
    NUM_CONFIGS = 7
    IS_1D_KEY = '1D'
    IS_CLASSICAL_KEY = 'classical'
    CONFIG_DIR = Path("Config-filesIII")

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
            (configs[f"C{i}"]["numerical_method_information"].getint('order'),
            results[f"S{i}"], results[f"S{i+1}"])
            for i in range(1, NUM_CONFIGS, 2)
        ]

        SWE = results[f"S{NUM_CONFIGS}"]

        for i in range(1,NUM_CONFIGS+1):
            headers = ['x', 'h', 'um', 'a1', 'a2']
            filename = f"Numerical_Results/TestIII/data_set{i}.csv"
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
        for Ni, ref, rswe in data_sets:
            filename = f"Numerical_Results/TestIII/Errors_eps_{Ni}.csv"
            rswe_h = np.linalg.norm(ref[:, 1] - rswe[:, 1], 1)/np.linalg.norm(ref[:, 1], 1)
            rswe_u = np.linalg.norm(ref[:, 2] - rswe[:, 2], 1)/np.linalg.norm(ref[:, 2] , 1)
            swe_h = np.linalg.norm(ref[:,1] - SWE[:,1],1)/np.linalg.norm(ref[:, 1] , 1)
            swe_u = np.linalg.norm(ref[:,2] - SWE[:,2],1)/np.linalg.norm(ref[:, 2] , 1)

            rswe_a1 = np.linalg.norm(ref[:, 3] - rswe[:, 3], 1)/np.linalg.norm(ref[:, 3] , 1)
            rswe_a2 = np.linalg.norm(ref[:, 4] - rswe[:, 4], 1)/np.linalg.norm(ref[:, 4] , 1)
            swe_a1 = 1
            swe_a2 =  1
            rswe_a3 = 0
            swe_a3 = 0
            rswe_a4 = 0
            swe_a4 = 0
            rswe_a5 = 0
            swe_a5 = 0
            rswe_a6 = 0
            swe_a6 = 0
            
            if Ni>2:
                rswe_a3 = np.linalg.norm(ref[:, 5] - rswe[:, 5], 1)/np.linalg.norm(ref[:, 5] , 1)
                swe_a3 = 1
            if Ni>3:            
                rswe_a4 = np.linalg.norm(ref[:, 6] - rswe[:, 6], 1)/np.linalg.norm(ref[:, 6] , 1)
                swe_a4 = 1
            if Ni>4:
                rswe_a5 = np.linalg.norm(ref[:, 7] - rswe[:, 7], 1)/np.linalg.norm(ref[:,7] , 1)
                swe_a5 = 1
            if Ni>5:            
                rswe_a6 = np.linalg.norm(ref[:, 8] - rswe[:, 8], 1)/np.linalg.norm(ref[:, 6] , 1)
                swe_a6 = 1

            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([f"for N = {Ni}"])
                writer.writerow([f"Error for SWE in h: {swe_h}"])
                writer.writerow([f"Error for SWE in u: {swe_u}"])
                writer.writerow([f"Error for RSWE in h: {rswe_h}"])
                writer.writerow([f"Error for RSWE in u: {rswe_u}"])
                writer.writerow([f"Error for SWE in a1: {swe_a1}"])
                writer.writerow([f"Error for SWE in a2: {swe_a2}"])
                writer.writerow([f"Error for RSWE in a1: {rswe_a1}"])
                writer.writerow([f"Error for RSWE in a2: {rswe_a2}"])
                writer.writerow([f"Error for SWE in a3: {swe_a3}"])
                writer.writerow([f"Error for SWE in a4: {swe_a4}"])
                writer.writerow([f"Error for RSWE in a3: {rswe_a3}"])
                writer.writerow([f"Error for RSWE in a4: {rswe_a4}"])
                writer.writerow([f"Error for SWE in a5: {swe_a5}"])
                writer.writerow([f"Error for SWE in a6: {swe_a6}"])
                writer.writerow([f"Error for RSWE in a5: {rswe_a5}"])
                writer.writerow([f"Error for RSWE in a6: {rswe_a6}\n"])
          
          
       
        h_min, h_max, u_min, u_max, a1_min, a1_max, a2_min, a2_max = [], [], [], [], [], [],[],[]

        for i, (Ni, swme, r_swme) in enumerate(data_sets):
            h_min.append(min(np.min(swme[:,1]),np.min(SWE[:,1]),np.min(r_swme[:,1])))
            h_max.append(max(np.max(swme[:,1]),np.max(SWE[:,1]),np.max(r_swme[:,1])))
            u_min.append(min(np.min(swme[:,2]),np.min(SWE[:,2]),np.min(r_swme[:,2])))
            u_max.append(max(np.max(swme[:,2]),np.max(SWE[:,2]),np.max(r_swme[:,2])))
            a1_min.append(min(np.min(swme[:,3]),np.min(r_swme[:,3])))
            a1_max.append(max(np.max(swme[:,3]),np.max(r_swme[:,3])))
            a2_min.append(min(np.min(swme[:,4]),np.min(r_swme[:,4])))
            a2_max.append(max(np.max(swme[:,4]),np.max(r_swme[:,4])))
        hmin, hmax = min(h_min), max(h_max)
        umin, umax = min(u_min), max(u_max)
        a1min, a1max =min(a1_min), max(a1_max)
        a2min, a2max =min(a2_min), max(a2_max)
        margea1 = abs(a1max-a1min)
        margea2 = abs(a2max-a2min)
        margeu = abs((umax)-(umin))
        margeh = abs(hmax-hmin)
      
        # fig, axes = plt.subplots(1, 2, figsize=(15, 8))
       
        # axes[0].plot(
        #         meshes['M2'].cell_center_positions, SWE[:,1],
        #         label="SWE", linewidth=2.5, color='red'
        # ) 
            
        # axes[0].plot(
        #         meshes['M3'].cell_center_positions, results["S2"][:,1],
        #         label="RSWME, N = 2", linewidth=2.5, color='limegreen'
        #     )
        # axes[0].plot(
        #         meshes['M3'].cell_center_positions, results["S4"][:,1],
        #         label="RSWME, N = 4", linewidth=2.5, color='orange', linestyle='-.'
        #     )

        # axes[0].plot(
        #         meshes['M1'].cell_center_positions, results["S6"][:,1],
        #         label="RSWME, N = 6", linewidth=2.5, color='darkslategray', linestyle=':'
        #     )  

        # axes[0].plot(
        #         meshes['M1'].cell_center_positions, results["S1"][:,1],
        #         label="SWME, N = 2", linewidth=2.5,  color='blue'
        #     )

        # axes[0].plot(
        #         meshes['M1'].cell_center_positions, results["S3"][:,1],
        #         label="SWME, N = 4", linewidth=2.5, color='cyan', linestyle='-.'
        #     ) 

        # axes[0].plot(
        #         meshes['M1'].cell_center_positions, results["S5"][:,1],
        #         label="SWME, N = 6", linewidth=2.5,  color='magenta', linestyle=':'
        #     )

          

        # axes[0].set_ylabel("$h$", fontsize=16, fontweight='bold')
        # axes[0].set_xlabel("$x$", fontsize=16, fontweight='bold')
        # axes[0].tick_params(axis='both', which='major', labelsize=14)
        # axes[0].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

        # axes[1].plot(
        #         meshes['M2'].cell_center_positions, SWE[:,2],
        #         label="SWE", linewidth=2.5, color='red'
        # ) 
            
        # axes[1].plot(
        #         meshes['M3'].cell_center_positions, results["S2"][:,2],
        #         label="RSWME, N = 2", linewidth=2.5, color='limegreen'
        #     )
        # axes[1].plot(
        #         meshes['M3'].cell_center_positions, results["S4"][:,2],
        #         label="RSWME, N = 4", linewidth=2.5, color='orange', linestyle='-.'
        #     )
        # axes[1].plot(
        #         meshes['M1'].cell_center_positions, results["S6"][:,2],
        #         label="RSWME, N = 6", linewidth=2.5, color='darkslategray' , linestyle=':'
        #     )   

        # axes[1].plot(
        #         meshes['M1'].cell_center_positions, results["S1"][:,2],
        #         label="SWME, N = 2", linewidth=2.5,  color='blue',
        #     )

        # axes[1].plot(
        #         meshes['M1'].cell_center_positions, results["S3"][:,2],
        #         label="SWME, N = 4", linewidth=2.5, color='cyan' ,linestyle='-.'
        #     )    

        # axes[1].plot(
        #         meshes['M1'].cell_center_positions, results["S5"][:,2],
        #         label="SWME, N = 6", linewidth=2.5,  color='magenta', linestyle=':'
        #     )

         

        # axes[1].set_ylabel("$u_m$", fontsize=16, fontweight='bold')
        # axes[1].set_xlabel("$x$", fontsize=16, fontweight='bold')
        # axes[1].tick_params(axis='both', which='major', labelsize=14)
        # axes[1].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

        
        # lines, labels = axes[0].get_legend_handles_labels()
        # rswme_handles = [lines[1], lines[2], lines[3]]
        # rswme_labels = [labels[1], labels[2], labels[3]]

        # swme_handles = [lines[4], lines[5], lines[6]]
        # swme_labels = [labels[4], labels[5], labels[6]]

        # swe_handles = [lines[0]]
        # swe_labels = [labels[0]]
        # #all_handles = rswme_handles + swme_handles + swe_handles
        # #all_labels = rswme_labels + swme_labels + swe_labels
        # #fig.legend(all_handles, all_labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15),
        # #   fontsize=16, framealpha=1, edgecolor='black', shadow=True)
        
        # legend1 = fig.legend(
        #     rswme_handles, rswme_labels,
        #     loc='upper left',
        #     bbox_to_anchor=(0.20, 1.15),
        #     fontsize=16,
        #     framealpha=1,
        #     edgecolor='black',
        #     shadow=True
        # )

        # legend2 = fig.legend(
        #     swme_handles, swme_labels,
        #     loc='upper center',
        #     bbox_to_anchor=(0.5, 1.15),
        #     fontsize=16,
        #     framealpha=1,
        #     edgecolor='black',
        #     shadow=True
        # )

        # legend3 = fig.legend(
        #     swe_handles, swe_labels,
        #     loc='upper right',
        #     bbox_to_anchor=(0.80, 1.15),
        #     fontsize=16,
        #     framealpha=1,
        #     edgecolor='black',
        #     shadow=True
        # )

        # # Add all legends to the axes
        # axes[0].add_artist(legend1)
        # axes[0].add_artist(legend2)
        # axes[0].add_artist(legend3)
        # plt.tight_layout()
        # plt.savefig("OverviewplotsCaseIIIb.pdf", bbox_inches='tight')
        # plt.close()
    
        fig2, axes2 = plt.subplots(2, 3, figsize=(15, 8))
        z = np.linspace(0,1,100)
        P = int(meshes['M3'].resolution/2)
        
        velocity_profile_swme2 = pdes['pde1'].compute_vertical_velocity_profile(2,
                                                                      results["S1"],
                                                                      z)[P]

        velocity_profile_swme4 = pdes['pde3'].compute_vertical_velocity_profile(4,
                                                                      results["S3"],
                                                                      z)[P]
        velocity_profile_swme6 = pdes['pde5'].compute_vertical_velocity_profile(6,
                                                                      results["S5"],
                                                                      z)[P]
                                                                                                                                 
        velocity_profile_rswme2 = pdes['pde2'].compute_vertical_velocity_profile(2,
                                                                      results["S2"],
                                                                      z)[P]
        velocity_profile_rswme4 = pdes['pde4'].compute_vertical_velocity_profile(4,
                                                                      results["S4"],
                                                                      z)[P]
                                                                    
        velocity_profile_rswme6 = pdes['pde6'].compute_vertical_velocity_profile(4,
                                                                      results["S6"],
                                                                      z)[P]                                                            
        
        velocity_profile_swe = pdes['pde7'].compute_vertical_velocity_profile(0,
                                                                      results["S7"],
                                                                      z)[P]
        
        #fig1, axes1 = plt.subplots(1, 1, figsize=(15, 8))
        fig1 = plt.figure(figsize=(12, 6))
        axes1 = fig1.add_subplot(111)
        
        axes1.plot(
                velocity_profile_swe,z,
                label="SWE", linewidth=2.5, color='red'
            )

        axes1.plot(
                velocity_profile_rswme2,z,
                label="RSWME2", linewidth=2.5, color='limegreen'
            )
        axes1.plot(
                velocity_profile_rswme4,z,
                label="RSWME4", linewidth=2.5, color='orange', linestyle='-.'
            )

        axes1.plot(
                 velocity_profile_rswme6,z,
                label="RSWME6", linewidth=2.5, color='darkslategray' , linestyle=':'
            )

        axes1.plot(
                velocity_profile_swme2,z,
                label="SWME2", linewidth=2.5, color='blue'
            )
        axes1.plot(
                 velocity_profile_swme4,z,
                label="SWME4", linewidth=2.5, color='cyan' , linestyle='-.'
            )
            
        axes1.plot(
                velocity_profile_swme6,z,
                label="SWME6", linewidth=2.5, color='magenta', linestyle=':'
            )
        
        #axes1.set_title("Velocity profiles", fontsize=18, fontweight='bold')
        axes1.set_ylabel("$\zeta$", fontsize=18, fontweight='bold', labelpad=10)
        axes1.set_xlabel("$u(x_0,t,\zeta)$", fontsize=18, fontweight='bold', labelpad=10)
        axes1.set_xlim(0.10,0.18)
        axes1.tick_params(axis='both', which='major', labelsize=14)
        axes1.grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

        
        lines, labels = axes1.get_legend_handles_labels()
        labels = ['SWE', 'RSWME2', 'RSWME4', 'RSWME6', 'SWME2', 'SWME4', 'SWME6']

        rswme_handles = [lines[1], lines[2], lines[3]]
        rswme_labels = [labels[1], labels[2], labels[3]]

        swme_handles = [lines[4], lines[5], lines[6]]
        swme_labels = [labels[4], labels[5], labels[6]]

        swe_handles = [lines[0]]
        swe_labels = [labels[0]]
        #all_handles = rswme_handles + swme_handles + swe_handles
        #all_labels = rswme_labels + swme_labels + swe_labels
        #fig1.legend(all_handles, all_labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.2), 
        #         fontsize=16, framealpha=1, edgecolor='black', shadow=True)


        legend1 = fig1.legend(
            rswme_handles, rswme_labels,
            loc='upper left',
            bbox_to_anchor=(0.2, 1.2),
            fontsize=16,
            framealpha=1,
            edgecolor='black',
            shadow=True
        )

        legend2 = fig1.legend(
            swme_handles, swme_labels,
            loc='upper center',
            bbox_to_anchor=(0.5, 1.2),
            fontsize=16,
            framealpha=1,
            edgecolor='black',
            shadow=True
        )

        legend3 = fig1.legend(
            swe_handles, swe_labels,
            loc='upper right',
            bbox_to_anchor=(0.80, 1.2),
            fontsize=16,
            framealpha=1,
            edgecolor='black',
            shadow=True
        )

        # Add all legends to the axes
        axes1.add_artist(legend1)
        axes1.add_artist(legend2)
        axes1.add_artist(legend3)
        plt.tight_layout()
        plt.savefig("VelocityProfileIIIb.pdf", bbox_inches='tight', dpi=300)
        plt.close()


        # fig2, axes2 = plt.subplots(3, 2, figsize=(15, 14))
       
      
        # axes2[0,0].plot(
        #         meshes['M2'].cell_center_positions, 0.0*results["S2"][:,3],
        #         label="SWE", linewidth=2.5, color='red'
        # ) 
            
        # axes2[0,0].plot(
        #         meshes['M3'].cell_center_positions, results["S2"][:,3],
        #         label="RSWME, N = 2", linewidth=2.5, color='limegreen'
        #     )
        # axes2[0,0].plot(
        #         meshes['M3'].cell_center_positions, results["S4"][:,3],
        #         label="RSWME, N = 4", linewidth=2.5, color='orange', linestyle='-.'
        #     )

        
        # axes2[0,0].plot(
        #         meshes['M1'].cell_center_positions, results["S6"][:,3],
        #         label="RSWME, N = 6", linewidth=2.5, color='darkslategray' , linestyle=':'
        #     )

        # axes2[0,0].plot(
        #         meshes['M1'].cell_center_positions, results["S1"][:,3],
        #         label="SWME, N = 2", linewidth=2.5,  color='blue'
        #     )

        # axes2[0,0].plot(
        #         meshes['M1'].cell_center_positions, results["S3"][:,3],
        #         label="SWME, N = 4", linewidth=2.5, color='cyan' , linestyle='-.'
        #     )
           
        # axes2[0,0].plot(
        #         meshes['M1'].cell_center_positions, results["S5"][:,3],
        #         label="SWME, N = 6", linewidth=2.5,  color='magenta', linestyle=':'
        #     )

        # axes2[0, 0].set_ylabel(r'$\alpha_1$', fontsize=16, fontweight='bold')
        # axes2[0, 0].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
        # axes2[0, 0].tick_params(axis='both', which='major', labelsize=14)
        # axes2[0, 0].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

        # axes2[0,1].plot(
        #         meshes['M2'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="SWE", linewidth=2.5, color='red'
        # ) 
            
        # axes2[0,1].plot(
        #         meshes['M3'].cell_center_positions, results["S2"][:,4],
        #         label="RSWME, N = 2", linewidth=2.5, color='limegreen'
        #     )
        # axes2[0,1].plot(
        #         meshes['M3'].cell_center_positions, results["S4"][:,4],
        #         label="RSWME, N = 4", linewidth=2.5, color='orange', linestyle='-.'
        #     )
        # axes2[0,1].plot(
        #         meshes['M1'].cell_center_positions, results["S6"][:,4],
        #         label="RSWME, N = 6", linewidth=2.5, color='darkslategray' , linestyle=':'
        #     )
           

        # axes2[0,1].plot(
        #         meshes['M1'].cell_center_positions, results["S1"][:,4],
        #         label="SWME, N = 2", linewidth=2.5,  color='blue'
        #     )

        # axes2[0,1].plot(
        #         meshes['M1'].cell_center_positions, results["S3"][:,4],
        #         label="SWME, N = 4", linewidth=2.5, color='cyan' , linestyle='-.'
        #     )

        # axes2[0,1].plot(
        #         meshes['M1'].cell_center_positions, results["S5"][:,4],
        #         label="SWME, N = 6", linewidth=2.5,  color='magenta', linestyle=':'
        #     )

       
        # axes2[0, 1].set_ylabel(r'$\alpha_2$', fontsize=16, fontweight='bold')
        # axes2[0, 1].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
        # axes2[0, 1].tick_params(axis='both', which='major', labelsize=14)
        # axes2[0, 1].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
        
        # axes2[1,0].plot(
        #         meshes['M2'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="SWE", linewidth=2.5, color='red'
        # ) 
            
        # axes2[1,0].plot(
        #         meshes['M3'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="RSWME, N = 2", linewidth=2.5, color='limegreen', 
        #     )
        # axes2[1,0].plot(
        #         meshes['M3'].cell_center_positions, results["S4"][:,5],
        #         label="RSWME, N = 4", linewidth=2.5, color='orange', linestyle='-.'
        #     )

        # axes2[1,0].plot(
        #         meshes['M1'].cell_center_positions, results["S6"][:,5],
        #         label="RSWME, N = 6", linewidth=2.5, color='darkslategray' , linestyle=':'
        #     )

        # axes2[1,0].plot(
        #         meshes['M1'].cell_center_positions, 0.0*results["S1"][:,4],
        #         label="SWME, N = 2", linewidth=2.5,  color='blue'
        # )

        # axes2[1,0].plot(
        #         meshes['M1'].cell_center_positions, results["S3"][:,5],
        #         label="SWME, N = 4", linewidth=2.5, color='cyan' , linestyle='-.'
        #     )
        
        # axes2[1,0].plot(
        #         meshes['M1'].cell_center_positions, results["S5"][:,5],
        #         label="SWME, N = 6", linewidth=2.5,  color='magenta', linestyle=':'
        #     )

        
           
        # axes2[1, 0].set_ylabel(r'$\alpha_3$', fontsize=16, fontweight='bold')
        # axes2[1, 0].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
        # axes2[1, 0].tick_params(axis='both', which='major', labelsize=14)
        # axes2[1, 0].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
        
        # axes2[1,1].plot(
        #         meshes['M3'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="SWE", linewidth=2.5, color='red'
        #     )
        # axes2[1,1].plot(
        #         meshes['M3'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="RSWME, N = 2", linewidth=2.5, color='limegreen'
        #     )
        # axes2[1,1].plot(
        #         meshes['M3'].cell_center_positions, results["S4"][:,6],
        #         label="RSWME, N = 4", linewidth=2.5, color='orange', linestyle='-.'
        #     )

        # axes2[1,1].plot(
        #         meshes['M1'].cell_center_positions, results["S6"][:,6],
        #         label="RSWME, N = 6", linewidth=2.5, color='darkslategray' , linestyle=':'
        #     )
           

        # axes2[1,1].plot(
        #         meshes['M1'].cell_center_positions, 0.0*results["S1"][:,4],
        #         label="SWME, N = 2", linewidth=2.5,  color='blue'
        #     )

        # axes2[1,1].plot(
        #         meshes['M1'].cell_center_positions, results["S3"][:,6],
        #         label="SWME, N = 4", linewidth=2.5, color='cyan' , linestyle='-.'
        #     )

        # axes2[1,1].plot(
        #         meshes['M1'].cell_center_positions, results["S5"][:,6],
        #         label="SWME, N = 6", linewidth=2.5,  color='magenta', linestyle=':'
        #     )

        
        # axes2[1, 1].set_ylabel(r'$\alpha_4$', fontsize=16, fontweight='bold')
        # axes2[1, 1].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
        # axes2[1, 1].tick_params(axis='both', which='major', labelsize=14)
        # axes2[1, 1].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)

        # axes2[2,0].plot(
        #         meshes['M2'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="SWE", linewidth=2.5, color='red'
        # ) 
            
        # axes2[2,0].plot(
        #         meshes['M3'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="RSWME, N = 2", linewidth=2.5, color='limegreen'
        #     )
        # axes2[2,0].plot(
        #         meshes['M3'].cell_center_positions, 0.0*results["S4"][:,5],
        #         label="RSWME, N = 4", linewidth=2.5, color='orange', linestyle='-.'
        #     )

        # axes2[2,0].plot(
        #         meshes['M1'].cell_center_positions, results["S6"][:,7],
        #         label="RSWME, N = 6", linewidth=2.5, color='darkslategray' , linestyle=':'
        #     )
           

        # axes2[2,0].plot(
        #         meshes['M1'].cell_center_positions, 0.0*results["S1"][:,4],
        #         label="SWME, N = 2", linewidth=2.5,  color='blue'
        #     )

        # axes2[2,0].plot(
        #         meshes['M1'].cell_center_positions, 0*results["S3"][:,5],
        #         label="SWME, N = 4", linewidth=2.5, color='cyan', linestyle='-.'
        #     )
        
        # axes2[2,0].plot(
        #      meshes['M1'].cell_center_positions, results["S5"][:,7],
        #         label="SWME, N = 6", linewidth=2.5,  color='magenta', linestyle=':'
        #     )

        
        # axes2[2, 0].set_ylabel(r'$\alpha_5$', fontsize=16, fontweight='bold')
        # axes2[2, 0].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
        # axes2[2, 0].tick_params(axis='both', which='major', labelsize=14)
        # axes2[2, 0].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
        
        # axes2[2,1].plot(
        #         meshes['M3'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="SWE", linewidth=2.5, color='red'
        #     )
        # axes2[2,1].plot(
        #         meshes['M3'].cell_center_positions, 0.0*results["S2"][:,4],
        #         label="RSWME, N = 2", linewidth=2.5, color='limegreen', linestyle="--"
        #     )
        # axes2[2,1].plot(
        #         meshes['M3'].cell_center_positions, 0.0*results["S4"][:,6],
        #         label="RSWME, N = 4", linewidth=2.5, color='orange', linestyle='-.'
        #     )
        # axes2[2,1].plot(
        #         meshes['M1'].cell_center_positions, results["S6"][:,8],
        #         label="RSWME, N = 6", linewidth=2.5, color='darkslategray' , linestyle=':'
        #     )

        # axes2[2,1].plot(
        #         meshes['M1'].cell_center_positions, 0.0*results["S1"][:,4],
        #         label="SWME, N = 2", linewidth=2.5,  color='blue'
        #     )

        # axes2[2,1].plot(
        #         meshes['M1'].cell_center_positions, 0.0*results["S3"][:,6],
        #         label="SWME, N = 4", linewidth=2.5, color='cyan' , linestyle='-.'
        #     )

        # axes2[2,1].plot(
        #         meshes['M1'].cell_center_positions, results["S5"][:,8],
        #         label="SWME, N = 6", linewidth=2.5,  color='magenta', linestyle=':'
        #     )

       
           
        # axes2[2, 1].set_ylabel(r'$\alpha_6$', fontsize=16, fontweight='bold')
        # axes2[2, 1].set_xlabel(r"$x$", fontsize=16, fontweight='bold')
        # axes2[2, 1].tick_params(axis='both', which='major', labelsize=14)
        # axes2[2, 1].grid(True, linestyle='--', alpha=0.5, linewidth=0.7)
        
        
        # lines, labels = axes2[0, 0].get_legend_handles_labels()
        # rswme_handles = [lines[1], lines[2], lines[3]]
        # rswme_labels = [labels[1], labels[2], labels[3]]

        # swme_handles = [lines[4], lines[5], lines[6]]
        # swme_labels = [labels[4], labels[5], labels[6]]

        # swe_handles = [lines[0]]
        # swe_labels = [labels[0]]
        # #all_handles = rswme_handles + swme_handles + swe_handles
        # #all_labels = rswme_labels + swme_labels + swe_labels
        # #fig2.legend(all_handles, all_labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.08),
        # #    fontsize=16, framealpha=1, edgecolor='black', shadow=True)
        
        # legend1 = fig2.legend(
        #     rswme_handles, rswme_labels,
        #     loc='upper left',
        #     bbox_to_anchor=(0.2, 1.08),
        #     fontsize=16,
        #     framealpha=1,
        #     edgecolor='black',
        #     shadow=True
        # )

        # legend2 = fig2.legend(
        #     swme_handles, swme_labels,
        #     loc='upper center',
        #     bbox_to_anchor=(0.5, 1.08),
        #     fontsize=16,
        #     framealpha=1,
        #     edgecolor='black',
        #     shadow=True
        # )

        # legend3 = fig2.legend(
        #     swe_handles, swe_labels,
        #     loc='upper right',
        #     bbox_to_anchor=(0.80, 1.08),
        #     fontsize=16,
        #     framealpha=1,
        #     edgecolor='black',
        #     shadow=True
        # )


        # axes2[0, 0].add_artist(legend1)
        # axes2[0, 0].add_artist(legend2)
        # axes2[0, 0].add_artist(legend3)

        # plt.tight_layout()
        # plt.savefig("OverviewThirdPlotMomentsIIIb.pdf",format='pdf', bbox_inches='tight')
        #plt.close()

    else:
         print('2D not implemented yet')

  

if __name__ == '__main__':
    main()