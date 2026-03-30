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
import csv


def main():
    #print("withimplicitFast!")
    NUM_CONFIGS = 12
    IS_1D_KEY = '1D'
    IS_CLASSICAL_KEY = 'classical'
    CONFIG_DIR = Path("Config-filesH")
    
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
            results, times = run_simulations(configs, simulations, NUM_CONFIGS)
        else:
            "Error no classical methods"
        
        data_sets = [
            (1 / configs[f"C{i}"]['pde_information'].getfloat('slipLength'),
            results[f"S{i}"], results[f"S{i+1}"], results[f"S{i+2}"], results[f"S{i+3}"])
            for i in range(1, NUM_CONFIGS+1, 4)
        ]

        
        for i in range(1,NUM_CONFIGS+1):
            headers = ['x', 'h', 'um', 'a1', 'a2'] 
            filename = f"Numerical_Results/Test1_Hyperbolicity/EFastdata_set{i}.csv"
            eps_val =1/(float(configs[f"C{i}"]['pde_information']['slipLength']))
            time = times[f"S{i}"]
            order = configs[f"C{i}"]["numerical_method_information"].getfloat('order')
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
        results_list = []
        for eps, ref, swe, rswe, hrswe in data_sets:
            filename = f"Numerical_Results/Test1_Hyperbolicity/Errors_eps_{eps}.csv"
            swe_h = np.linalg.norm(ref[:, 1] - swe[:, 1], 1)/np.linalg.norm(ref[:, 1] , 1)
            swe_u = np.linalg.norm(ref[:, 2] - swe[:, 2], 1)/np.linalg.norm(ref[:, 2] , 1)
            rswe_h = np.linalg.norm(ref[:, 1] - rswe[:, 1], 1)/np.linalg.norm(ref[:, 1] , 1)
            rswe_u = np.linalg.norm(ref[:, 2] - rswe[:, 2], 1)/np.linalg.norm(ref[:, 2] , 1)
            hrswe_h = np.linalg.norm(ref[:, 1] - hrswe[:, 1], 1)/np.linalg.norm(ref[:, 1] , 1)
            hrswe_u = np.linalg.norm(ref[:, 2] - hrswe[:, 2], 1)/np.linalg.norm(ref[:, 2] , 1)

            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([f"for eps = {eps}"])
                writer.writerow([f"Error for SWE in h: {swe_h}"])
                writer.writerow([f"Error for SWE in u: {swe_u}"])
                writer.writerow([f"Error for RSWE in h: {rswe_h}"])
                writer.writerow([f"Error for RSWE in u: {rswe_u}"])
                writer.writerow([f"Error for HRSWE in h: {hrswe_h}"])
                writer.writerow([f"Error for HRSWE in u: {hrswe_u}"])
           

    else:
        print('2D not implemented yet')

  

if __name__ == '__main__':
    main()