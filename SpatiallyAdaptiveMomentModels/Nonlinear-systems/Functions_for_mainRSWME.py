import simulation
import pde
import mesh
import spatialDiscretization
import timeIntegration
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import configparser
import timeit
from pathlib import Path
from typing import Dict, List, Any, Tuple


def create_pde(pde_info):
        pde_type = pde_info['pde_type']
        initial_condition = pde_info['initialCondition']
        viscosity = pde_info.getfloat('viscosity')
        slip_length = pde_info.getfloat('slipLength')

        pde_map = {
            'SWME1D': pde.SWME1D(initial_condition, viscosity, slip_length, hyperbolic=False),
            'HSWME1D': pde.SWME1D(initial_condition, viscosity, slip_length, hyperbolic=True),
            'Reduced_SWME1D': pde.Reduced_SWME1D(initial_condition, viscosity, slip_length, hyperbolic=False ),
            'HReduced_SWME1D': pde.Reduced_SWME1D(initial_condition, viscosity, slip_length, hyperbolic=True ),
        }

        if pde_type in pde_map:
            return pde_map[pde_type]
        else:
            print(f'PDE_type "{pde_type}" is not implemented yet')
            return None


def create_spatial(method_info):
        fvm_type = method_info['fvm_type']
        pvm = method_info['pvm']

        if fvm_type == 'PVM':
            pvm_map = {
                'PRICE': spatialDiscretization.PRICE(),
                'LF': spatialDiscretization.LF(),
            }
            if pvm in pvm_map:
                return pvm_map[pvm]
            else:
                print(f'PVM method "{pvm}" is not implemented')
        else:
            print(f'Finite volume type "{fvm_type}" is not implemented')
            return None



def create_time_integration(method_info, pde_info):
        int_type = method_info["timeIntegrator"]
        N = method_info.getint("order")
        nu = float(pde_info['viscosity'])
        lamda = float(pde_info['slipLength'])

        integrator_map = {
            'ImplicitEuler': timeIntegration.ImplicitEuler(),
            'ExplicitEuler': timeIntegration.ExplicitEuler(),
            'ImplicitDirect': timeIntegration.ImplicitDirect(nu, lamda,N),
            'ImplicitFast': timeIntegration.ImplicitFast(nu, lamda,N),
        }

        if int_type in integrator_map:
            return integrator_map[int_type]
        else:
            print(f"Time integrator '{int_type}' is not implemented.")
            return None

def create_mesh(method_info):
    grid_info = method_info["grid_information"]
    return mesh.UniformRectangularMesh1D([grid_info.getfloat('x1boundary'),grid_info.getfloat('x2boundary')],
                                               grid_info.getint('resolutionX'))

def load_configs(config_files: List[str], num_configs: int) -> Dict[str, Dict[str, Any]]:
    return {
        f'C{i}': {
            'pde_information': config['pde_information'],
            'grid_information': config['grid_information'],
            'numerical_method_information': config['numerical_method_information']
        }
        for i, file in enumerate(config_files, start=1)
        if (config := configparser.ConfigParser()).read(file)
    }

def create_pdes(configs: Dict[str, Dict[str, Any]], num_configs: int) -> Dict[str, Any]:
    return {
        f'pde{i}': create_pde(configs[f'C{i}']['pde_information'])
        for i in range(1, num_configs + 1)
    }

def create_spatial_discretizations(configs: Dict[str, Dict[str, Any]], num_configs: int) -> Dict[str, Any]:
    return {
        f'spatialDiscretization{i}': create_spatial(configs[f'C{i}']['numerical_method_information']) 
        for i in range(1, num_configs + 1)
    }

def create_time_integrations(configs: Dict[str, Dict[str, Any]], num_configs:int) -> Dict[str, Any]:
    return {
        f'time_integration{i}': create_time_integration(configs[f'C{i}']['numerical_method_information'], configs[f'C{i}']["pde_information"])
        for i in range(1, num_configs + 1)
    }


def extract_is_1d(configs: Dict[str, Dict[str, Any]], num_configs:int, key:str)-> List[bool]:
    return [
        configs[f'C{i}']['pde_information'].getboolean(key)
        for i in range(1, num_configs + 1)
    ]

def extract_is_classical(configs: Dict[str, Dict[str, Any]], num_configs:int, key:str) -> List[bool]:
    return [
        configs[f'C{i}']['numerical_method_information']['method'] == key
        for i in range(1, num_configs + 1)
    ]

def create_meshes(configs: Dict[str, Dict[str, Any]], num_configs:int) -> Dict[str, Any]:
    return {
        f'M{i}': create_mesh(configs[f'C{i}'])
        for i in range(1, num_configs + 1)
    }

def create_simulations(configs: Dict[str, Dict[str, Any]], pdes: Dict[str, Any],
                      meshes: Dict[str, Any], spatial_discretizations: Dict[str, Any],
                      time_integrations: Dict[str, Any], num_configs: int) -> Dict[str, Any]:
    return {
        f"S{i}": simulation.ClassicalSimulation1D(
            configs[f"C{i}"]["numerical_method_information"].getint('order'),
            pdes[f"pde{i}"],
            meshes[f"M{i}"],
            configs[f"C{i}"]["numerical_method_information"]['boundaryCondition'],
            configs[f"C{i}"]['pde_information']['initialCondition'],
            spatial_discretizations[f"spatialDiscretization{i}"],
            time_integrations[f"time_integration{i}"],
            configs[f"C{i}"]['pde_information'].getfloat('slipLength'),
            configs[f"C{i}"]['pde_information'].getfloat('viscosity'),
            configs[f"C{i}"]['pde_information']['pde_type']
        )
        for i in range(1, num_configs + 1)
    }

def run_simulations(configs: Dict[str, Dict[str,Any]], simulations:Dict[str, Any], num_configs:int) -> Tuple[Dict[str, Any], Dict[str, float]]:
    results, timings = {}, {}
    for i in range(1, num_configs + 1):
        start = timeit.default_timer()
        data_array = simulations[ f"S{i}"].run_simulation(configs[f"C{i}"]["numerical_method_information"].getfloat('t_end'))

        timings[ f"S{i}"] = timeit.default_timer() - start
        results[ f"S{i}"] = data_array
    return results, timings


def Godunovmin(h: np.array, dx: float) -> np.array:
    n = len(h)

    hg = np.zeros(n + 2)
    hg[1:-1] = h
    hg[0] = h[-1]    
    hg[-1] = h[0]    

    Dh = np.zeros(h.shape)

    for i in range(1, n + 1):
        if -hg[i] >= 0 and -hg[i + 1] >= 0:
            F_right = -hg[i]**4
        elif -hg[i] <= 0 and -hg[i + 1] <= 0:
            F_right = -hg[i + 1]**4
        else:
            F_right = 0 

        if -hg[i - 1] >= 0 and - hg[i] >= 0:
            F_left = -hg[i - 1]**4
        elif -hg[i - 1] <= 0 and -hg[i] <= 0:
            F_left = -hg[i]**4
        else:
            F_left = 0  # Sonic point

        Dh[i - 1] = (F_right - F_left) / dx

    return Dh


