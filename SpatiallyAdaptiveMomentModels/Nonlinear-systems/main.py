# Packages & Local Imports
import simulation
import pde
import mesh
import spatialDiscretization
import timeIntegration
import plotting
import pandas as pd
import configparser
import timeit

# Recharge Specific Imports. Everything is written as a try-except block to
# avoid import errors in case the recharge module never merges with the main
# branch. 
try: 
    import os
    from recharge.recharge_pde import RechargeSWME1D
    from recharge.laws import HortonInfiltration

    # Define a global boolean flag that establishes if the recharge module 
    # was introduced at runtime.
    HAS_RECHARGE = True

    # Make a post-processing directory if one doesn't exist. 
    # Compartmentalize the recharge results in a separate folder but in the same
    # Results/ directory.
    os.makedirs('Data-processing/Results/Recharge', exist_ok=True)

except ImportError:
    # Silent fail, don't print an if statement
    HAS_RECHARGE = False

def main():

    config = configparser.ConfigParser()
    config.read('Config-files/config.txt')
    pde_information = config['pde_information']
    grid_information = config['grid_information']
    numerical_method_information = config['numerical_method_information']

    linear_source = pde_information.getboolean('linear_source')
    time_integrator = numerical_method_information['timeIntegrator']
    linear_source_implicit = linear_source and time_integrator == 'ImplicitEuler'
    exact_source_computation = time_integrator == 'Exact'

    if pde_information['pde_type'] == 'SWME1D':
        _pde = pde.SWME1D(pde_information['initialCondition'],
                        pde_information.getfloat('viscosity'),
                        pde_information.getfloat('slipLength'),
                        False,
                        linear_source_implicit)
    elif pde_information['pde_type'] == 'HSWME1D':
        _pde = pde.SWME1D(pde_information['initialCondition'],
                        pde_information.getfloat('viscosity'),
                        pde_information.getfloat('slipLength'),
                        True,
                        linear_source_implicit)
        
    elif pde_information['pde_type'] == 'VegetationSWME1D':
        _pde = pde.VegetationSWME1D(pde_information['initialCondition'],
                                pde_information.getfloat('viscosity'),
                                pde_information.getfloat('slipLength'),
                                False,
                                linear_source_implicit,
                                0.008,
                                0.97,
                                800,
                                0.4)
    elif pde_information['pde_type'] == 'HME':
        _pde = pde.HermiteMomentEquations(
                        pde_information['initialCondition'],
                        pde_information.getfloat('relaxation_time'),
                        True,
                        True,
                        exact_source_computation)
    elif pde_information['pde_type'] == 'Grad':
        _pde = pde.HermiteMomentEquations(
                        pde_information['initialCondition'],
                        pde_information.getfloat('relaxation_time'),
                        False,
                        True,
                        exact_source_computation)
    elif pde_information['pde_type'] == 'RechargeSWME1D':
        if not HAS_RECHARGE:
            raise ImportError(
                "Config requests pde_type='RechargeSWME1D' but the recharge module  "
                "is not available in this branch/environment."
            )
        
        # Build the infiltration object
        infiltration_model = HortonInfiltration(
            f0 = pde_information.getfloat('horton_f0'),
            fc = pde_information.getfloat('horton_fc'),
            k = pde_information.getfloat('horton_k'),
        )

        # Build the pde object
        _pde = RechargeSWME1D(
            pde_information['initialCondition'],
            pde_information.getfloat('viscosity'),
            pde_information.getfloat('slipLength'),
            pde_information.getboolean('hyperbolic', fallback=False),    
            pde_information.getboolean('linear_source', fallback=False),
            pde_information.getfloat('rainfall_rate'),
            infiltration_model,
            pde_information.getfloat('f_R', fallback=1.0),
            pde_information.getfloat('f_I', fallback=0.0),
            )
    else:
        print('PDE_type is not implemented yet')
    
    ##########################################################################

    if numerical_method_information['fvm_type'] == 'PVM':
        if numerical_method_information['pvm'] == 'PRICE':
            _spatialDiscretization = spatialDiscretization.PRICE()
        elif numerical_method_information['pvm'] == 'LF':
            _spatialDiscretization = spatialDiscretization.LF()
        elif numerical_method_information['pvm'] == 'Roe':
            _spatialDiscretization = spatialDiscretization.Roe()
        elif numerical_method_information['pvm'] == 'Osher':
            _spatialDiscretization = spatialDiscretization.Osher()
        else:
            print('this pvm method is not implemented yet')
    else:
        print('this finite volume type is not implemented yet')

    if numerical_method_information['timeIntegrator'] == 'ImplicitEuler':
        _time_integration = timeIntegration.ImplicitEuler(linear_source)
    elif numerical_method_information['timeIntegrator'] == 'ExplicitEuler':
        _time_integration = timeIntegration.ExplicitEuler()
    elif numerical_method_information['timeIntegrator'] == 'Exact':
        _time_integration = timeIntegration.Exact()

    #########################################################################

    if pde_information.getboolean('1D'):

        _mesh = mesh.UniformRectangularMesh1D([grid_information.getfloat('x1boundary'),grid_information.getfloat('x2boundary')],
                                               grid_information.getint('resolutionX')) #TODO: Implement different grids

        if numerical_method_information['method'] == 'spatially_adaptive':
            start_order = int(numerical_method_information['start_order'])

            if numerical_method_information['coupling'] == 'nonconservative':
                _simulation = simulation.NonConservativeAdaptiveSimulation1D(
                    start_order,
                    _pde,
                    _mesh,
                    numerical_method_information['boundaryCondition'],
                    pde_information['initialCondition'],
                    pde_information['breakdown_criterion'],
                    _spatialDiscretization,
                    _time_integration
                )
            elif numerical_method_information['coupling'] == 'conservative':
                _simulation = simulation.ConservativeAdaptiveSimulation1D(
                    start_order,
                    _pde,
                    _mesh,
                    numerical_method_information['boundaryCondition'],
                    pde_information['initialCondition'],
                    pde_information['breakdown_criterion'],
                    _spatialDiscretization,
                    _time_integration
                )
        elif numerical_method_information['method'] == 'smoothedAdaptive':
            if numerical_method_information['coupling'] == 'nonconservative':
                start_order = int(numerical_method_information['start_order'])
                _simulation = simulation.SmoothedConsAdaptiveSimulation1D(
                    start_order,
                    _pde,
                    _mesh,
                    numerical_method_information['boundaryCondition'],
                    pde_information['initialCondition'],
                    pde_information['breakdown_criterion'],
                    _spatialDiscretization,
                    _time_integration)
            elif numerical_method_information['coupling'] == 'nonconservative':
                start_order = int(numerical_method_information['start_order'])
                _simulation = simulation.SmoothedNonConsAdaptiveSimulation1D(
                    start_order,
                    _pde,
                    _mesh,
                    numerical_method_information['boundaryCondition'],
                    pde_information['initialCondition'],
                    pde_information['breakdown_criterion'],
                    _spatialDiscretization,
                    _time_integration)                
        elif numerical_method_information['method'] == 'interpolatedAdaptive':
            start_order = int(numerical_method_information['start_order'])
            _simulation = simulation.InterpolatedAdaptiveSimulation1D(
                start_order,
                _pde,
                _mesh,
                numerical_method_information['boundaryCondition'],
                pde_information['initialCondition'],
                pde_information['breakdown_criterion'],
                _spatialDiscretization,
                _time_integration) 
        elif numerical_method_information['method'] == 'classical':
            _simulation = simulation.ClassicalSimulation1D(
                numerical_method_information.getint('order'),
                _pde,
                _mesh,
                numerical_method_information['boundaryCondition'],
                pde_information['initialCondition'],
                _spatialDiscretization,
                _time_integration)
            
        elif numerical_method_information['method'] == 'micro_macro':
            _simulation = simulation.Micro_macro(
                [int(order) for order in numerical_method_information['orders'].split(',')],
                _pde,
                _mesh,
                numerical_method_information['boundaryCondition'],
                pde_information['initialCondition'],
                _spatialDiscretization,
                _time_integration)
        
        # Modify the old plotting code to be slighly more readable.
        swme_plot_types = ['SWME1D', 'HSWME1D']
        if HAS_RECHARGE:
            swme_plot_types.append('RechargeSWME1D')

        if pde_information['pde_type'] in swme_plot_types:
            if numerical_method_information['method'] in [
                'spatially_adaptive',
                'smoothedAdaptive',
                'interpolatedAdaptive'
            ]:
                _plotting = plotting.SWME1DPlotAdaptive(_pde, _mesh, _simulation)
            elif numerical_method_information['method'] == 'classical':
                _plotting = plotting.SWME1DPlotClassical(_pde, _mesh, _simulation)
        elif pde_information['pde_type'] in ['HME', 'Grad']:
            if numerical_method_information['method'] in [
                'spatially_adaptive',
                'smoothedAdaptive',
                'interpolatedAdaptive'
            ]:
                _plotting = plotting.HME1DPlotAdaptive(_pde, _mesh, _simulation)
            elif numerical_method_information['method'] == 'classical':
                _plotting = plotting.HME1DPlotClassical(_pde, _mesh, _simulation)
    
        start = timeit.default_timer()
        if (
            HAS_RECHARGE 
            and pde_information['pde_type'] == 'RechargeSWME1D'
            and numerical_method_information['method'] == 'classical'
        ):
            # Store the history in the simulation object
            _simulation.store_history = True
            # Store every 1 time step. Adjust to larger values to reduce storage.
            _simulation.history_stride = 1

        data_array = _simulation.run_simulation(numerical_method_information.getfloat('t_end'))

        # Recharge specific post-processing
        if HAS_RECHARGE and pde_information['pde_type'] == 'RechargeSWME1D':
            final_df = pd.DataFrame(data_array, columns=["x", "h", "u_m", "a1"])
            final_df.to_csv(
                "Data-processing/Results/Recharge/recharge_results.csv",
                index=False,
            )

            if hasattr(_simulation, "history") and len(_simulation.history) > 0:
                field_rows = []
                summary_rows = []

                for entry in _simulation.history:
                    step = entry["step"]
                    time = entry["time"]
                    snapshot = entry["data"]

                    for row in snapshot:
                        field_rows.append({
                            "step": step,
                            "time": time,
                            "x": row[0],
                            "h": row[1],
                            "u_m": row[2],
                            "a1": row[3],
                        })

                    summary_rows.append({
                        "step": step,
                        "time": time,
                        "mean_h": snapshot[:, 1].mean(),
                        "mean_u_m": snapshot[:, 2].mean(),
                        "mean_a1": snapshot[:, 3].mean(),
                        "min_h": snapshot[:, 1].min(),
                        "max_h": snapshot[:, 1].max(),
                    })
                
                pd.DataFrame(field_rows).to_csv(
                    "Data-processing/Results/Recharge/recharge_field_history.csv",
                    index=False,
                )
                pd.DataFrame(summary_rows).to_csv(
                    "Data-processing/Results/Recharge/recharge_summary_history.csv",
                    index=False,
                )

        stop = timeit.default_timer()
        print('Time: ', stop - start)
        data_frame = pd.DataFrame(data_array)

        # Making the plotting call safe
        if '_plotting' in locals():
            _plotting.plot(data_array)
        else: 
            print("No plotting class defined for pde_type =", 
                  pde_information['pde_type'])
        # data_frame.to_csv('Data-processing/Output/test.csv', index=False,header=False)
        # data_frame.to_csv(
        #     'Data-processing/Results/KineticMomentEquations/smoothAndShockTube_order10_relaxation0.1_time1.0_3000.csv',
        #     index=False,
        #     header=False)
    else:
        print('2D not implemented yet')

if __name__ == '__main__':
    main()