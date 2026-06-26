import simulation
import pde
import mesh
import spatialDiscretization
import timeIntegration
import plotting
import pandas as pd
import configparser
import timeit
from pathlib import Path
import numpy as np

def main():

    config = configparser.ConfigParser()
    config.read('Config-files/config.txt')
    pde_information = config['pde_information']
    grid_information = config['grid_information']
    numerical_method_information = config['numerical_method_information']
    adaptive_simulation_information = config['adaptive_simulation_information']

    linear_source = pde_information.getboolean('linear_source')
    time_integrator = numerical_method_information['timeIntegrator']
    CFL_number = numerical_method_information.getfloat('CFL')
    linear_source_implicit = linear_source and time_integrator == 'ImplicitEuler'
    exact_source_computation = time_integrator == 'Exact'

    eigenstructure_available = False
    compute_eigenvalues_and_eigenvectors = None

    type_model_error_estimator = adaptive_simulation_information['type_model_error_estimator']
    time_step_splitting = adaptive_simulation_information.getboolean('time_step_splitting')

    tol_coarsening_model_difference = adaptive_simulation_information.getfloat('tol_coarsening_model_difference')
    tol_refinement_model_difference = adaptive_simulation_information.getfloat('tol_refinement_model_difference')

    tol_coarsening_heur = adaptive_simulation_information.getfloat('tol_coarsening_heur')
    tol_refinement_heur_source = adaptive_simulation_information.getfloat('tol_refinement_heur_source')
    tol_refinement_heur_grad = adaptive_simulation_information.getfloat('tol_refinement_heur_grad') 

    tol_res1 = adaptive_simulation_information.getfloat('tol_res1')
    tol_res2 = adaptive_simulation_information.getfloat('tol_res2')    

    if type_model_error_estimator == 'model_difference':
        tols_coarsening = [tol_coarsening_model_difference]
        tols_refinement = [tol_refinement_model_difference]
    elif type_model_error_estimator == 'heuristics_plus_discretization':
        tols_coarsening = [tol_res1,tol_res2,tol_coarsening_heur]
        tols_refinement = [tol_refinement_heur_source,tol_refinement_heur_grad,tol_refinement_heur_grad]

    smoothing = adaptive_simulation_information.getboolean('smoothing')
    smooth_par = adaptive_simulation_information.getint('smooth_par')
    path_conservation = adaptive_simulation_information.getboolean('path_conservation')
    start_order = adaptive_simulation_information.getint('start_order')
    interpolation = adaptive_simulation_information.getboolean('interpolation')
    spatial_discretization_predictor = adaptive_simulation_information['spatial_discretization_predictor']
    spatial_discretization_interface = adaptive_simulation_information['spatial_discretization_interface']
    two_step_domain_decomposition_evaluation = adaptive_simulation_information.getboolean('two_step_domain_decomposition_evaluation')
    hierarchical = adaptive_simulation_information['hierarchical']

    min_order = 0
    order_diff = 1

    if pde_information['pde_type'] == 'SWME1D':
        min_order = 0
        order_diff = 1
        _pde = pde.SWME1D(pde_information['initialCondition'],
                        pde_information.getfloat('viscosity'),
                        pde_information.getfloat('slipLength'),
                        False,
                        linear_source_implicit,
                        type_model_error_estimator)
    elif pde_information['pde_type'] == 'HSWME1D':
        min_order = 0
        order_diff = 1
        _pde = pde.SWME1D(pde_information['initialCondition'],
                        pde_information.getfloat('viscosity'),
                        pde_information.getfloat('slipLength'),
                        True,
                        linear_source_implicit,
                        type_model_error_estimator)
        
    elif pde_information['pde_type'] == 'VegetationSWME1D':
        min_order = 0
        order_diff = 1
        _pde = pde.VegetationSWME1D(pde_information['initialCondition'],
                                pde_information.getfloat('viscosity'),
                                pde_information.getfloat('slipLength'),
                                False,
                                linear_source_implicit,
                                type_model_error_estimator,
                                pde_information.getfloat('diameter'),
                                pde_information.getfloat('CD'),
                                pde_information.getfloat('surface_density'),
                                pde_information.getfloat('h_v'))
    elif pde_information['pde_type'] == 'HME':
        min_order = 2
        order_diff = 2
        _pde = pde.HermiteMomentEquations(
                        pde_information['initialCondition'],
                        pde_information.getfloat('relaxation_time'),
                        True,
                        True,
                        exact_source_computation,
                        type_model_error_estimator)
        eigenstructure_available = True
        compute_eigenvalues_and_eigenvectors = _pde.compute_eigenvalues_and_eigenvectors
    elif pde_information['pde_type'] == 'Grad':
        min_order = 2
        order_diff = 2
        _pde = pde.HermiteMomentEquations(
                        pde_information['initialCondition'],
                        pde_information.getfloat('relaxation_time'),
                        False,
                        True,
                        exact_source_computation,
                        type_model_error_estimator)
    else:
        print('PDE_type is not implemented yet')

    ##########################################################################

    if numerical_method_information['fvm_type'] == 'PVM':
        nr_of_quadrature_points = numerical_method_information.getint('nr_of_quadrature_points')
        if numerical_method_information['pvm'] == 'PRICE':
            _spatialDiscretization = spatialDiscretization.PRICE(nr_of_quadrature_points)
        elif numerical_method_information['pvm'] == 'LF':
            _spatialDiscretization = spatialDiscretization.LF(nr_of_quadrature_points)
        elif numerical_method_information['pvm'] == 'Roe':
            _spatialDiscretization = spatialDiscretization.Roe(nr_of_quadrature_points)
        elif numerical_method_information['pvm'] == 'Osher':
            _spatialDiscretization = spatialDiscretization.Osher(nr_of_quadrature_points,eigenstructure_available,compute_eigenvalues_and_eigenvectors)
        else:
            print('this pvm method is not implemented yet')

        if spatial_discretization_predictor == 'PRICE':
            _spatialDiscretizationPredictor = spatialDiscretization.PRICE(nr_of_quadrature_points)
        elif spatial_discretization_predictor == 'LF':
            _spatialDiscretizationPredictor = spatialDiscretization.LF(nr_of_quadrature_points)
        elif spatial_discretization_predictor == 'Roe':
            _spatialDiscretizationPredictor = spatialDiscretization.Roe(nr_of_quadrature_points)
        elif spatial_discretization_predictor == 'Osher':
            _spatialDiscretizationPredictor = spatialDiscretization.Osher(nr_of_quadrature_points,eigenstructure_available,compute_eigenvalues_and_eigenvectors)
        else:
            print('this pvm method is not implemented yet')

        if spatial_discretization_interface == 'PRICE':
            _spatialDiscretizationInterface = spatialDiscretization.PRICE(nr_of_quadrature_points)
        elif spatial_discretization_interface == 'LF':
            _spatialDiscretizationInterface = spatialDiscretization.LF(nr_of_quadrature_points)
        elif spatial_discretization_interface == 'Roe':
            _spatialDiscretizationInterface = spatialDiscretization.Roe(nr_of_quadrature_points)
        elif spatial_discretization_interface == 'Osher':
            _spatialDiscretizationInterface = spatialDiscretization.Osher(nr_of_quadrature_points,eigenstructure_available,compute_eigenvalues_and_eigenvectors)
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


        if numerical_method_information['method'] == 'modelAdaptiveMomentSimulation1D':
            
            start_order = int(adaptive_simulation_information['start_order'])
            _simulation = simulation.ModelAdaptiveMomentSimulation1D(
            start_order,
            min_order,
            _pde,
            _mesh,
            CFL_number,
            numerical_method_information['boundaryCondition'],
            pde_information['initialCondition'],
            smoothing,
            smooth_par,
            interpolation,
            path_conservation,
            _spatialDiscretization,
            _spatialDiscretizationInterface,
            _spatialDiscretizationPredictor,
            _time_integration,
            two_step_domain_decomposition_evaluation,
            hierarchical,
            type_model_error_estimator,
            time_step_splitting,
            order_diff,
            tols_coarsening,
            tols_refinement) 
        elif numerical_method_information['method'] == 'classical':
            _simulation = simulation.ClassicalSimulation1D(
                numerical_method_information.getint('order'),
                _pde,
                _mesh,
                numerical_method_information['boundaryCondition'],
                pde_information['initialCondition'],
                _spatialDiscretization,
                _time_integration,
                CFL_number)
            
        elif numerical_method_information['method'] == 'micro_macro':
            _simulation = simulation.Micro_macro(
                [int(order) for order in numerical_method_information['orders'].split(',')],
                _pde,
                _mesh,
                numerical_method_information['boundaryCondition'],
                pde_information['initialCondition'],
                _spatialDiscretization,
                _time_integration)

        if pde_information['pde_type'] == 'SWME1D' or pde_information['pde_type'] == 'HSWME1D'\
            or pde_information['pde_type'] == 'VegetationSWME1D':
            if numerical_method_information['method'] == 'modelAdaptiveMomentSimulation1D':
                _plotting = plotting.SWME1DPlotAdaptive(_pde,_mesh,_simulation,type_model_error_estimator)
            elif numerical_method_information['method'] == 'classical' or\
                    numerical_method_information['method'] == 'micro_macro':
                _plotting = plotting.SWME1DPlotClassical(_pde,_mesh,_simulation)
        elif pde_information['pde_type'] == 'HME' or pde_information['pde_type'] == 'Grad':
            if numerical_method_information['method'] == 'modelAdaptiveMomentSimulation1D':
                _plotting = plotting.HME1DPlotAdaptive(_pde,_mesh,_simulation)
            elif numerical_method_information['method'] == 'classical' or\
                    numerical_method_information['method'] == 'micro_macro':
                _plotting = plotting.HME1DPlotClassical(_pde,_mesh,_simulation)
    
        start = timeit.default_timer()
        data_array = _simulation.run_simulation(numerical_method_information.getfloat('t_end'))
        stop = timeit.default_timer()
        time = np.full(1,stop-start)
        print('Time: ', time)
        data_frame = pd.DataFrame(data_array)
        _plotting.plot(data_array)
        data_frame_time = pd.DataFrame(time)

        output_information = config['output_information']
        if output_information.getboolean('export_data'):

            base_folder = Path("Data-processing/Output")

            foldername = output_information['output_folder_name']                               
            output_dir = base_folder / foldername

            # Create folders if they don't exist
            output_dir.mkdir(parents=True, exist_ok=True)

            # filename = 'smoothPlusDam'+'.csv'
            filename = output_information['output_file_name']+'.csv'
            
            filename_time = 'time_'+filename

            outputname = output_dir / filename
            outputname_time = output_dir / filename_time
            data_frame.to_csv(outputname,index=False,header=False)
            data_frame_time.to_csv(outputname_time,index=False,header=False)
    else:
        print('2D not implemented yet')

if __name__ == '__main__':
    main()