import simulation
import pde
import mesh
import spatialDiscretization
import timeIntegration
import pandas as pd
import timeit
import numpy as np
from pathlib import Path

def main():

    spatial_discretizations_interior = ["Osher","PRICE"]
    spatial_discretizations_predictor = ["LF","PRICE"]
    spatial_discretizations_interface = ["Osher","PRICE"]
    resolutions = [1000,6000]
    smooth_pars = [5,10,20,40]
    relaxation_times = [0.05,0.5]
    relaxation_times_str = ["0p05","0p5"]
    tolerances_decrease = [0.0001,0.001,0.01]
    tolerances_increase = [0.0002,0.002,0.02]
    tolerances_decrease_str = ["0p0001","0p001","0p01"]
    tolerances_increase_str = ["0p0002","0p002","0p02"]

    _time_integration = timeIntegration.Exact()

    for i in range(len(spatial_discretizations_interior)):
        for j in range(len(spatial_discretizations_predictor)):
            for k in range(len(spatial_discretizations_interface)):
                for m in range(len(smooth_pars)):
                    for n in range(len(relaxation_times)):
                        for o in range(len(tolerances_decrease)):
                            for p in range(len(tolerances_increase)):
                                spatial_discretization_interior = spatial_discretizations_interior[i]
                                resolution = resolutions[i]
                                spatial_discretization_predictor = spatial_discretizations_predictor[j]
                                spatial_discretization_interface = spatial_discretizations_interface[k]
                                smooth_par = smooth_pars[m]
                                relaxation_time = relaxation_times[n]
                                relaxation_time_str = relaxation_times_str[n]
                                tolerance_decrease = tolerances_decrease[o]
                                tolerance_decrease_str = tolerances_decrease_str[o]
                                tolerance_increase = tolerances_increase[p]
                                tolerance_increase_str = tolerances_increase_str[p]

                                if tolerance_increase > tolerance_decrease:
                                    if not (spatial_discretization_interior == 'PRICE' and spatial_discretization_interface == 'Osher'):
                                        _pde = pde.HermiteMomentEquations(
                                                        "shockTube_noVelocity",
                                                        relaxation_time,
                                                        True,
                                                        True,
                                                        True,
                                                        tolerance_decrease,
                                                        tolerance_increase)
                                        compute_eigenvalues_and_eigenvectors = _pde.compute_eigenvalues_and_eigenvectors
                                        ##########################################################################

                                        nr_of_quadrature_points = 3
                                        if spatial_discretization_interior == 'PRICE':
                                            _spatialDiscretization_interior = spatialDiscretization.PRICE(nr_of_quadrature_points)
                                            _spatialDiscretization_interface = spatialDiscretization.PRICE(nr_of_quadrature_points)
                                        elif spatial_discretization_interior == 'Osher':
                                            _spatialDiscretization_interior = spatialDiscretization.Osher(nr_of_quadrature_points,
                                                                                                        True,
                                                                                                        compute_eigenvalues_and_eigenvectors)
                                            if spatial_discretization_interface == 'Osher':
                                                _spatialDiscretization_interface = spatialDiscretization.Osher(nr_of_quadrature_points,
                                                                                                            True,
                                                                                                            compute_eigenvalues_and_eigenvectors)
                                            elif spatial_discretization_interface == 'PRICE':
                                                _spatialDiscretization_interface = spatialDiscretization.PRICE(nr_of_quadrature_points)                                            
                                        if spatial_discretization_predictor == 'PRICE':
                                            _spatialDiscretization_predictor = spatialDiscretization.PRICE(nr_of_quadrature_points)
                                        elif spatial_discretization_predictor == 'LF':
                                            _spatialDiscretization_predictor = spatialDiscretization.LF(nr_of_quadrature_points)                                        



                                        _mesh = mesh.UniformRectangularMesh1D([-1.5,1.5], resolution)

                                        _simulation = simulation.SmoothedModelAdaptiveSimulation1D(
                                            12,
                                            _pde,
                                            _mesh,
                                            "INFLOW_OUTFLOW",
                                            "shockTube_noVelocity",
                                            "nothing",
                                            smooth_par,
                                            _spatialDiscretization_interior,
                                            _spatialDiscretization_interface,
                                            _spatialDiscretization_predictor,
                                            _time_integration) 
                                    
                                        start = timeit.default_timer()
                                        data_array = _simulation.run_simulation(0.3)
                                        stop = timeit.default_timer()
                                        data_frame = pd.DataFrame(data_array)
                                        time = np.full(1,stop-start)
                                        data_frame_time = pd.DataFrame(time)

                                        foldername = 'tau'+relaxation_time_str+'_'+spatial_discretization_interior                                
                                        output_dir = Path(foldername)

                                        # Create folders if they don't exist
                                        output_dir.mkdir(parents=True, exist_ok=True)

                                        filename = 'pred'+spatial_discretization_predictor+\
                                            '_interface'+spatial_discretization_interface+\
                                                '_smoothpar'+str(smooth_par)+\
                                                '_toldown'+tolerance_decrease_str+\
                                                '_tolup'+tolerance_increase_str+'.csv'
                                        
                                        filename_time = 'time'+filename

                                        outputname = output_dir / filename
                                        outputname_time = output_dir / filename_time
                                        data_frame.to_csv(outputname,index=False,header=False)
                                        data_frame_time.to_csv(outputname_time,index=False,header=False)


if __name__ == '__main__':
    main()