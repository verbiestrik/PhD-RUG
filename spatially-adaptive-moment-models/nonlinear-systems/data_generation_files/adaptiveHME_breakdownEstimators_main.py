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

    resolution = 2000
    smooth_pars = [25,50,100,200]
    relaxation_times = [0.05,0.5]
    relaxation_times_str = ["0p05","0p5"]
    tolerances_decrease = [0.0001,0.00025]
    tolerances_increase = [0.0002,0.0005]
    tolerances_decrease_str = ["0p0001","0p00025"]
    tolerances_increase_str = ["0p0002","0p0005"]
    end_times = [0.1,0.25,0.5,0.75]
    end_times_str = ['0p1','0p25','0p5','0p75']

    nr_of_quadrature_points = 3

    _time_integration = timeIntegration.Exact()

    for i in range(len(end_times)):
        for m in range(len(smooth_pars)):
            for n in range(len(relaxation_times)):
                for o in range(len(tolerances_decrease)):
                    for p in range(len(tolerances_increase)):
                        end_time = end_times[i]
                        end_time_str = end_times_str[i]
                        smooth_par = smooth_pars[m]
                        relaxation_time = relaxation_times[n]
                        relaxation_time_str = relaxation_times_str[n]
                        tolerance_decrease = tolerances_decrease[o]
                        tolerance_decrease_str = tolerances_decrease_str[o]
                        tolerance_increase = tolerances_increase[p]
                        tolerance_increase_str = tolerances_increase_str[p]

                        if tolerance_increase > tolerance_decrease:
                            _pde = pde.HermiteMomentEquations(
                                            "smooth_and_shockTube",
                                            relaxation_time,
                                            True,
                                            True,
                                            True,
                                            tolerance_decrease,
                                            tolerance_increase)
                            compute_eigenvalues_and_eigenvectors = _pde.compute_eigenvalues_and_eigenvectors
                            ##########################################################################


                            _spatialDiscretization_interior = spatialDiscretization.Osher(nr_of_quadrature_points,
                                                                                        True,
                                                                                        compute_eigenvalues_and_eigenvectors)
                            _spatialDiscretization_interface = spatialDiscretization.PRICE(nr_of_quadrature_points)                                            
                            _spatialDiscretization_predictor = spatialDiscretization.PRICE(1)                                  

                            _mesh = mesh.UniformRectangularMesh1D([-2.75,3.25], resolution)

                            _simulation = simulation.SmoothedModelAdaptiveSimulationWithInterpolation1D(
                                12,
                                _pde,
                                _mesh,
                                "INFLOW_OUTFLOW",
                                "smooth_and_shockTube",
                                "nothing",
                                smooth_par,
                                _spatialDiscretization_interior,
                                _spatialDiscretization_interface,
                                _spatialDiscretization_predictor,
                                _time_integration) 
                        
                            start = timeit.default_timer()
                            data_array = _simulation.run_simulation(end_time)
                            coarsening_estimator,refinement_estimator = _simulation.get_breakdown_estimators()
                            stop = timeit.default_timer()
                            data_frame = pd.DataFrame(data_array)
                            time = np.full(1,stop-start)
                            data_frame_time = pd.DataFrame(time)
                            data_frame_coarsening_estimator = pd.DataFrame(coarsening_estimator)
                            data_frame_refinement_estimator = pd.DataFrame(refinement_estimator)

                            foldername = 'shock_plus_smooth'+'tau'+relaxation_time_str+'_'+'Osher'                                
                            output_dir = Path(foldername)

                            # Create folders if they don't exist
                            output_dir.mkdir(parents=True, exist_ok=True)

                            filename = 'pred'+'PRICE'+\
                                '_interface'+'PRICE'+\
                                    '_smoothpar'+str(smooth_par)+\
                                    '_toldown'+tolerance_decrease_str+\
                                    '_tolup'+tolerance_increase_str+'_t'+str(end_time_str)+'_Kn'+relaxation_time_str+'.csv'
                            
                            filename_time = 'time'+filename
                            filename_coarsening_estimator = 'estimMin_'+filename
                            filename_refinement_estimator = 'estimPlus_'+filename

                            outputname = output_dir / filename
                            outputname_time = output_dir / filename_time
                            outputname_coarsening_estimator = output_dir / filename_coarsening_estimator
                            outputname_refinement_estimator = output_dir / filename_refinement_estimator
                            data_frame.to_csv(outputname,index=False,header=False)
                            data_frame_time.to_csv(outputname_time,index=False,header=False)
                            data_frame_coarsening_estimator.to_csv(outputname_coarsening_estimator,index=False,header=False)
                            data_frame_refinement_estimator.to_csv(outputname_refinement_estimator,index=False,header=False)

if __name__ == '__main__':
    main()