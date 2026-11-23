import simulation_HSWMEDataGeneration
import pde_HSWMEDataGeneration
import mesh
import spatialDiscretization
import timeIntegration
import pandas as pd
import numpy as np
import timeit
from pathlib import Path

def main_SWME_errorData():

    viscosities = [0.05,0.5,5.0]
    slip_lengths = [0.05,0.5,5.0]
    orders = [0,1,2,3,4,5,6]
    initialHeight_names = ["hDamLarge","hDamMid", "hDamSmall", "hSmoothSteep", "hSmoothMid","hSmoothSmall"]
    initialVelocity_names = ["uConst", "uVarXnotvarZ", "uVarZnotvarX","uVarXvarZ"]
    initialVelocityMagnitude_names = ["uFast","uMid","uSlow"]
    num_scheme_names = ["Roe","PRICE"]
    for i in range(len(viscosities)):
        for j in range(len(slip_lengths)):
            for k in range(len(orders)):
                for m in range(len(initialHeight_names)):
                    for n in range(len(initialVelocity_names)):
                        for o in range(len(initialVelocityMagnitude_names)):
                            for p in range(len(num_scheme_names)):
                                linear_source = True
                                time_integrator = "ImplicitEuler"
                                linear_source_implicit = linear_source and time_integrator == 'ImplicitEuler'
                                exact_source_computation = time_integrator == 'Exact'

                                _pde = pde_HSWMEDataGeneration.SWME1D(initialHeight_names[m],
                                                                    initialVelocity_names[n],
                                                                    initialVelocityMagnitude_names[o],
                                                                    viscosities[i],
                                                                    slip_lengths[j],
                                                                    True,
                                                                    linear_source_implicit)
                                
                                if num_scheme_names[p] == "Roe":
                                    _spatialDiscretization = spatialDiscretization.Roe()
                                elif num_scheme_names[p] == "PRICE":
                                    _spatialDiscretization = spatialDiscretization.PRICE()

                                _time_integration = timeIntegration.ImplicitEuler(linear_source)

                                _mesh = mesh.UniformRectangularMesh1D([-1,1],2000)

                                _simulation = simulation_HSWMEDataGeneration.ClassicalSimulation1D(
                                    orders[k],
                                    _pde,
                                    _mesh,
                                    "INFLOW_OUTFLOW",
                                    initialHeight_names[m],
                                    initialVelocity_names[n],
                                    initialVelocityMagnitude_names[o],
                                    _spatialDiscretization,
                                    _time_integration)
                                
                                start = timeit.default_timer()
                                data_array,covariates = _simulation.run_simulation(0.2)
                                stop = timeit.default_timer()
                                # print('Time: ', stop - start)
                                covariates[16] = viscosities[i]
                                covariates[17] = slip_lengths[j]
                                data_frame = pd.DataFrame(data_array)
                                covariates_frame = pd.DataFrame(covariates)

                                foldername = initialHeight_names[m]+'_'\
                                    +initialVelocity_names[n]+'_'+initialVelocityMagnitude_names[o]
                                
                                base_dir = Path('error_trainingData')
                                output_dir = base_dir / foldername

                                # Create folders if they don't exist
                                output_dir.mkdir(parents=True, exist_ok=True)

                                filename = 'lambda'+str(slip_lengths[j])+'_viscosity'+str(viscosities[i])+'_order'\
                                    +str(orders[k])+'_FVM'+num_scheme_names[p]+'_2000cells'+'.csv'
                                filename_covs = 'covs_'+filename

                                outputname = output_dir / filename
                                data_frame.to_csv(
                                    outputname,
                                    index=False,
                                    header=False)
                                covariates_frame.to_csv(
                                    output_dir / filename_covs,
                                    index = False,
                                    header = False)

if __name__ == '__main__':
    main_SWME_errorData()