import pandas as pd
import numpy as np
from pathlib import Path

viscosities = [0.05,0.5,5.0]
slip_lengths = [0.05,0.5,5.0]
orders = [0,1,2,3,4,5]
initialHeight_names = ["hDamLarge","hDamMid", "hDamSmall", "hSmoothSteep", "hSmoothMid","hSmoothSmall"]
initialVelocity_names = ["uConst", "uVarXnotvarZ", "uVarZnotvarX","uVarXvarZ"]
initialVelocityMagnitude_names = ["uFast","uMid","uSlow"]
num_scheme_name = "Roe"

mesh_resolution = 400
number_of_moments = 6

number_of_rows = len(viscosities)*len(slip_lengths)*len(orders)*len(initialHeight_names)*len(initialVelocity_names)*len(initialVelocityMagnitude_names)

error_output = np.zeros((number_of_rows))
design_matrix = np.zeros((number_of_rows,19))

base_dir = Path('Data-processing')
dir_output = base_dir / 'Output' / 'HonoursProject-Cyril'
base_dir = base_dir / 'Output' / 'HonoursProject-Cyril' / 'error_trainingData'

row = 0
for i in range(len(viscosities)):
    for j in range(len(slip_lengths)):
        for k in range(len(orders)):
            for m in range(len(initialHeight_names)):
                for n in range(len(initialVelocity_names)):
                    for o in range(len(initialVelocityMagnitude_names)):

                        foldername = initialHeight_names[m]+'_'+initialVelocity_names[n]+'_'+initialVelocityMagnitude_names[o]

                        input_dir = base_dir / foldername

                        filename = 'lambda'+str(slip_lengths[j])+'_viscosity'+str(viscosities[i])+'_order'+str(orders[k])+'_FVM'+num_scheme_name+'.csv'
                        filename_ref = 'lambda'+str(slip_lengths[j])+'_viscosity'+str(viscosities[i])+'_order'+str(6)+'_FVM'+num_scheme_name+'.csv'
                        filename_covs = 'covs_'+filename
                        filename_covs_ref = 'covs_'+filename_ref

                        inputname = input_dir / filename
                        inputname_ref = input_dir / filename_ref
                        inputname_covs = input_dir / filename_covs
                        inputname_covs_ref = input_dir / filename_covs_ref

                        data_matrix = np.loadtxt(inputname, delimiter=',')
                        data_matrix_ref = np.loadtxt(inputname_ref, delimiter=',')
                        covs_vector = np.loadtxt(inputname_covs,delimiter=',')
                        covs_vector_ref = np.loadtxt(inputname_covs_ref,delimiter=',')
                        full_data_matrix = np.zeros((mesh_resolution,number_of_moments+3))
                        full_data_matrix[:,:data_matrix.shape[1]] = data_matrix

                        err_h = np.linalg.norm(full_data_matrix[:,1]-data_matrix_ref[:,1])/np.linalg.norm(data_matrix_ref[:,1])
                        err_u = np.linalg.norm(full_data_matrix[:,2]-data_matrix_ref[:,2])/np.linalg.norm(data_matrix_ref[:,2])
                        # err_alpha1 = np.linalg.norm(full_data_matrix[:,3]-data_matrix_ref[:,3])/np.linalg.norm(data_matrix_ref[:,3])
                        # err_alpha2 = np.linalg.norm(full_data_matrix[:,4]-data_matrix_ref[:,4])/np.linalg.norm(data_matrix_ref[:,4])
                        # err_alpha3 = np.linalg.norm(full_data_matrix[:,5]-data_matrix_ref[:,5])/np.linalg.norm(data_matrix_ref[:,5])
                        # err_alpha4 = np.linalg.norm(full_data_matrix[:,6]-data_matrix_ref[:,6])/np.linalg.norm(data_matrix_ref[:,6])
                        # err_alpha5 = np.linalg.norm(full_data_matrix[:,7]-data_matrix_ref[:,7])/np.linalg.norm(data_matrix_ref[:,7])
                        # err_alpha6 = np.linalg.norm(full_data_matrix[:,8]-data_matrix_ref[:,8])/np.linalg.norm(data_matrix_ref[:,8])

                        error_output[row] = np.sqrt(err_h**2+err_u**2)
                        design_matrix[row,:-1] = covs_vector_ref[:-1]
                        design_matrix[row,-1] = orders[k]
                        
                        row += 1

column_names = ["h", "h*um", "h*alpha1","h*alpha2","h*alpha3","h*alpha4","h*alpha5","h*alpha6",\
                "d_h", "d_h*um", "d_h*alpha1","d_h*alpha2","d_h*alpha3","d_h*alpha4","d_h*alpha5",\
                    "d_h*alpha6","viscosity","slip length","order"]

column_names_with_output = ["output","h", "h*um", "h*alpha1","h*alpha2","h*alpha3","h*alpha4","h*alpha5","h*alpha6",\
                "d_h", "d_h*um", "d_h*alpha1","d_h*alpha2","d_h*alpha3","d_h*alpha4","d_h*alpha5",\
                    "d_h*alpha6","viscosity","slip length","order"]

error_and_designMatrix_combined = np.zeros((number_of_rows,20))
error_and_designMatrix_combined[:,0] = error_output
error_and_designMatrix_combined[:,1:] = design_matrix

data_frame_errorOutput = pd.DataFrame(error_output)
data_frame_designMatrix = pd.DataFrame(design_matrix,columns=column_names)
data_frame_error_and_designMatrix_combined = pd.DataFrame(error_and_designMatrix_combined,columns=column_names_with_output)

filename_errorOutput = 'error_values.csv'
filename_designMatrix = 'design_matrix.csv'
filename_error_and_designMatrix_combined = 'error_and_designMatrix_combined.csv'

outputname_errorOutput = dir_output / filename_errorOutput
outputname_designMatrix = dir_output / filename_designMatrix
outputname_error_and_designMatrix_combined = dir_output / filename_error_and_designMatrix_combined

data_frame_errorOutput.to_csv(
    outputname_errorOutput,
    index=False,
    header=False)
data_frame_designMatrix.to_csv(
    outputname_designMatrix,
    index = False,
    header = True)
data_frame_error_and_designMatrix_combined.to_csv(
    outputname_error_and_designMatrix_combined,
    index = False,
    header = True)