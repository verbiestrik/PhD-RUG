import numpy as np
from pathlib import Path

mesh_resolution = 400
number_of_moments = 6

viscosities = [0.05]
slip_lengths = [0.05]
orders = [0]
initialHeight_names = ["hDamLarge"]
initialVelocity_names = ["uConst"]
initialVelocityMagnitude_names = ["uSlow"]
num_scheme_names = ["Roe"]

foldername = initialHeight_names[0]+'_'+initialVelocity_names[0]+'_'+initialVelocityMagnitude_names[0]

base_dir = Path('Data-processing')
output_dir = base_dir / 'Output' / 'HonoursProject-Cyril'
base_dir = base_dir / 'Output' / 'HonoursProject-Cyril' / 'error_trainingData'
input_dir = base_dir / foldername

# # Create folders if they don't exist
# output_dir.mkdir(parents=True, exist_ok=True)

filename = 'lambda'+str(slip_lengths[0])+'_viscosity'+str(viscosities[0])+'_order'+str(orders[0])+'_FVM'+num_scheme_names[0]+'.csv'
filename_ref = 'lambda'+str(slip_lengths[0])+'_viscosity'+str(viscosities[0])+'_order'+str(6)+'_FVM'+num_scheme_names[0]+'.csv'
filename_covs = 'covs_'+filename

inputname = input_dir / filename
inputname_ref = input_dir / filename_ref
inputname_covs = input_dir / filename_covs

data_matrix = np.loadtxt(inputname, delimiter=',')
data_matrix_ref = np.loadtxt(inputname_ref, delimiter=',')
covs_vector = np.loadtxt(inputname_covs,delimiter=',')
full_data_matrix = np.zeros((mesh_resolution,number_of_moments+3))
full_data_matrix[:,:data_matrix.shape[1]] = data_matrix

err_h = np.linalg.norm(full_data_matrix[:,1]-data_matrix_ref[:,1])/np.linalg.norm(data_matrix_ref[:,1])
err_u = np.linalg.norm(full_data_matrix[:,2]-data_matrix_ref[:,2])/np.linalg.norm(data_matrix_ref[:,2])
err_alpha1 = np.linalg.norm(full_data_matrix[:,3]-data_matrix_ref[:,3])/np.linalg.norm(data_matrix_ref[:,3])
err_alpha2 = np.linalg.norm(full_data_matrix[:,4]-data_matrix_ref[:,4])/np.linalg.norm(data_matrix_ref[:,4])
err_alpha3 = np.linalg.norm(full_data_matrix[:,5]-data_matrix_ref[:,5])/np.linalg.norm(data_matrix_ref[:,5])
err_alpha4 = np.linalg.norm(full_data_matrix[:,6]-data_matrix_ref[:,6])/np.linalg.norm(data_matrix_ref[:,6])
err_alpha5 = np.linalg.norm(full_data_matrix[:,7]-data_matrix_ref[:,7])/np.linalg.norm(data_matrix_ref[:,7])
err_alpha6 = np.linalg.norm(full_data_matrix[:,8]-data_matrix_ref[:,8])/np.linalg.norm(data_matrix_ref[:,8])

error_output = np.zeros((3888))
design_matrix = np.zeros((1,19))

design_matrix[0,:] = covs_vector
print(design_matrix[0,:])

print(err_u)