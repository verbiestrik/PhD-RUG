import numpy as np

N_A = 1500
N_B = 50

pde_type    = 'SWLME1D'
distr       = 'uniform'
order       = 1
mu          = 0.1
sigma       = 0.05
slip_length = 0.1
IC          = 'smoothWave_linearVelocity'
T           = 2.0

data_array_C = np.load("Data/data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type, distr, order, N_A+N_B, mu, sigma, slip_length, IC, T))

data_array_A = data_array_C[0:N_A,:]
data_array_B = data_array_C[N_A:N_A+N_B,:]
np.save("Data/data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type, distr, order, N_A, mu, sigma, slip_length, IC, T), data_array_A)
np.save("Data/data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type, distr, order, N_B, mu, sigma, slip_length, IC, T), data_array_B)