import numpy as np

N_A = 20
N_B = 60

pde_type    = 'SWME1D'
distr       = 'uniform'
order       = 1
mu          = 0.1
sigma       = 0.05
slip_length = 0.1
IC          = 'lowDamBreak_linearVelocity'
T           = 0.2
integrator  = 'ImplicitEuler'

data_array_A = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type, distr, order, N_A, mu, sigma, slip_length, IC, T, integrator))
data_array_B = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type, distr, order, N_B, mu, sigma, slip_length, IC, T, integrator))

data_array_C = np.concatenate((data_array_A, data_array_B), axis=0)
np.save("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type, distr, order, N_A+N_B, mu, sigma, slip_length, IC, T, integrator), data_array_C)