import Visualisation
import PDE
import numpy as np

IC = 'highDamBreak_withVelocity'
pde_type = 'SWME1D'


spatially_adaptive = False
orders = 1, 5, 1
viscosity = 0.1
slip_length = 0.1

monte_carlo = True
distr = "uniform"
order = 1
n_MC = 10
mu = 1.0
sigma = 0.1

stochastic_Galerkin = False
mom_order = 1
SG_order = 1

if pde_type == 'SWME1D':
    _pde = PDE.SWME1D(IC, viscosity, slip_length, hyperbolic=False)
elif pde_type == 'HSWME1D':
    _pde = PDE.SWME1D(IC, viscosity, slip_length, hyperbolic=True)
elif pde_type == 'VegetationSWME1D':
    _pde = PDE.VegetationSWME1D(IC, viscosity, slip_length, False, 1, 1, 1)
elif pde_type == 'SGSWME1D' and stochastic_Galerkin and not spatially_adaptive and not monte_carlo:
    _pde = PDE.SGSWME1D(IC, distr, mu, sigma, slip_length, hyperbolic=False)
elif pde_type == 'SGSWME1D':
    print("pde_type can only be SGSWME1D if stochastic_Galerkin is True and spatially_adaptive and monte_carlo are False")
else:
    print('This pde_type is not implemented yet')


if spatially_adaptive:
    data_array = np.load("Data\data_{0}_orders={1}_nu={2}_lambda={3}_IC={4}.npy".format(pde_type, orders, viscosity, slip_length, IC))
    Visualisation.visualisation(_pde, True, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC, data_array, orders)

elif monte_carlo and pde_type == 'SWME1D' and not stochastic_Galerkin and not spatially_adaptive:
    data_array = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}.npy".format(pde_type, distr, order, n_MC, mu, sigma, slip_length, IC))
    Visualisation.visualisation(_pde, False, True, False, len(data_array[:,0]), data_array[0,0,0], data_array[0,-1,0], IC, data_array, n_MC = n_MC, order = order)

elif stochastic_Galerkin:
    data_array = np.load("Data\data_{0}_{1}_MO={2},SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}.npy".format(pde_type, distr, mom_order, SG_order, mu, sigma, slip_length, IC))
    Visualisation.visualisation(_pde, False, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC, data_array, mom_order = mom_order, SG_order = SG_order)

elif monte_carlo:
    print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

else:
    data_array = np.load("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}.npy".format(pde_type, order, viscosity, slip_length, IC))
    Visualisation.visualisation(_pde, False, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC, data_array, order = order)