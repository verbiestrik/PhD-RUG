import Visualisation
import PDE
import numpy as np
import matplotlib.pyplot as plt

IC_A = 'lowDamBreak_withVelocity'
IC_B = 'lowDamBreak_withVelocity'
IC_C = 'lowDamBreak_withVelocity'
pde_type_A = 'SGSWME1D'
pde_type_B = 'SWME1D'
pde_type_C = 'SGSWME1D'

spatially_adaptive_A = False
spatially_adaptive_B = False
spatially_adaptive_C = False
orders_A = 1, 5, 1
orders_B = 1, 5, 1
orders_C = 1, 5, 1
viscosity_A = 0.1
viscosity_B = 0.1
viscosity_C = 0.1
slip_length_A = 0.1
slip_length_B = 0.1
slip_length_C = 0.1

linear_source_implicit_A = True
linear_source_implicit_B = True
linear_source_implicit_C = True
monte_carlo_A = False
monte_carlo_B = False
monte_carlo_C = False
distr_A = "uniform"
distr_B = "uniform"
distr_C = "uniform"
order_A = 1
order_B = 1
order_C = 1
n_MC_A = 10
n_MC_B = 10
n_MC_C = 10
mu_A = 0.1
mu_B = 0.1
mu_C = 0.1
sigma_A = 0.01
sigma_B = 0.01
sigma_C = 0.01

stochastic_Galerkin_A = True
stochastic_Galerkin_B = True
stochastic_Galerkin_C = True
mom_order_A = 1
mom_order_B = 1
mom_order_C = 1
SG_order_A = 0
SG_order_B = 1
SG_order_C = 2

title = 'Dam Break of Low Dam'
t_end = 0.2

label_A = 'Test Case'
label_B = 'K=1'
label_C = 'K=2'
color_1_A = 'tab:green'
color_2_A = 'greenyellow'
color_3_A = 'forestgreen'
color_1_B = 'tab:red'
color_2_B = 'lightcoral'
color_3_B = 'firebrick'
color_1_C = 'tab:blue'
color_2_C = 'lightskyblue'
color_3_C = 'mediumblue'

if pde_type_A == 'SWME1D':
    _pde_A = PDE.SWME1D(IC_A, viscosity_A, slip_length_A, hyperbolic=False, linear_source=linear_source_implicit_A)
elif pde_type_A == 'HSWME1D':
    _pde_A = PDE.SWME1D(IC_A, viscosity_A, slip_length_A, hyperbolic=True, linear_source=linear_source_implicit_A)
elif pde_type_A == 'VegetationSWME1D':
    _pde_A = PDE.VegetationSWME1D(IC_A, viscosity_A, slip_length_A, False, linear_source_implicit_A, 0.008, 1, 264)
elif pde_type_A == 'SGSWME1D' and stochastic_Galerkin_A and not spatially_adaptive_A and not monte_carlo_A:
    _pde_A = PDE.SGSWME1D(IC_A, distr_A, mu_A, sigma_A, slip_length_A, hyperbolic=False)
elif pde_type_A == 'HSGSWME1D' and stochastic_Galerkin_A and not spatially_adaptive_A and not monte_carlo_A:
    _pde_A = PDE.SGSWME1D(IC_A, distr_A, mu_A, sigma_A, slip_length_A, hyperbolic=True)
elif pde_type_A == 'SGSWME1D' or pde_type_A == 'HSGSWME1D':
    print("pde_type can only be SGSWME1D if stochastic_Galerkin is True and spatially_adaptive and monte_carlo are False")
else:
    print('This pde_type is not implemented yet')


if pde_type_B == 'SWME1D':
    _pde_B = PDE.SWME1D(IC_B, viscosity_B, slip_length_B, hyperbolic=False, linear_source=linear_source_implicit_B)
elif pde_type_B == 'HSWME1D':
    _pde_B = PDE.SWME1D(IC_B, viscosity_B, slip_length_B, hyperbolic=True, linear_source=linear_source_implicit_B)
elif pde_type_B == 'VegetationSWME1D':
    _pde_B = PDE.VegetationSWME1D(IC_B, viscosity_B, slip_length_B, False, linear_source_implicit_B, 0.008, 1, 264)
elif pde_type_B == 'SGSWME1D' and stochastic_Galerkin_B and not spatially_adaptive_B and not monte_carlo_B:
    _pde_B = PDE.SGSWME1D(IC_B, distr_B, mu_B, sigma_B, slip_length_B, hyperbolic=False)
elif pde_type_B == 'HSGSWME1D' and stochastic_Galerkin_B and not spatially_adaptive_B and not monte_carlo_B:
    _pde_B = PDE.SGSWME1D(IC_B, distr_B, mu_B, sigma_B, slip_length_B, hyperbolic=True)
elif pde_type_B == 'SGSWME1D' or pde_type_B == 'HSGSWME1D':
    print("pde_type can only be SGSWME1D if stochastic_Galerkin is True and spatially_adaptive and monte_carlo are False")
else:
    print('This pde_type is not implemented yet')

if pde_type_C == 'SWME1D':
    _pde_C = PDE.SWME1D(IC_C, viscosity_C, slip_length_C, hyperbolic=False, linear_source=linear_source_implicit_C)
elif pde_type_C == 'HSWME1D':
    _pde_C = PDE.SWME1D(IC_C, viscosity_C, slip_length_C, hyperbolic=True, linear_source=linear_source_implicit_C)
elif pde_type_C == 'VegetationSWME1D':
    _pde_C = PDE.VegetationSWME1D(IC_C, viscosity_C, slip_length_C, False, linear_source_implicit_C, 0.008, 1, 264)
elif pde_type_C == 'SGSWME1D' and stochastic_Galerkin_C and not spatially_adaptive_C and not monte_carlo_C:
    _pde_C = PDE.SGSWME1D(IC_C, distr_C, mu_C, sigma_C, slip_length_C, hyperbolic=False)
elif pde_type_C == 'HSGSWME1D' and stochastic_Galerkin_C and not spatially_adaptive_C and not monte_carlo_C:
    _pde_C = PDE.SGSWME1D(IC_C, distr_C, mu_C, sigma_C, slip_length_C, hyperbolic=True)
elif pde_type_C == 'SGSWME1D' or pde_type_C == 'HSGSWME1D':
    print("pde_type can only be SGSWME1D if stochastic_Galerkin is True and spatially_adaptive and monte_carlo are False")
else:
    print('This pde_type is not implemented yet')


if spatially_adaptive_A:
    data_array = np.load("Data\data_{0}_orders={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_A, orders_A, viscosity_A, slip_length_A, IC_A, t_end))
    Visualisation.visualisation(_pde_A, True, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_A, data_array, title, label_A, orders = orders_A, color_1 = color_1_A, color_2 = color_2_A, color_3 = color_3_A)

elif monte_carlo_A and pde_type_A == 'SWME1D' and not stochastic_Galerkin_A and not spatially_adaptive_A:
    data_array = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_A, distr_A, order_A, n_MC_A, mu_A, sigma_A, slip_length_A, IC_A, t_end))
    ax_arr = Visualisation.visualisation(_pde_A, False, True, False, len(data_array[0,:,0]), data_array[0,0,0], data_array[0,-1,0], IC_A, data_array, title, label_A, n_MC = n_MC_A, order = order_A, color_1 = color_1_A, color_2 = color_2_A, color_3 = color_3_A)

elif stochastic_Galerkin_A:
    data_array = np.load("Data\data_{0}_{1}_MO={2},SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_A, distr_A, mom_order_A, SG_order_A, mu_A, sigma_A, slip_length_A, IC_A, t_end))
    ax_arr = Visualisation.visualisation(_pde_A, False, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_A, data_array, title, label_A, mom_order = mom_order_A, SG_order = SG_order_A, color_1 = color_1_A, color_2 = color_2_A, color_3 = color_3_A)

elif monte_carlo_A:
    print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

else:
    data_array = np.load("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_A, order_A, viscosity_A, slip_length_A, IC_A, t_end))
    ax_arr = Visualisation.visualisation(_pde_A, False, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_A, data_array, title, label_A, order = order_A, color_2 = color_2_A, color_3 = color_3_A)


if spatially_adaptive_B:
    data_array = np.load("Data\data_{0}_orders={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_B, orders_B, viscosity_B, slip_length_B, IC_B, t_end))
    ax_arr = Visualisation.visualisation_add(ax_arr,_pde_B, True, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_B, data_array, label_B, orders = orders_B, color_1 = color_1_B, color_2 = color_2_B, color_3 = color_3_B)

elif monte_carlo_B and pde_type_B == 'SWME1D' and not stochastic_Galerkin_B and not spatially_adaptive_B:
    data_array = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_B, distr_B, order_B, n_MC_B, mu_B, sigma_B, slip_length_B, IC_B, t_end))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, False, True, False, len(data_array[:,0]), data_array[0,0,0], data_array[0,-1,0], IC_B, data_array, label_B, n_MC = n_MC_B, order = order_B, color_1 = color_1_B, color_2 = color_2_B, color_3 = color_3_B)

elif stochastic_Galerkin_B:
    data_array = np.load("Data\data_{0}_{1}_MO={2},SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_B, distr_B, mom_order_B, SG_order_B, mu_B, sigma_B, slip_length_B, IC_B, t_end))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, False, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_B, data_array, label_B, mom_order = mom_order_B, SG_order = SG_order_B, color_1 = color_1_B, color_2 = color_2_B, color_3 = color_3_B)

elif monte_carlo_B:
    print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

else:
    data_array = np.load("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_B, order_B, viscosity_B, slip_length_B, IC_B, t_end))
    ax_arr = Visualisation.visualisation(ax_arr, _pde_B, False, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_B, data_array, label_B, order = order_B, color_1 = color_1_B, color_2 = color_2_B, color_3 = color_3_B)


if spatially_adaptive_C:
    data_array = np.load("Data\data_{0}_orders={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_C, orders_C, viscosity_C, slip_length_C, IC_C, t_end))
    ax_arr = Visualisation.visualisation_add(ax_arr,_pde_C, True, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_C, data_array, label_C, orders = orders_C, color_1 = color_1_C, color_2 = color_2_C, color_3 = color_3_C)

elif monte_carlo_C and pde_type_C == 'SWME1D' and not stochastic_Galerkin_C and not spatially_adaptive_C:
    data_array = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_C, distr_C, order_C, n_MC_C, mu_C, sigma_C, slip_length_C, IC_C, t_end))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, False, True, False, len(data_array[:,0]), data_array[0,0,0], data_array[0,-1,0], IC_B, data_array, label_C, n_MC = n_MC_C, order = order_B, color_1 = color_1_C, color_2 = color_2_C, color_3 = color_3_C)

elif stochastic_Galerkin_B:
    data_array = np.load("Data\data_{0}_{1}_MO={2},SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_C, distr_C, mom_order_C, SG_order_C, mu_C, sigma_C, slip_length_C, IC_C, t_end))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, False, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_C, data_array, label_C, mom_order = mom_order_C, SG_order = SG_order_C, color_1 = color_1_C, color_2 = color_2_C, color_3 = color_3_C)

elif monte_carlo_B:
    print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

else:
    data_array = np.load("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_C, order_C, viscosity_C, slip_length_C, IC_C, t_end))
    ax_arr = Visualisation.visualisation(ax_arr, _pde_C, False, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_C, data_array, label_C, order = order_C, color_1 = color_1_C, color_2 = color_2_C, color_3 = color_3_C)

plt.show()