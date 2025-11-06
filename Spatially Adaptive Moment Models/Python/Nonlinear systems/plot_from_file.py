import Visualisation
import PDE
import numpy as np
import matplotlib.pyplot as plt

IC_A = 'lowDamBreak_linearVelocity'
IC_B = 'lowDamBreak_linearVelocity'
IC_C = 'lowDamBreak_linearVelocity'
IC_D = 'lowDamBreak_linearVelocity'
pde_type_A = 'SWME1D'
pde_type_B = 'SWME1D'
pde_type_C = 'SGSWME1D'
pde_type_D = 'SGSWME1D'
integrator_A = 'ImplicitEuler'
integrator_B = 'ImplicitEuler'
integrator_C = 'ImplicitEuler'
integrator_D = 'ImplicitEuler'

method_A = 'classical'
method_B = 'classical'
method_C = 'classical'
method_D = 'classical'
start_order_A = 5
start_order_B = 5
start_order_C = 5
start_order_D = 5
max_order_A = 5
max_order_B = 5
max_order_C = 5
max_order_D = 5
viscosity_A = 0.1
viscosity_B = 0.1
viscosity_C = 0.1
viscosity_D = 0.1
slip_length_A = 0.1
slip_length_B = 0.1
slip_length_C = 0.1
slip_length_D = 0.1

linear_source_implicit_A = False
linear_source_implicit_B = False
linear_source_implicit_C = False
linear_source_implicit_D = False
monte_carlo_A = True
monte_carlo_B = True
monte_carlo_C = False
monte_carlo_D = False
distr_A = "uniform"
distr_B = "uniform"
distr_C = "uniform"
distr_D = "uniform"
order_A = 1
order_B = 1
order_C = 1
order_D = 1
n_MC_A = 80
n_MC_B = 90
n_MC_C = 110
n_MC_D = 120
mu_A = 0.1
mu_B = 0.1
mu_C = 0.1
mu_D = 0.1
sigma_A = 0.05
sigma_B = 0.05
sigma_C = 0.05
sigma_D = 0.05

stochastic_Galerkin_A = False
stochastic_Galerkin_B = False
stochastic_Galerkin_C = True
stochastic_Galerkin_D = True
mom_order_A = 1
mom_order_B = 1
mom_order_C = 1
mom_order_D = 1
SG_order_A = 0
SG_order_B = 0
SG_order_C = 1
SG_order_D = 2

title = 'Low Dambreak with Linear Velocity'
t_end = 0.2

label_A = 'Monte Carlo N=110'
label_B = 'Monte Carlo N=120'
label_C = 'Galerkin K=1'
label_D = 'Galerkin K=2'
color_1_A = 'tab:green'
color_2_A = 'greenyellow'
color_3_A = 'forestgreen'
color_1_B = 'tab:red'
color_2_B = 'lightcoral'
color_3_B = 'firebrick'
color_1_C = 'tab:blue'
color_2_C = 'lightskyblue'
color_3_C = 'mediumblue'
color_1_D = 'saddlebrown'
color_2_D = 'olive'
color_3_D = 'slategrey'
fill_A = False
fill_B = False
fill_C = False
fill_D = False
#ls_A = '-'
#ls_B = '--'
#ls_C = '-.'
#ls_D = ':'
ls_A = (0, (3, 9))
ls_B = (3, (3, 9))
ls_C = (6, (3, 9))
ls_D = (9, (3, 9))
lw_A = 2
lw_B = 2
lw_C = 2
lw_D = 2

if pde_type_A == 'SWME1D':
    _pde_A = PDE.SWME1D(IC_A, viscosity_A, slip_length_A, hyperbolic=False, linear_source=linear_source_implicit_A)
elif pde_type_A == 'HSWME1D':
    _pde_A = PDE.SWME1D(IC_A, viscosity_A, slip_length_A, hyperbolic=True, linear_source=linear_source_implicit_A)
elif pde_type_A == 'VegetationSWME1D':
    _pde_A = PDE.VegetationSWME1D(IC_A, viscosity_A, slip_length_A, False, linear_source_implicit_A, 0.008, 1, 264)
elif pde_type_A == 'SGSWME1D' and stochastic_Galerkin_A and not method_A == 'spatially_adaptive' and not method_A == 'micro_macro' and not monte_carlo_A:
    _pde_A = PDE.SGSWME1D(IC_A, distr_A, mu_A, sigma_A, slip_length_A, hyperbolic=False)
elif pde_type_A == 'HSGSWME1D' and stochastic_Galerkin_A and not method_A == 'spatially_adaptive' and not method_A == 'micro_macro' and not monte_carlo_A:
    _pde_A = PDE.SGSWME1D(IC_A, distr_A, mu_A, sigma_A, slip_length_A, hyperbolic=True)
elif pde_type_A == 'SGSWME1D' or pde_type_A == 'HSGSWME1D':
    print("pde_type can only be SGSWME1D if stochastic_Galerkin is True and spatially_adaptive, micro_macro and monte_carlo are False")
else:
    print('This pde_type is not implemented yet')


if pde_type_B == 'SWME1D':
    _pde_B = PDE.SWME1D(IC_B, viscosity_B, slip_length_B, hyperbolic=False, linear_source=linear_source_implicit_B)
elif pde_type_B == 'HSWME1D':
    _pde_B = PDE.SWME1D(IC_B, viscosity_B, slip_length_B, hyperbolic=True, linear_source=linear_source_implicit_B)
elif pde_type_B == 'VegetationSWME1D':
    _pde_B = PDE.VegetationSWME1D(IC_B, viscosity_B, slip_length_B, False, linear_source_implicit_B, 0.008, 1, 264)
elif pde_type_B == 'SGSWME1D' and stochastic_Galerkin_B and not method_B == 'spatially_adaptive' and not method_B == 'micro_macro' and not monte_carlo_B:
    _pde_B = PDE.SGSWME1D(IC_B, distr_B, mu_B, sigma_B, slip_length_B, hyperbolic=False)
elif pde_type_B == 'HSGSWME1D' and stochastic_Galerkin_B and not method_B == 'spatially_adaptive' and not method_B == 'micro_macro' and not monte_carlo_B:
    _pde_B = PDE.SGSWME1D(IC_B, distr_B, mu_B, sigma_B, slip_length_B, hyperbolic=True)
elif pde_type_B == 'SGSWME1D' or pde_type_B == 'HSGSWME1D':
    print("pde_type can only be SGSWME1D if stochastic_Galerkin is True and spatially_adaptive, micro_macro and monte_carlo are False")
else:
    print('This pde_type is not implemented yet')

if pde_type_C == 'SWME1D':
    _pde_C = PDE.SWME1D(IC_C, viscosity_C, slip_length_C, hyperbolic=False, linear_source=linear_source_implicit_C)
elif pde_type_C == 'HSWME1D':
    _pde_C = PDE.SWME1D(IC_C, viscosity_C, slip_length_C, hyperbolic=True, linear_source=linear_source_implicit_C)
elif pde_type_C == 'VegetationSWME1D':
    _pde_C = PDE.VegetationSWME1D(IC_C, viscosity_C, slip_length_C, False, linear_source_implicit_C, 0.008, 1, 264)
elif pde_type_C == 'SGSWME1D' and stochastic_Galerkin_C and not method_C == 'spatially_adaptive' and not method_C == 'micro_macro' and not monte_carlo_C:
    _pde_C = PDE.SGSWME1D(IC_C, distr_C, mu_C, sigma_C, slip_length_C, hyperbolic=False)
elif pde_type_C == 'HSGSWME1D' and stochastic_Galerkin_C and not method_C == 'spatially_adaptive' and not method_C == 'micro_macro' and not monte_carlo_C:
    _pde_C = PDE.SGSWME1D(IC_C, distr_C, mu_C, sigma_C, slip_length_C, hyperbolic=True)
elif pde_type_C == 'SGSWME1D' or pde_type_C == 'HSGSWME1D':
    print("pde_type can only be SGSWME1D if stochastic_Galerkin is True and spatially_adaptive, micro_macro and monte_carlo are False")
else:
    print('This pde_type is not implemented yet')

if pde_type_D == 'SWME1D':
    _pde_D = PDE.SWME1D(IC_D, viscosity_D, slip_length_D, hyperbolic=False, linear_source=linear_source_implicit_D)
elif pde_type_D == 'HSWME1D':
    _pde_D = PDE.SWME1D(IC_D, viscosity_D, slip_length_D, hyperbolic=True, linear_source=linear_source_implicit_D)
elif pde_type_D == 'VegetationSWME1D':
    _pde_D = PDE.VegetationSWME1D(IC_D, viscosity_D, slip_length_D, False, linear_source_implicit_D, 0.008, 1, 264)
elif pde_type_D == 'SGSWME1D' and stochastic_Galerkin_D and not method_D == 'spatially_adaptive' and not method_D == 'micro_macro' and not monte_carlo_D:
    _pde_D = PDE.SGSWME1D(IC_D, distr_D, mu_D, sigma_D, slip_length_D, hyperbolic=False)
elif pde_type_D == 'HSGSWME1D' and stochastic_Galerkin_D and not method_D == 'spatially_adaptive' and not method_D == 'micro_macro' and not monte_carlo_D:
    _pde_D = PDE.SGSWME1D(IC_D, distr_D, mu_D, sigma_D, slip_length_D, hyperbolic=True)
elif pde_type_D == 'SGSWME1D' or pde_type_D == 'HSGSWME1D':
    print("pde_type can only be SGSWME1D if stochastic_Galerkin is True and spatially_adaptive, micro_macro and monte_carlo are False")
else:
    print('This pde_type is not implemented yet')


if method_A == 'spatially_adaptive':
    data_array = np.load("Data\data_{0}_start_order={1}_nu={2}_lambda={3}_IC={4}_T={5}_integrator={6}.npy".format(pde_type_A, start_order_A, viscosity_A, slip_length_A, IC_A, t_end, integrator_A))
    Visualisation.visualisation(_pde_A, method_A, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_A, data_array, title, label_A, fill_A, ls_A, lw_A, max_order = max_order_A, color_1 = color_1_A, color_2 = color_2_A, color_3 = color_3_A)

elif method_A == 'micro_macro':
    print("Plotting is not yet implemented for micro-macro.")

elif monte_carlo_A and pde_type_A == 'SWME1D' and not stochastic_Galerkin_A and not method_A == 'spatially_adaptive' and not method_A == 'micro_macro':
    data_array = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type_A, distr_A, order_A, n_MC_A, mu_A, sigma_A, slip_length_A, IC_A, t_end, integrator_A))
    ax_arr = Visualisation.visualisation(_pde_A, method_A, True, False, len(data_array[0,:,0]), data_array[0,0,0], data_array[0,-1,0], IC_A, data_array, title, label_A, fill_A, ls_A, lw_A, n_MC = n_MC_A, order = order_A, color_1 = color_1_A, color_2 = color_2_A, color_3 = color_3_A)

elif stochastic_Galerkin_A:
    data_array = np.load("Data\data_{0}_{1}_MO={2}_SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type_A, distr_A, mom_order_A, SG_order_A, mu_A, sigma_A, slip_length_A, IC_A, t_end, integrator_A))
    ax_arr = Visualisation.visualisation(_pde_A, method_A, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_A, data_array, title, label_A, fill_A, ls_A, lw_A, mom_order = mom_order_A, SG_order = SG_order_A, color_1 = color_1_A, color_2 = color_2_A, color_3 = color_3_A)

elif monte_carlo_A:
    print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

else:
    data_array = np.load("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}_integrator={6}.npy".format(pde_type_A, order_A, viscosity_A, slip_length_A, IC_A, t_end, integrator_A))
    ax_arr = Visualisation.visualisation(_pde_A, method_A, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_A, data_array, title, label_A, fill_A, ls_A, lw_A, order = order_A, color_2 = color_2_A, color_3 = color_3_A)


if method_B == 'spatially_adaptive':
    data_array = np.load("Data\data_{0}_start_order={1}_nu={2}_lambda={3}_IC={4}_T={5}_integrator={6}.npy".format(pde_type_B, start_order_B, viscosity_B, slip_length_B, IC_B, t_end, integrator_B))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, method_B, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_B, data_array, label_B, fill_B, ls_B, lw_B, max_order = max_order_B, color_1 = color_1_B, color_2 = color_2_B, color_3 = color_3_B)

elif method_B == 'micro_macro':
    print("Plotting is not yet implemented for micro-macro.")

elif monte_carlo_B and pde_type_B == 'SWME1D' and not stochastic_Galerkin_B and not method_B == 'spatially_adaptive' and not method_B == 'micro_macro':
    data_array = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type_B, distr_B, order_B, n_MC_B, mu_B, sigma_B, slip_length_B, IC_B, t_end, integrator_B))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, method_B, True, False, len(data_array[:,0]), data_array[0,0,0], data_array[0,-1,0], IC_B, data_array, label_B, fill_B, ls_B, lw_B, n_MC = n_MC_B, order = order_B, color_1 = color_1_B, color_2 = color_2_B, color_3 = color_3_B)

elif stochastic_Galerkin_B:
    data_array = np.load("Data\data_{0}_{1}_MO={2}_SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type_B, distr_B, mom_order_B, SG_order_B, mu_B, sigma_B, slip_length_B, IC_B, t_end, integrator_B))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, method_B, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_B, data_array, label_B, fill_B, ls_B, lw_B, mom_order = mom_order_B, SG_order = SG_order_B, color_1 = color_1_B, color_2 = color_2_B, color_3 = color_3_B)

elif monte_carlo_B:
    print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

else:
    data_array = np.load("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}_integrator={6}.npy".format(pde_type_B, order_B, viscosity_B, slip_length_B, IC_B, t_end, integrator_B))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, method_B, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_B, data_array, label_B, fill_B, ls_B, lw_B, order = order_B, color_1 = color_1_B, color_2 = color_2_B, color_3 = color_3_B)


if method_C == 'spatially_adaptive':
    data_array = np.load("Data\data_{0}_start_order={1}_nu={2}_lambda={3}_IC={4}_T={5}_integrator={6}.npy".format(pde_type_C, start_order_C, viscosity_C, slip_length_C, IC_C, t_end, integrator_C))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, method_C, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_C, data_array, label_C, fill_C, ls_C, lw_C, max_order = max_order_C, color_1 = color_1_C, color_2 = color_2_C, color_3 = color_3_C)

elif method_C == 'micro_macro':
    print("Plotting is not yet implemented for micro-macro.")

elif monte_carlo_C and pde_type_C == 'SWME1D' and not stochastic_Galerkin_C and not method_C == 'spatially_adaptive' and not method_C == 'micro_macro':
    data_array = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type_C, distr_C, order_C, n_MC_C, mu_C, sigma_C, slip_length_C, IC_C, t_end, integrator_C))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, method_C, True, False, len(data_array[:,0]), data_array[0,0,0], data_array[0,-1,0], IC_C, data_array, label_C, fill_C, ls_C, lw_C, n_MC = n_MC_C, order = order_C, color_1 = color_1_C, color_2 = color_2_C, color_3 = color_3_C)

elif stochastic_Galerkin_C:
    data_array = np.load("Data\data_{0}_{1}_MO={2}_SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type_C, distr_C, mom_order_C, SG_order_C, mu_C, sigma_C, slip_length_C, IC_C, t_end, integrator_C))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, method_C, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_C, data_array, label_C, fill_C, ls_C, lw_C, mom_order = mom_order_C, SG_order = SG_order_C, color_1 = color_1_C, color_2 = color_2_C, color_3 = color_3_C)

elif monte_carlo_C:
    print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

else:
    data_array = np.load("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}_integrator={6}.npy".format(pde_type_C, order_C, viscosity_C, slip_length_C, IC_C, t_end, integrator_C))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, method_C, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_C, data_array, label_C, fill_C, ls_C, lw_C, order = order_C, color_1 = color_1_C, color_2 = color_2_C, color_3 = color_3_C)


if method_D == 'spatially_adaptive':
    data_array = np.load("Data\data_{0}_start_order={1}_nu={2}_lambda={3}_IC={4}_T={5}_integrator={6}.npy".format(pde_type_D, start_order_D, viscosity_D, slip_length_D, IC_D, t_end, integrator_D))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_D, method_D, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_D, data_array, label_D, fill_D, ls_D, lw_D, max_order = max_order_D, color_1 = color_1_D, color_2 = color_2_D, color_3 = color_3_D)

elif method_D == 'micro_macro':
    print("Plotting is not yet implemented for micro-macro.")

elif monte_carlo_D and pde_type_D == 'SWME1D' and not stochastic_Galerkin_D and not method_D == 'spatially_adaptive' and not method_D == 'micro_macro':
    data_array = np.load("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type_D, distr_D, order_D, n_MC_D, mu_D, sigma_D, slip_length_D, IC_D, t_end, integrator_D))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_D, method_D, True, False, len(data_array[:,0]), data_array[0,0,0], data_array[0,-1,0], IC_D, data_array, label_D, fill_D, ls_D, lw_D, n_MC = n_MC_D, order = order_D, color_1 = color_1_D, color_2 = color_2_D, color_3 = color_3_D)

elif stochastic_Galerkin_D:
    data_array = np.load("Data\data_{0}_{1}_MO={2}_SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}_integrator={9}.npy".format(pde_type_D, distr_D, mom_order_D, SG_order_D, mu_D, sigma_D, slip_length_D, IC_D, t_end, integrator_D))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_D, method_D, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_D, data_array, label_D, fill_D, ls_D, lw_D, mom_order = mom_order_D, SG_order = SG_order_D, color_1 = color_1_D, color_2 = color_2_D, color_3 = color_3_D)

elif monte_carlo_D:
    print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

else:
    data_array = np.load("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}_integrator={6}.npy".format(pde_type_D, order_D, viscosity_D, slip_length_D, IC_D, t_end, integrator_D))
    ax_arr = Visualisation.visualisation_add(ax_arr, _pde_D, method_D, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], IC_D, data_array, label_D, fill_D, ls_D, lw_D, order = order_D, color_1 = color_1_D, color_2 = color_2_D, color_3 = color_3_D)

plt.show()