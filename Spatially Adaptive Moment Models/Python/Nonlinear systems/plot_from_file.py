import Visualisation
import PDE
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
no_lines = 3    # Must be between 1-4
crop_l   = 0.2  # Fraction to crop left side of domain by
crop_r   = 0.8  # Fraction to crop right side of domain by
show_only_one = 0 # Set 0 if you want to display all axes
show_diff     = False
display_title = False

IC_A = 'lowDamBreak_linearVelocity'
IC_B = 'lowDamBreak_linearVelocity'
IC_C = 'lowDamBreak_linearVelocity'
IC_D = 'lowDamBreak_linearVelocity'
pde_type_A = 'SWLME1D'
pde_type_B = 'SGSWLME1D'
pde_type_C = 'SGSWLME1D'
pde_type_D = 'SGSWLME1D'

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
monte_carlo_B = False
monte_carlo_C = False
monte_carlo_D = False
distr_A = "uniform"
distr_B = "uniform"
distr_C = "uniform"
distr_D = "uniform"
order_A = 2
order_B = 2
order_C = 2
order_D = 2
n_MC_A = 50
n_MC_B = 100
n_MC_C = 150
n_MC_D = 200
mu_A = 0.1
mu_B = 0.1
mu_C = 0.1
mu_D = 0.1
sigma_A = 0.05
sigma_B = 0.05
sigma_C = 0.05
sigma_D = 0.05

stochastic_Galerkin_A = False
stochastic_Galerkin_B = True
stochastic_Galerkin_C = True
stochastic_Galerkin_D = True
mom_order_A = 2
mom_order_B = 2
mom_order_C = 2
mom_order_D = 2
SG_order_A = 0
SG_order_B = 0
SG_order_C = 1
SG_order_D = 2

title = 'Low Dam Break with Linear Velocity, SWLME N=1 MC'
t_end = 0.2

label_A = 'MC S=200'
label_B = 'SG K=0'
label_C = 'SG K=1'
label_D = 'SG K=2'
color_1_A = 'black'
color_2_A = 'gray'
color_3_A = 'gray'
color_1_B = 'firebrick'
color_2_B = 'indianred'
color_3_B = 'indianred'
color_1_C = 'royalblue'
color_2_C = 'cornflowerblue'
color_3_C = 'cornflowerblue'
color_1_D = 'seagreen'
color_2_D = 'mediumseagreen'
color_3_D = 'mediumseagreen'

fill_A = False
fill_B = False
fill_C = False
fill_D = False
ls_A = ":"
ls_B = "-."
ls_C = "--"
ls_D = "-"
lw_A = 3
lw_B = 1.5
lw_C = 1.5
lw_D = 1.5


if no_lines > 0:
    if pde_type_A == 'SWME1D':
        _pde_A = PDE.SWME1D(IC_A, viscosity_A, slip_length_A, hyperbolic=False, linearised=False, linear_source=linear_source_implicit_A)
    elif pde_type_A == 'HSWME1D':
        _pde_A = PDE.SWME1D(IC_A, viscosity_A, slip_length_A, hyperbolic=True, linearised=False, linear_source=linear_source_implicit_A)
    elif pde_type_A == 'SWLME1D':
        _pde_A = PDE.SWME1D(IC_A, viscosity_A, slip_length_A, hyperbolic=False, linearised=True, linear_source=linear_source_implicit_A)
    elif pde_type_A == 'VegetationSWME1D':
        _pde_A = PDE.VegetationSWME1D(IC_A, viscosity_A, slip_length_A, False, False, linear_source_implicit_A, 0.008, 1, 264)
    elif pde_type_A == 'SGSWLME1D' and stochastic_Galerkin_A and not method_A == 'spatially_adaptive' and not method_A == 'micro_macro' and not monte_carlo_A:
        _pde_A = PDE.SGSWLME1D(IC_A, distr_A, mu_A, sigma_A, slip_length_A, hyperbolic=False)
    elif pde_type_A == 'HSGSWLME1D' and stochastic_Galerkin_A and not method_A == 'spatially_adaptive' and not method_A == 'micro_macro' and not monte_carlo_A:
        _pde_A = PDE.SGSWLME1D(IC_A, distr_A, mu_A, sigma_A, slip_length_A, hyperbolic=True)
    elif pde_type_A == 'SGSWLME1D' or pde_type_A == 'HSGSWLME1D':
        print("pde_type can only be SGSWLME1D if stochastic_Galerkin is True and spatially_adaptive, micro_macro and monte_carlo are False")
    else:
        print('This pde_type is not implemented yet')

if no_lines > 1:
    if pde_type_B == 'SWME1D':
        _pde_B = PDE.SWME1D(IC_B, viscosity_B, slip_length_B, hyperbolic=False, linearised=False, linear_source=linear_source_implicit_B)
    elif pde_type_B == 'HSWME1D':
        _pde_B = PDE.SWME1D(IC_B, viscosity_B, slip_length_B, hyperbolic=True, linearised=False, linear_source=linear_source_implicit_B)
    elif pde_type_B == 'SWLME1D':
        _pde_B = PDE.SWME1D(IC_B, viscosity_B, slip_length_B, hyperbolic=False, linearised=True, linear_source=linear_source_implicit_B)
    elif pde_type_B == 'VegetationSWME1D':
        _pde_B = PDE.VegetationSWME1D(IC_B, viscosity_B, slip_length_B, False, False, linear_source_implicit_B, 0.008, 1, 264)
    elif pde_type_B == 'SGSWLME1D' and stochastic_Galerkin_B and not method_B == 'spatially_adaptive' and not method_B == 'micro_macro' and not monte_carlo_B:
        _pde_B = PDE.SGSWLME1D(IC_B, distr_B, mu_B, sigma_B, slip_length_B, hyperbolic=False)
    elif pde_type_B == 'HSGSWLME1D' and stochastic_Galerkin_B and not method_B == 'spatially_adaptive' and not method_B == 'micro_macro' and not monte_carlo_B:
        _pde_B = PDE.SGSWLME1D(IC_B, distr_B, mu_B, sigma_B, slip_length_B, hyperbolic=True)
    elif pde_type_B == 'SGSWLME1D' or pde_type_B == 'HSGSWLME1D':
        print("pde_type can only be SGSWLME1D if stochastic_Galerkin is True and spatially_adaptive, micro_macro and monte_carlo are False")
    else:
        print('This pde_type is not implemented yet')

if no_lines > 2:
    if pde_type_C == 'SWME1D':
        _pde_C = PDE.SWME1D(IC_C, viscosity_C, slip_length_C, hyperbolic=False, linearised=False, linear_source=linear_source_implicit_C)
    elif pde_type_C == 'HSWME1D':
        _pde_C = PDE.SWME1D(IC_C, viscosity_C, slip_length_C, hyperbolic=True, linearised=False, linear_source=linear_source_implicit_C)
    elif pde_type_C == 'SWLME1D':
        _pde_C = PDE.SWME1D(IC_C, viscosity_C, slip_length_C, hyperbolic=False, linearised=True, linear_source=linear_source_implicit_C)
    elif pde_type_C == 'VegetationSWME1D':
        _pde_C = PDE.VegetationSWME1D(IC_C, viscosity_C, slip_length_C, False, False, linear_source_implicit_C, 0.008, 1, 264)
    elif pde_type_C == 'SGSWLME1D' and stochastic_Galerkin_C and not method_C == 'spatially_adaptive' and not method_C == 'micro_macro' and not monte_carlo_C:
        _pde_C = PDE.SGSWLME1D(IC_C, distr_C, mu_C, sigma_C, slip_length_C, hyperbolic=False)
    elif pde_type_C == 'HSGSWLME1D' and stochastic_Galerkin_C and not method_C == 'spatially_adaptive' and not method_C == 'micro_macro' and not monte_carlo_C:
        _pde_C = PDE.SGSWLME1D(IC_C, distr_C, mu_C, sigma_C, slip_length_C, hyperbolic=True)
    elif pde_type_C == 'SGSWLME1D' or pde_type_C == 'HSGSWLME1D':
        print("pde_type can only be SGSWLME1D if stochastic_Galerkin is True and spatially_adaptive, micro_macro and monte_carlo are False")
    else:
        print('This pde_type is not implemented yet')

if no_lines > 3:
    if pde_type_D == 'SWME1D':
        _pde_D = PDE.SWME1D(IC_D, viscosity_D, slip_length_D, hyperbolic=False, linearised=False, linear_source=linear_source_implicit_D)
    elif pde_type_D == 'HSWME1D':
        _pde_D = PDE.SWME1D(IC_D, viscosity_D, slip_length_D, hyperbolic=True, linearised=False, linear_source=linear_source_implicit_D)
    elif pde_type_D == 'SWLME1D':
        _pde_D = PDE.SWME1D(IC_D, viscosity_D, slip_length_D, hyperbolic=False, linearised=True, linear_source=linear_source_implicit_D)
    elif pde_type_D == 'VegetationSWME1D':
        _pde_D = PDE.VegetationSWME1D(IC_D, viscosity_D, slip_length_D, False, False, linear_source_implicit_D, 0.008, 1, 264)
    elif pde_type_D == 'SGSWLME1D' and stochastic_Galerkin_D and not method_D == 'spatially_adaptive' and not method_D == 'micro_macro' and not monte_carlo_D:
        _pde_D = PDE.SGSWLME1D(IC_D, distr_D, mu_D, sigma_D, slip_length_D, hyperbolic=False)
    elif pde_type_D == 'HSGSWLME1D' and stochastic_Galerkin_D and not method_D == 'spatially_adaptive' and not method_D == 'micro_macro' and not monte_carlo_D:
        _pde_D = PDE.SGSWLME1D(IC_D, distr_D, mu_D, sigma_D, slip_length_D, hyperbolic=True)
    elif pde_type_D == 'SGSWLME1D' or pde_type_D == 'HSGSWLME1D':
        print("pde_type can only be SGSWLME1D if stochastic_Galerkin is True and spatially_adaptive, micro_macro and monte_carlo are False")
    else:
        print('This pde_type is not implemented yet')


if no_lines > 0:
    if method_A == 'spatially_adaptive':
        data_array = np.load("Data/data_{0}_start_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_A, start_order_A, viscosity_A, slip_length_A, IC_A, t_end))
        Visualisation.visualisation(_pde_A, method_A, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, title, label_A, fill_A, ls_A, lw_A, max_order=max_order_A, color_1=color_1_A, color_2=color_2_A, color_3=color_3_A, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif method_A == 'micro_macro':
        print("Plotting is not yet implemented for micro-macro.")

    elif monte_carlo_A and pde_type_A == 'SWLME1D' and not stochastic_Galerkin_A and not method_A == 'spatially_adaptive' and not method_A == 'micro_macro':
        data_array = np.load("Data/data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_A, distr_A, order_A, n_MC_A, mu_A, sigma_A, slip_length_A, IC_A, t_end))
        ax_arr = Visualisation.visualisation(_pde_A, method_A, True, False, len(data_array[0,:,0]), data_array[0,0,0], data_array[0,-1,0], data_array, title, label_A, fill_A, ls_A, lw_A, n_MC=n_MC_A, order=order_A, color_1=color_1_A, color_2=color_2_A, color_3=color_3_A, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif stochastic_Galerkin_A:
        data_array = np.load("Data/data_{0}_{1}_MO={2}_SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_A, distr_A, mom_order_A, SG_order_A, mu_A, sigma_A, slip_length_A, IC_A, t_end))
        ax_arr = Visualisation.visualisation(_pde_A, method_A, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, title, label_A, fill_A, ls_A, lw_A, mom_order=mom_order_A, SG_order=SG_order_A, color_1=color_1_A, color_2=color_2_A, color_3=color_3_A, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif monte_carlo_A:
        print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

    else:
        data_array = np.load("Data/data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_A, order_A, viscosity_A, slip_length_A, IC_A, t_end))
        ax_arr = Visualisation.visualisation(_pde_A, method_A, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, title, label_A, fill_A, ls_A, lw_A, order=order_A, color_2=color_2_A, color_3=color_3_A, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

if no_lines > 1:
    if method_B == 'spatially_adaptive':
        data_array = np.load("Data/data_{0}_start_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_B, start_order_B, viscosity_B, slip_length_B, IC_B, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, method_B, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_B, fill_B, ls_B, lw_B, max_order=max_order_B, color_1=color_1_B, color_2=color_2_B, color_3=color_3_B, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif method_B == 'micro_macro':
        print("Plotting is not yet implemented for micro-macro.")

    elif monte_carlo_B and pde_type_B == 'SWLME1D' and not stochastic_Galerkin_B and not method_B == 'spatially_adaptive' and not method_B == 'micro_macro':
        data_array = np.load("Data/data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_B, distr_B, order_B, n_MC_B, mu_B, sigma_B, slip_length_B, IC_B, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, method_B, True, False, len(data_array[0,:,0]), data_array[0,0,0], data_array[0,-1,0], data_array, label_B, fill_B, ls_B, lw_B, n_MC=n_MC_B, order=order_B, color_1=color_1_B, color_2=color_2_B, color_3=color_3_B, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif stochastic_Galerkin_B:
        data_array = np.load("Data/data_{0}_{1}_MO={2}_SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_B, distr_B, mom_order_B, SG_order_B, mu_B, sigma_B, slip_length_B, IC_B, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, method_B, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_B, fill_B, ls_B, lw_B, mom_order=mom_order_B, SG_order=SG_order_B, color_1=color_1_B, color_2=color_2_B, color_3=color_3_B, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif monte_carlo_B:
        print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

    else:
        data_array = np.load("Data/data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_B, order_B, viscosity_B, slip_length_B, IC_B, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_B, method_B, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_B, fill_B, ls_B, lw_B, order=order_B, color_1=color_1_B, color_2=color_2_B, color_3=color_3_B, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

if no_lines > 2:
    if method_C == 'spatially_adaptive':
        data_array = np.load("Data/data_{0}_start_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_C, start_order_C, viscosity_C, slip_length_C, IC_C, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, method_C, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_C, fill_C, ls_C, lw_C, max_order=max_order_C, color_1=color_1_C, color_2=color_2_C, color_3=color_3_C, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif method_C == 'micro_macro':
        print("Plotting is not yet implemented for micro-macro.")

    elif monte_carlo_C and pde_type_C == 'SWLME1D' and not stochastic_Galerkin_C and not method_C == 'spatially_adaptive' and not method_C == 'micro_macro':
        data_array = np.load("Data/data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_C, distr_C, order_C, n_MC_C, mu_C, sigma_C, slip_length_C, IC_C, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, method_C, True, False, len(data_array[0,:,0]), data_array[0,0,0], data_array[0,-1,0], data_array, label_C, fill_C, ls_C, lw_C, n_MC=n_MC_C, order=order_C, color_1=color_1_C, color_2=color_2_C, color_3=color_3_C, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif stochastic_Galerkin_C:
        data_array = np.load("Data/data_{0}_{1}_MO={2}_SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_C, distr_C, mom_order_C, SG_order_C, mu_C, sigma_C, slip_length_C, IC_C, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, method_C, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_C, fill_C, ls_C, lw_C, mom_order=mom_order_C, SG_order=SG_order_C, color_1=color_1_C, color_2=color_2_C, color_3=color_3_C, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif monte_carlo_C:
        print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

    else:
        data_array = np.load("Data/data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_C, order_C, viscosity_C, slip_length_C, IC_C, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_C, method_C, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_C, fill_C, ls_C, lw_C, order=order_C, color_1=color_1_C, color_2=color_2_C, color_3=color_3_C, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

if no_lines > 3:
    if method_D == 'spatially_adaptive':
        data_array = np.load("Data/data_{0}_start_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_D, start_order_D, viscosity_D, slip_length_D, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_D, method_D, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_D, fill_D, ls_D, lw_D, max_order=max_order_D, color_1=color_1_D, color_2=color_2_D, color_3=color_3_D, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif method_D == 'micro_macro':
        print("Plotting is not yet implemented for micro-macro.")

    elif monte_carlo_D and pde_type_D == 'SWLME1D' and not stochastic_Galerkin_D and not method_D == 'spatially_adaptive' and not method_D == 'micro_macro':
        data_array = np.load("Data/data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_D, distr_D, order_D, n_MC_D, mu_D, sigma_D, slip_length_D, IC_D, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_D, method_D, True, False, len(data_array[0,:,0]), data_array[0,0,0], data_array[0,-1,0], data_array, label_D, fill_D, ls_D, lw_D, n_MC=n_MC_D, order=order_D, color_1=color_1_D, color_2=color_2_D, color_3=color_3_D, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif stochastic_Galerkin_D:
        data_array = np.load("Data/data_{0}_{1}_MO={2}_SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_type_D, distr_D, mom_order_D, SG_order_D, mu_D, sigma_D, slip_length_D, IC_D, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_D, method_D, False, True, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_D, fill_D, ls_D, lw_D, mom_order=mom_order_D, SG_order=SG_order_D, color_1=color_1_D, color_2=color_2_D, color_3=color_3_D, display_title=display_title, crop_l=crop_l, crop_r=crop_r)

    elif monte_carlo_D:
        print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")

    else:
        data_array = np.load("Data/data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_type_D, order_D, viscosity_D, slip_length_D, IC_D, t_end))
        ax_arr = Visualisation.visualisation_add(ax_arr, _pde_D, method_D, False, False, len(data_array[:,0]), data_array[0,0], data_array[-1,0], data_array, label_D, fill_D, ls_D, lw_D, order=order_D, color_1=color_1_D, color_2=color_2_D, color_3=color_3_D, display_title=display_title, crop_l=crop_l, crop_r=crop_r)


if show_only_one == 1:
    ax_arr[1].remove()
    if mom_order_A > 0 or order_A > 0:
        ax_arr[2].remove()
    if mom_order_A > 1 or order_A > 1:
        ax_arr[3].remove()
    ax_arr[0].set_position([0.12,0.12,0.87,0.87])

elif show_only_one == 2:
    ax_arr[0].remove()
    if mom_order_A > 0 or order_A > 0:
        ax_arr[2].remove()
    if mom_order_A > 1 or order_A > 1:
        ax_arr[3].remove()
    ax_arr[1].set_position([0.12,0.12,0.87,0.87])

elif show_only_one == 3:
    ax_arr[0].remove()
    ax_arr[1].remove()
    if mom_order_A > 1 or order_A > 1:
        ax_arr[3].remove()
    ax_arr[2].set_position([0.12,0.12,0.87,0.87])

elif show_only_one == 4:
    ax_arr[0].remove()
    ax_arr[1].remove()
    ax_arr[2].remove()
    ax_arr[3].set_position([0.12,0.12,0.87,0.87])


if show_diff == True:
    h_1_exp = np.array(ax_arr[0].lines[0].get_ydata())
    h_1_var = np.array(ax_arr[0].lines[1].get_ydata()) - h_1_exp
    h_2_exp = np.array(ax_arr[0].lines[3].get_ydata())
    h_2_var = np.array(ax_arr[0].lines[4].get_ydata()) - h_2_exp
    h_x = ax_arr[0].lines[0].get_xdata()
    h_2_exp_diff = (h_2_exp - h_1_exp)/h_1_exp
    h_2_var_diff = (h_2_var - h_1_var)/h_1_exp

    um_1_exp = np.array(ax_arr[1].lines[0].get_ydata())
    um_1_var = np.array(ax_arr[1].lines[1].get_ydata()) - um_1_exp
    um_2_exp = np.array(ax_arr[1].lines[3].get_ydata())
    um_2_var = np.array(ax_arr[1].lines[4].get_ydata()) - um_2_exp
    um_x = ax_arr[1].lines[0].get_xdata()
    um_2_exp_diff = (um_2_exp - um_1_exp)/um_1_exp
    um_2_var_diff = (um_2_var - um_1_var)/um_1_var

    if order_A > 0 or mom_order_A > 0:
        alpha1_1_exp = np.array(ax_arr[2].lines[0].get_ydata())
        alpha1_1_var = np.array(ax_arr[2].lines[1].get_ydata()) - alpha1_1_exp
        alpha1_2_exp = np.array(ax_arr[2].lines[3].get_ydata())
        alpha1_2_var = np.array(ax_arr[2].lines[4].get_ydata()) - alpha1_2_exp
        alpha1_x = ax_arr[2].lines[0].get_xdata()
        alpha1_2_exp_diff = (alpha1_2_exp - alpha1_1_exp)/alpha1_1_exp
        alpha1_2_var_diff = (alpha1_2_var - alpha1_1_var)/alpha1_1_var
    
    if order_A > 1 or mom_order_A > 1:
        alpha2_1_exp = np.array(ax_arr[3].lines[0].get_ydata())
        alpha2_1_var = np.array(ax_arr[3].lines[1].get_ydata()) - alpha2_1_exp
        alpha2_2_exp = np.array(ax_arr[3].lines[3].get_ydata())
        alpha2_2_var = np.array(ax_arr[3].lines[4].get_ydata()) - alpha2_2_exp
        alpha2_x = ax_arr[3].lines[0].get_xdata()
        alpha2_2_exp_diff = (alpha2_2_exp - alpha2_1_exp)/alpha2_1_exp
        alpha2_2_var_diff = (alpha2_2_var - alpha2_1_var)/alpha2_1_var

    if no_lines > 2:
        h_3_exp = np.array(ax_arr[0].lines[6].get_ydata())
        h_3_var = np.array(ax_arr[0].lines[7].get_ydata()) - h_3_exp
        h_3_exp_diff = (h_3_exp - h_1_exp)/h_1_exp
        h_3_var_diff = (h_3_var - h_1_var)/h_1_exp

        um_3_exp = np.array(ax_arr[1].lines[6].get_ydata())
        um_3_var = np.array(ax_arr[1].lines[7].get_ydata()) - um_3_exp
        um_3_exp_diff = (um_3_exp - um_1_exp)/um_1_exp
        um_3_var_diff = (um_3_var - um_1_var)/um_1_var

        if order_A > 0 or mom_order_A > 0:
            alpha1_3_exp = np.array(ax_arr[2].lines[6].get_ydata())
            alpha1_3_var = np.array(ax_arr[2].lines[7].get_ydata()) - alpha1_3_exp
            alpha1_3_exp_diff = (alpha1_3_exp - alpha1_1_exp)/alpha1_1_exp
            alpha1_3_var_diff = (alpha1_3_var - alpha1_1_var)/alpha1_1_var
        
        if order_A > 1 or mom_order_A > 1:
            alpha2_3_exp = np.array(ax_arr[3].lines[6].get_ydata())
            alpha2_3_var = np.array(ax_arr[3].lines[7].get_ydata()) - alpha2_3_exp
            alpha2_3_exp_diff = (alpha2_3_exp - alpha2_1_exp)/alpha2_1_exp
            alpha2_3_var_diff = (alpha2_3_var - alpha2_1_var)/alpha2_1_var
    
    if no_lines > 3:
        h_4_exp = np.array(ax_arr[0].lines[9].get_ydata())
        h_4_var = np.array(ax_arr[0].lines[10].get_ydata()) - h_4_exp
        h_4_exp_diff = (h_4_exp - h_1_exp)/h_1_exp
        h_4_var_diff = (h_4_var - h_1_var)/h_1_exp

        um_4_exp = np.array(ax_arr[1].lines[9].get_ydata())
        um_4_var = np.array(ax_arr[1].lines[10].get_ydata()) - um_4_exp
        um_4_exp_diff = (um_4_exp - um_1_exp)/um_1_exp
        um_4_var_diff = (um_4_var - um_1_var)/um_1_var

        if order_A > 0 or mom_order_A > 0:
            alpha1_4_exp = np.array(ax_arr[2].lines[9].get_ydata())
            alpha1_4_var = np.array(ax_arr[2].lines[10].get_ydata()) - alpha1_4_exp
            alpha1_4_exp_diff = (alpha1_4_exp - alpha1_1_exp)/alpha1_1_exp
            alpha1_4_var_diff = (alpha1_4_var - alpha1_1_var)/alpha1_1_var
        
        if order_A > 1 or mom_order_A > 1:
            alpha2_4_exp = np.array(ax_arr[3].lines[9].get_ydata())
            alpha2_4_var = np.array(ax_arr[3].lines[10].get_ydata()) - alpha2_4_exp
            alpha2_4_exp_diff = (alpha2_4_exp - alpha2_1_exp)/alpha2_1_exp
            alpha2_4_var_diff = (alpha2_4_var - alpha2_1_var)/alpha2_1_var

plt.show()


if show_diff == True:
    order = min(order_A, mom_order_A)
    if order == 0:
        fig, (ax1, ax2) = plt.subplots(1,2)

    elif order == 1:
        fig, (ax1, ax2, ax3) = plt.subplots(1,3)

    elif order == 2:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

    ax1.plot(h_x, h_2_exp_diff, label='{}, mean'.format(label_B), color=color_1_B, linestyle=ls_B, linewidth=lw_B)
    if SG_order_B > 0:
        ax1.plot(h_x, h_2_var_diff, label='{}, std'.format(label_B), color=color_2_B, linestyle=ls_B, linewidth=lw_B)
    if no_lines > 2:
        ax1.plot(h_x, h_3_exp_diff, label='{}, mean'.format(label_C), color=color_1_C, linestyle=ls_C, linewidth=lw_C)
        ax1.plot(h_x, h_3_var_diff, label='{}, std'.format(label_C), color=color_2_C, linestyle=ls_C, linewidth=lw_C)
    if no_lines > 3:
        ax1.plot(h_x, h_4_exp_diff, label='{}, mean'.format(label_D), color=color_1_D, linestyle=ls_D, linewidth=lw_D)
        ax1.plot(h_x, h_4_var_diff, label='{}, std'.format(label_D), color=color_2_D, linestyle=ls_D, linewidth=lw_D)
    ax1.set_xlabel('x')
    ax1.legend()
    ax1.grid()

    ax2.plot(um_x, um_2_exp_diff, label='{}, mean'.format(label_B), color=color_1_B, linestyle=ls_B, linewidth=lw_B)
    if SG_order_B > 0:
        ax2.plot(um_x, um_2_var_diff, label='{}, std'.format(label_B), color=color_2_B, linestyle=ls_B, linewidth=lw_B)
    if no_lines > 2:
        ax2.plot(um_x, um_3_exp_diff, label='{}, mean'.format(label_C), color=color_1_C, linestyle=ls_C, linewidth=lw_C)
        ax2.plot(um_x, um_3_var_diff, label='{}, std'.format(label_C), color=color_2_C, linestyle=ls_C, linewidth=lw_C)
    if no_lines > 3:
        ax2.plot(um_x, um_4_exp_diff, label='{}, mean'.format(label_D), color=color_1_D, linestyle=ls_D, linewidth=lw_D)
        ax2.plot(um_x, um_4_var_diff, label='{}, std'.format(label_D), color=color_2_D, linestyle=ls_D, linewidth=lw_D)
    ax2.set_xlabel('x')
    ax2.legend()
    ax2.grid()
    
    if order > 0:
        ax3.plot(alpha1_x, alpha1_2_exp_diff, label='{}, mean'.format(label_B), color=color_1_B, linestyle=ls_B, linewidth=lw_B)
        if SG_order_B > 0:
            ax3.plot(alpha1_x, alpha1_2_var_diff, label='{}, std'.format(label_B), color=color_2_B, linestyle=ls_B, linewidth=lw_B)
        if no_lines > 2:
            ax3.plot(alpha1_x, alpha1_3_exp_diff, label='{}, mean'.format(label_C), color=color_1_C, linestyle=ls_C, linewidth=lw_C)
            ax3.plot(alpha1_x, alpha1_3_var_diff, label='{}, std'.format(label_C), color=color_2_C, linestyle=ls_C, linewidth=lw_C)
        if no_lines > 3:
            ax3.plot(alpha1_x, alpha1_4_exp_diff, label='{}, mean'.format(label_D), color=color_1_D, linestyle=ls_D, linewidth=lw_D)
            ax3.plot(alpha1_x, alpha1_4_var_diff, label='{}, std'.format(label_D), color=color_2_D, linestyle=ls_D, linewidth=lw_D)
        ax3.set_xlabel('x')
        ax3.legend()
        ax3.grid()
    
    if order > 1:
        ax4.plot(alpha2_x, alpha2_2_exp_diff, label='{}, mean'.format(label_B), color=color_1_B, linestyle=ls_B, linewidth=lw_B)
        if SG_order_B > 0:
            ax4.plot(alpha2_x, alpha2_2_var_diff, label='{}, std'.format(label_B), color=color_2_B, linestyle=ls_B, linewidth=lw_B)
        if no_lines > 2:
            ax4.plot(alpha2_x, alpha2_3_exp_diff, label='{}, mean'.format(label_C), color=color_1_C, linestyle=ls_C, linewidth=lw_C)
            ax4.plot(alpha2_x, alpha2_3_var_diff, label='{}, std'.format(label_C), color=color_2_C, linestyle=ls_C, linewidth=lw_C)
        if no_lines > 3:
            ax4.plot(alpha2_x, alpha2_4_exp_diff, label='{}, mean'.format(label_D), color=color_1_D, linestyle=ls_D, linewidth=lw_D)
            ax4.plot(alpha2_x, alpha2_4_var_diff, label='{}, std'.format(label_D), color=color_2_D, linestyle=ls_D, linewidth=lw_D)
        ax4.set_xlabel('x')
        ax4.legend()
        ax4.grid()
    
    if show_only_one == 1:
        ax2.remove()
        if mom_order_A > 0 or order_A > 0:
            ax3.remove()
        if mom_order_A > 1 or order_A > 1:
            ax4.remove()
        ax1.set_position([0.18,0.12,0.81,0.87])

    elif show_only_one == 2:
        ax1.remove()
        if mom_order_A > 0 or order_A > 0:
            ax3.remove()
        if mom_order_A > 1 or order_A > 1:
            ax4.remove()
        ax2.set_position([0.15,0.12,0.84,0.87])

    elif show_only_one == 3:
        ax1.remove()
        ax2.remove()
        if mom_order_A > 1 or order_A > 1:
            ax4.remove()
        ax3.set_position([0.15,0.12,0.84,0.87])

    elif show_only_one == 4:
        ax1.remove()
        ax2.remove()
        ax3.remove()
        ax4.set_position([0.12,0.12,0.87,0.87])

    plt.show()