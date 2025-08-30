import numpy as np
import matplotlib.pyplot as plt
import PDE

def visualisation(pde_type: PDE.PDE,
                  spatially_adaptive: bool,
                  monte_carlo: bool,
                  stochastic_galerkin: bool,
                  resolutionX: int,
                  x1boundary: float,
                  x2boundary: float,
                  IC: str,
                  data_array: np.array,
                  title: str,
                  label: str,
                  **kwargs):
    
    z = np.linspace(0,1,100)
    color_1 = kwargs["color_1"] if "color_1" in kwargs else 'tab:green'
    color_2 = kwargs["color_2"] if "color_2" in kwargs else 'greenyellow'
    color_3 = kwargs["color_3"] if "color_3" in kwargs else 'forestgreen'

    if not monte_carlo and not stochastic_galerkin:
        if spatially_adaptive:
            order = kwargs["orders"] if "orders" in kwargs else print("Orders should be specified for spatially adaptive plotting")
            max_order = np.max([int(order_i) for order_i in order.split(',')])
            velocity_profile = pde_type.compute_vertical_velocity_profile(max_order, data_array, z)
        
        else:
            order = kwargs["order"] if "order" in kwargs else print("Order should be specified for classical plotting")
            max_order = order
            velocity_profile = pde_type.compute_vertical_velocity_profile(order, data_array, z)
        
        if max_order == 0:
            fig, (ax1, ax2, ax3) = plt.subplots(1,3)
        elif max_order == 1:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2) 
        elif max_order == 2:
            fig, ((ax1, ax2, ax3), (ax4, ax5, _)) = plt.subplots(2,3)
        elif max_order == 3:
            fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3)
        elif max_order == 4:
            fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (_, ax7, _)) = plt.subplots(3,3)
        elif max_order == 5:
            fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, _)) = plt.subplots(2,3)
        elif max_order == 6:
            fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(2,3)
        else:
            print("Spatially adaptive and deterministic classical simulations cannot be plotted for orders > 6 at this moment")
        
        ax1.set_title('Vertical Velocity Profile at x = {0}, x = {1} and x = {2} for {3}'.format(x1boundary + (x2boundary - x1boundary)/4, x1boundary + (x2boundary - x1boundary)/2, x1boundary + 3*(x2boundary - x1boundary)/4, title))
        ax1.plot(velocity_profile[int(resolutionX/2),:], z)
        ax1.set_ylabel('z')
        
        ax2.set_title('$h$ for {}'.format(title))
        ax2.plot(data_array[:,0], data_array[:,1], '-', label='$h$ for {}'.format(label), color=color_1)
        ax2.set_xlabel('x')
        ax2.legend()
        
        ax3.set_title('$u_m$ for {}'.format(title))
        ax3.plot(data_array[:,0], data_array[:,2], '-', label='$u_m$ for {}'.format(label), color=color_1)
        ax3.set_xlabel('x')
        ax3.legend()

        if max_order > 0:
            ax4.set_title(r'$\alpha_1$ for {}'.format(title))
            ax4.plot(data_array[:,0], data_array[:,3], '-', label=r'$\alpha_1$ for {}'.format(label), color=color_1)
            ax4.set_xlabel('x')
            ax4.legend()
        
        if max_order > 1:
            ax5.set_title(r'$\alpha_2$ for {}'.format(title))
            ax5.plot(data_array[:,0], data_array[:,4], '-', label=r'$\alpha_2$ for {}'.format(label), color=color_1)
            ax5.set_xlabel('x')
            ax5.legend()
        
        if max_order > 2:
            ax6.set_title(r'$\alpha_3$ for {}'.format(title))
            ax6.plot(data_array[:,0], data_array[:,5], '-', label=r'$\alpha_3$ for {}'.format(label), color=color_1)
            ax6.set_xlabel('x')
            ax6.legend()
        
        if max_order > 3:
            ax6.set_title(r'$\alpha_4$ for {}'.format(title))
            ax6.plot(data_array[:,0], data_array[:,6], '-', label=r'$\alpha_4$ for {}'.format(label), color=color_1)
            ax6.set_xlabel('x')
            ax6.legend()

        if max_order > 4:
            ax6.set_title(r'$\alpha_5$ for {}'.format(title))
            ax6.plot(data_array[:,0], data_array[:,7], '-', label=r'$\alpha_5$ for {}'.format(label), color=color_1)
            ax6.set_xlabel('x')
            ax6.legend()

        if max_order > 5:
            ax6.set_title(r'$\alpha_6$ for {}'.format(title))
            ax6.plot(data_array[:,0], data_array[:,8], '-', label=r'$\alpha_6$ for {}'.format(label), color=color_1)
            ax6.set_xlabel('x')
            ax6.legend()
    
    elif monte_carlo or stochastic_galerkin and not spatially_adaptive:
        if monte_carlo and stochastic_galerkin:
            print("Stochastic Galerkin and a monte carlo loop cannot (and should not) be used at the same time")
        
        elif monte_carlo:
            n_MC = kwargs["n_MC"] if "n_MC" in kwargs else print("Number of Monte Carlo samples should be given for Monte Carlo plotting")
            order = kwargs["order"] if "order" in kwargs else print("Order should be specified for monte carlo loop plotting")
            h_exp = np.average([data_array[n,:,1] for n in range(n_MC)], axis=0)
            h_var = np.var([data_array[n,:,1] for n in range(n_MC)], axis=0)
            um_exp = np.average([data_array[n,:,2] for n in range(n_MC)], axis=0)
            um_var = np.var([data_array[n,:,2] for n in range(n_MC)], axis=0)
        
            velocity_profile_exp = np.zeros((len(um_exp), len(z)))
            velocity_profile_var = np.zeros((len(um_exp), len(z)))
            
            if order == 0:
                for i in range(len(um_exp)):
                    velocity_profile_exp[i,:] = um_exp[i]*(np.ones(len(z)))
                    velocity_profile_var[i,:] = um_var[i]*(np.ones(len(z)))
                fig, (ax1, ax2, ax3) = plt.subplots(1,3)

            elif order == 1:
                alpha1_exp = np.average([data_array[n,:,3] for n in range(n_MC)], axis=0)
                alpha1_var = np.var([data_array[n,:,3] for n in range(n_MC)], axis=0)
                for i in range(len(um_exp)):
                    velocity_profile_exp[i,:] = um_exp[i]*(np.ones(len(z))) + alpha1_exp[i]*(np.ones(len(z)) - 2*z)
                    velocity_profile_var[i,:] = um_var[i]*(np.ones(len(z))) + alpha1_var[i]*(np.ones(len(z)) - 2*z)**2
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
                
                ax4.set_title(r'$\alpha_1$ for {}'.format(title))
                ax4.plot(data_array[0,:,0], alpha1_exp, '-', label=r'$\alpha_1$ for {}'.format(label), color=color_1)
                ax4.fill_between(data_array[0,:,0], alpha1_exp - np.sqrt(alpha1_var), alpha1_exp + np.sqrt(alpha1_var),label='$1\sigma$ for {}'.format(label), color=color_1, alpha=0.2)
                ax4.set_xlabel('x')
                ax4.legend()
            
            else:
                print("This order is not implemented yet for the Monte Carlo loop")
            
            
            ax2.set_title('$h$ for {}'.format(title))
            ax3.set_title('$u_m$ for {}'.format(title))
            x = data_array[0,:,0]
        
        else:
            mom_order = kwargs["mom_order"] if "mom_order" in kwargs else print("Moment order should be specified for stochastic Galerkin projection plotting")
            SG_order  = kwargs["SG_order"]  if "SG_order"  in kwargs else print("Stochastic Galerkin order should be specified for stochastic Galerkin projection plotting")
            velocity_profile_exp, velocity_profile_var = pde_type.compute_vertical_velocity_profile(mom_order, SG_order, data_array[:,1:], z)
            h_exp, h_var, um_exp, um_var, alpha1_exp, alpha1_var = pde_type.compute_exp_and_var(mom_order, SG_order, data_array[:,1:], True)
            
            if mom_order == 0:
                fig, (ax1, ax2, ax3) = plt.subplots(1,3)
            
            elif mom_order == 1:
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
                ax4.set_title(r'$\alpha_1$ for {}'.format(title))
                ax4.plot(data_array[:,0], alpha1_exp, '-', label=r'$\alpha_1$ for {}'.format(label), color=color_1)
                ax4.fill_between(data_array[:,0], alpha1_exp - np.sqrt(alpha1_var), alpha1_exp + np.sqrt(alpha1_var),label='$1\sigma$ for {}'.format(label), color=color_1, alpha=0.2)
                ax4.set_xlabel('x')
                ax4.legend()
            
            else:
                print("This moment order is not implemented yet for the SGSWME1D")
            
            ax2.set_title('$h$ for {}'.format(title))
            ax3.set_title('$u_m$ for {}'.format(title))
            x = data_array[:,0]
                
        ax1.set_title('Vertical Velocity Profile at x = {0}, x = {1} and x = {2} for {3}'.format(x1boundary + (x2boundary - x1boundary)/4, x1boundary + (x2boundary - x1boundary)/2, x1boundary + 3*(x2boundary - x1boundary)/4, title))
        ax1.plot(velocity_profile_exp[int(resolutionX/4),:], z, label='$u$ at x=-1 for {}'.format(label), color=color_2)
        ax1.fill_betweenx(z, velocity_profile_exp[int(resolutionX/4),:] - np.sqrt(velocity_profile_var[int(resolutionX/4),:]), velocity_profile_exp[int(resolutionX/4),:] + np.sqrt(velocity_profile_var[int(resolutionX/4),:]), label='$1\sigma$ at x=-1 for {}'.format(label), color=color_2, alpha=0.2)
        ax1.plot(velocity_profile_exp[int(resolutionX/2),:], z, label='$u$ at x=0 for {}'.format(label), color=color_1)
        ax1.fill_betweenx(z, velocity_profile_exp[int(resolutionX/2),:] - np.sqrt(velocity_profile_var[int(resolutionX/2),:]), velocity_profile_exp[int(resolutionX/2),:] + np.sqrt(velocity_profile_var[int(resolutionX/2),:]), label='$1\sigma$ at x=0 for {}'.format(label), color=color_1, alpha=0.2)
        ax1.plot(velocity_profile_exp[int(3*resolutionX/4),:], z, label='$u$ at x=1 for {}'.format(label), color=color_3)
        ax1.fill_betweenx(z, velocity_profile_exp[int(3*resolutionX/4),:] - np.sqrt(velocity_profile_var[int(3*resolutionX/4),:]), velocity_profile_exp[int(3*resolutionX/4),:] + np.sqrt(velocity_profile_var[int(3*resolutionX/4),:]), label='$1\sigma$ at x=1 for {}'.format(label), color=color_3, alpha=0.2)
        ax1.set_ylabel('z')
        ax1.legend()
            
        ax2.plot(x, h_exp, '-', label='$h$ for {}'.format(label), color=color_1)
        ax2.fill_between(x, h_exp - np.sqrt(h_var), h_exp + np.sqrt(h_var), label='$1\sigma$ for {}'.format(label), color=color_1, alpha=0.2)
        ax2.set_xlabel('x')
        ax2.legend()
            
        ax3.plot(x, um_exp, '-', label='$u_m$ for {}'.format(label), color=color_1)
        ax3.fill_between(x, um_exp - np.sqrt(um_var), um_exp + np.sqrt(um_var), label='$1\sigma$ for {}'.format(label), color=color_1, alpha=0.2)
        ax3.set_xlabel('x')
        ax3.legend()

    
    else:
        print("Spatial adaptivity is not available for a monte carlo loop or stochastic Galerkin yet")
    
    try:
        return [ax1, ax2, ax3, ax4, ax5, ax6]
    except NameError:
        try:
            return [ax1, ax2, ax3, ax4, ax5]
        except NameError:
            try:
                return [ax1, ax2, ax3, ax4]
            except NameError:
                return [ax1, ax2, ax3]


def visualisation_add(ax_arr,
                      pde_type: PDE.PDE,
                      spatially_adaptive: bool,
                      monte_carlo: bool,
                      stochastic_galerkin: bool,
                      resolutionX: int,
                      x1boundary: float,
                      x2boundary: float,
                      IC: str,
                      data_array: np.array,
                      label: str,
                      **kwargs):
    z = np.linspace(0,1,100)
    color_1 = kwargs["color_1"] if "color_1" in kwargs else 'tab:red'
    color_2 = kwargs["color_2"] if "color_2" in kwargs else 'lightcoral'
    color_3 = kwargs["color_3"] if "color_3" in kwargs else 'firebrick'

    if not monte_carlo and not stochastic_galerkin:
        if spatially_adaptive:
            order = kwargs["orders"] if "orders" in kwargs else print("Orders should be specified for spatially adaptive plotting")
            max_order = np.max([int(order_i) for order_i in order.split(',')])
            velocity_profile = pde_type.compute_vertical_velocity_profile(max_order, data_array, z)
        
        else:
            order = kwargs["order"] if "order" in kwargs else print("Order should be specified for classical plotting")
            max_order = order
            velocity_profile = pde_type.compute_vertical_velocity_profile(order, data_array, z)
        
        
        ax_arr[0].plot(velocity_profile[int(resolutionX/2),:], z)
        ax_arr[1].plot(data_array[:,0], data_array[:,1], '-', label='$h$ for {}'.format(label), color=color_1)
        ax_arr[2].plot(data_array[:,0], data_array[:,2], '-', label='$u_m$ for {}'.format(label), color=color_1)

        if max_order > 0:
            ax_arr[3].plot(data_array[:,0], data_array[:,3], '-', label=r'$\alpha_1$ for {}'.format(label), color=color_1)
        
        if max_order > 1:
            ax_arr[4].plot(data_array[:,0], data_array[:,4], '-', label=r'$\alpha_2$ for {}'.format(label), color=color_1)
        
        if max_order > 2:
            ax_arr[5].plot(data_array[:,0], data_array[:,5], '-', label=r'$\alpha_3$ for {}'.format(label), color=color_1)
        
        if max_order > 3:
            ax_arr[6].plot(data_array[:,0], data_array[:,6], '-', label=r'$\alpha_4$ for {}'.format(label), color=color_1)

        if max_order > 4:
            ax_arr[7].plot(data_array[:,0], data_array[:,7], '-', label=r'$\alpha_5$ for {}'.format(label), color=color_1)

        if max_order > 5:
            ax_arr[8].plot(data_array[:,0], data_array[:,8], '-', label=r'$\alpha_6$ for {}'.format(label), color=color_1)
    
    elif monte_carlo or stochastic_galerkin and not spatially_adaptive:
        if monte_carlo and stochastic_galerkin:
            print("Stochastic Galerkin and a monte carlo loop cannot (and should not) be used at the same time")
        
        elif monte_carlo:
            n_MC = kwargs["n_MC"] if "n_MC" in kwargs else print("Number of Monte Carlo samples should be given for Monte Carlo plotting")
            order = kwargs["order"] if "order" in kwargs else print("Order should be specified for monte carlo loop plotting")
            h_exp = np.average([data_array[n,:,1] for n in range(n_MC)], axis=0)
            h_var = np.var([data_array[n,:,1] for n in range(n_MC)], axis=0)
            um_exp = np.average([data_array[n,:,2] for n in range(n_MC)], axis=0)
            um_var = np.var([data_array[n,:,2] for n in range(n_MC)], axis=0)
        
            velocity_profile_exp = np.zeros((len(um_exp), len(z)))
            velocity_profile_var = np.zeros((len(um_exp), len(z)))
            
            if order == 0:
                for i in range(len(um_exp)):
                    velocity_profile_exp[i,:] = um_exp[i]*(np.ones(len(z)))
                    velocity_profile_var[i,:] = um_var[i]*(np.ones(len(z)))

            elif order == 1:
                alpha1_exp = np.average([data_array[n,:,3] for n in range(n_MC)], axis=0)
                alpha1_var = np.var([data_array[n,:,3] for n in range(n_MC)], axis=0)
                for i in range(len(um_exp)):
                    velocity_profile_exp[i,:] = um_exp[i]*(np.ones(len(z))) + alpha1_exp[i]*(np.ones(len(z)) - 2*z)
                    velocity_profile_var[i,:] = um_var[i]*(np.ones(len(z))) + alpha1_var[i]*(np.ones(len(z)) - 2*z)**2
                
                ax_arr[3].plot(data_array[0,:,0], alpha1_exp, '-', label=r'$\alpha_1$ for {}'.format(label), color=color_1)
                ax_arr[3].fill_between(data_array[0,:,0], alpha1_exp - np.sqrt(alpha1_var), alpha1_exp + np.sqrt(alpha1_var),label='$1\sigma$ for {}'.format(label), color=color_1, alpha=0.2)
                ax_arr[3].legend()
            
            else:
                print("This order is not implemented yet for the Monte Carlo loop")
            
            x = data_array[0,:,0]
        
        else:
            mom_order = kwargs["mom_order"] if "mom_order" in kwargs else print("Moment order should be specified for stochastic Galerkin projection plotting")
            SG_order  = kwargs["SG_order"]  if "SG_order"  in kwargs else print("Stochastic Galerkin order should be specified for stochastic Galerkin projection plotting")
            velocity_profile_exp, velocity_profile_var = pde_type.compute_vertical_velocity_profile(mom_order, SG_order, data_array[:,1:], z)
            h_exp, h_var, um_exp, um_var, alpha1_exp, alpha1_var = pde_type.compute_exp_and_var(mom_order, SG_order, data_array[:,1:], True)
            
            if mom_order == 1:
                ax_arr[3].plot(data_array[:,0], alpha1_exp, '-', label=r'$\alpha_1$ for {}'.format(label), color=color_1)
                ax_arr[3].fill_between(data_array[:,0], alpha1_exp - np.sqrt(alpha1_var), alpha1_exp + np.sqrt(alpha1_var),label='$1\sigma$ for {}'.format(label), color=color_1, alpha=0.2)
                ax_arr[3].legend()

            else:
                print("This moment order is not implemented yet for the SGSWME1D")
            
            x = data_array[:,0]
                
        ax_arr[0].plot(velocity_profile_exp[int(resolutionX/4),:], z, '-', label='$u$ at x=-1 for {}'.format(label), color=color_2)
        ax_arr[0].fill_betweenx(z, velocity_profile_exp[int(resolutionX/4),:] - np.sqrt(velocity_profile_var[int(resolutionX/4),:]), velocity_profile_exp[int(resolutionX/4),:] + np.sqrt(velocity_profile_var[int(resolutionX/4),:]), label='$1\sigma$ at x=-1 for {}'.format(label), color=color_2, alpha=0.2)
        ax_arr[0].plot(velocity_profile_exp[int(resolutionX/2),:], z, '-', label='$u$ at x=0 for {}'.format(label), color=color_1)
        ax_arr[0].fill_betweenx(z, velocity_profile_exp[int(resolutionX/2),:] - np.sqrt(velocity_profile_var[int(resolutionX/2),:]), velocity_profile_exp[int(resolutionX/2),:] + np.sqrt(velocity_profile_var[int(resolutionX/2),:]), label='$1\sigma$ at x=0 for {}'.format(label), color=color_1, alpha=0.2)
        ax_arr[0].plot(velocity_profile_exp[int(3*resolutionX/4),:], z, '-', label='$u$ at x=1 for {}'.format(label), color=color_3)
        ax_arr[0].fill_betweenx(z, velocity_profile_exp[int(3*resolutionX/4),:] - np.sqrt(velocity_profile_var[int(3*resolutionX/4),:]), velocity_profile_exp[int(3*resolutionX/4),:] + np.sqrt(velocity_profile_var[int(3*resolutionX/4),:]), label='$1\sigma$ at x=1 for {}'.format(label), color=color_3, alpha=0.2)
        ax_arr[0].legend()
 
        ax_arr[1].plot(x, h_exp, '-', label='$h$ for {}'.format(label), color=color_1)
        ax_arr[1].fill_between(x, h_exp - np.sqrt(h_var), h_exp + np.sqrt(h_var), label='$1\sigma$ for {}'.format(label), color=color_1, alpha=0.2)
        ax_arr[1].legend()
            
        ax_arr[2].plot(x, um_exp, '-', label='$u_m$ for {}'.format(label), color=color_1)
        ax_arr[2].fill_between(x, um_exp - np.sqrt(um_var), um_exp + np.sqrt(um_var), label='$1\sigma$ for {}'.format(label), color=color_1, alpha=0.2)
        ax_arr[2].legend()

    
    else:
        print("Spatial adaptivity is not available for a monte carlo loop or stochastic Galerkin yet")
    
    return ax_arr