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
                  **kwargs):
    
    z = np.linspace(0,1,100)

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
        
        ax1.set_title('Vertical velocity profile at x = {}'.format(x1boundary + (x2boundary - x1boundary)/2))
        ax1.plot(velocity_profile[int(resolutionX/2),:], z)
        ax1.set_ylabel('z')
        
        ax2.set_title('$h$ for order = {0}, IC = {1}'.format(order, IC))
        ax2.plot(data_array[:,0], data_array[:,1], '-', label='$h$', color='tab:red')
        ax2.set_xlabel('x')
        ax2.legend()
        
        ax3.set_title('$u_m$ for order = {0}, IC = {1}'.format(order, IC))
        ax3.plot(data_array[:,0], data_array[:,2], '-', label='$u_m$', color='tab:green')
        ax3.set_xlabel('x')
        ax3.legend()

        if max_order > 0:
            ax4.set_title(r'$\alpha_1$ for order = {0}, IC = {1}'.format(order, IC))
            ax4.plot(data_array[:,0], data_array[:,3], '-', label=r'$\alpha_1$', color='tab:orange')
            ax4.set_xlabel('x')
            ax4.legend()
        
        if max_order > 1:
            ax5.set_title(r'$\alpha_2$ for order = {0}, IC = {1}'.format(order, IC))
            ax5.plot(data_array[:,0], data_array[:,4], '-', label=r'$\alpha_2$', color='tab:orange')
            ax5.set_xlabel('x')
            ax5.legend()
        
        if max_order > 2:
            ax6.set_title(r'$\alpha_3$ for order = {0}, IC = {1}'.format(order, IC))
            ax6.plot(data_array[:,0], data_array[:,5], '-', label=r'$\alpha_3$', color='tab:orange')
            ax6.set_xlabel('x')
            ax6.legend()
        
        if max_order > 3:
            ax6.set_title(r'$\alpha_4$ for order = {0}, IC = {1}'.format(order, IC))
            ax6.plot(data_array[:,0], data_array[:,6], '-', label=r'$\alpha_4$', color='tab:orange')
            ax6.set_xlabel('x')
            ax6.legend()

        if max_order > 4:
            ax6.set_title(r'$\alpha_5$ for order = {0}, IC = {1}'.format(order, IC))
            ax6.plot(data_array[:,0], data_array[:,7], '-', label=r'$\alpha_5$', color='tab:orange')
            ax6.set_xlabel('x')
            ax6.legend()

        if max_order > 5:
            ax6.set_title(r'$\alpha_6$ for order = {0}, IC = {1}'.format(order, IC))
            ax6.plot(data_array[:,0], data_array[:,8], '-', label=r'$\alpha_6$', color='tab:orange')
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
                
                ax4.set_title(r'$\alpha_1$ for order = {0}, IC = {1}'.format(order, IC))
                ax4.plot(data_array[0,:,0], alpha1_exp, '-', label=r'$\alpha_1$', color='tab:orange')
                ax4.fill_between(data_array[0,:,0], alpha1_exp - np.sqrt(alpha1_var), alpha1_exp + np.sqrt(alpha1_var),label='$1\sigma$', color='tab:orange', alpha=0.2)
                ax4.set_xlabel('x')
                ax4.legend()
            
            else:
                print("This order is not implemented yet for the Monte Carlo loop")
            
            
            ax2.set_title('$h$ for order = {0}, IC = {1}'.format(order, IC))
            ax3.set_title('$u_m$ for order = {0}, IC = {1}'.format(order, IC))
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
                ax4.set_title(r'$\alpha_1$ for MO = {0}, SO = {1}, IC = {2}'.format(mom_order, SG_order, IC))
                ax4.plot(data_array[:,0], alpha1_exp, '-', label=r'$\alpha_1$', color='tab:orange')
                ax4.fill_between(data_array[:,0], alpha1_exp - np.sqrt(alpha1_var), alpha1_exp + np.sqrt(alpha1_var),label='$1\sigma$', color='tab:orange', alpha=0.2)
                ax4.set_xlabel('x')
                ax4.legend()
            
            else:
                print("This moment order is not implemented yet for the SGSWME1D")
            
            ax2.set_title('$h$ for MO = {0}, SO = {1}, IC = {2}'.format(mom_order, SG_order, IC))
            ax3.set_title('$u_m$ for MO = {0}, SO = {1}, IC = {2}'.format(mom_order, SG_order, IC))
            x = data_array[:,0]
                
        ax1.set_title('Vertical velocity profile at x = {}'.format(x1boundary + (x2boundary - x1boundary)/2))
        ax1.plot(velocity_profile_exp[int(resolutionX/2),:], z)
        ax1.fill_betweenx(z, velocity_profile_exp[int(resolutionX/2),:] - np.sqrt(velocity_profile_var[int(resolutionX/2),:]), velocity_profile_exp[int(resolutionX/2),:] + np.sqrt(velocity_profile_var[int(resolutionX/2),:]), label='$1\sigma$', color='tab:blue', alpha=0.2)
        ax1.set_ylabel('z')
        ax1.legend()
            
        ax2.plot(x, h_exp, '-', label='$h$', color='tab:red')
        ax2.fill_between(x, h_exp - np.sqrt(h_var), h_exp + np.sqrt(h_var), label='$1\sigma$', color='tab:red', alpha=0.2)
        ax2.set_xlabel('x')
        ax2.legend()
            
        ax3.plot(x, um_exp, '-', label='$u_m$', color='tab:green')
        ax3.fill_between(x, um_exp - np.sqrt(um_var), um_exp + np.sqrt(um_var), label='$1\sigma$', color='tab:green', alpha=0.2)
        ax3.set_xlabel('x')
        ax3.legend()

    
    else:
        print("Spatial adaptivity is not available for a monte carlo loop or stochastic Galerkin yet")
    
    plt.show()