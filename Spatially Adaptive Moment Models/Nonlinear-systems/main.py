import simulation
import pde
import mesh
import spatialDiscretization
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import configparser
import timeit

def main():

    config = configparser.ConfigParser()
    config.read('/Users/vildarozek/Files/Fork/Untitled/Spatially Adaptive Moment Models/Nonlinear-systems/Config-files/config.txt')
    pde_information = config['pde_information']
    grid_information = config['grid_information']
    numerical_method_information = config['numerical_method_information']

    if pde_information['pde_type'] == 'SWME1D':
        _pde = pde.SWME1D(pde_information['initialCondition'],
                          pde_information.getfloat('viscosity'),
                          pde_information.getfloat('slipLength'),
                          hyperbolic=False)
    elif pde_information['pde_type'] == 'HSWME1D':
        _pde = pde.SWME1D(pde_information['initialCondition'],
                          pde_information.getfloat('viscosity'),
                          pde_information.getfloat('slipLength'),
                          hyperbolic=True)
    elif pde_information['pde_type'] == 'VegetationSWME1D':
        _pde = pde.VegetationSWME1D(pde_information['initialCondition'],
                                    pde_information.getfloat('viscosity'),
                                    pde_information.getfloat('slipLength'),
                                    False,
                                    0.008,
                                    1,
                                    264)
    else:
        print('PDE_type is not implemented yet')
    
    ##########################################################################

    if numerical_method_information['fvm_type'] == 'PVM':
        if numerical_method_information['pvm'] == 'PRICE':
            _spatialDiscretization = spatialDiscretization.PRICE()
        elif numerical_method_information['pvm'] == 'LF':
            _spatialDiscretization = spatialDiscretization.LF()
        else:
            print('this pvm method is not implemented yet')
    else:
        print('this finite volume type is not implemented yet')


    #########################################################################

    if pde_information.getboolean('1D'):

        _mesh = mesh.UniformRectangularMesh1D([grid_information.getfloat('x1boundary'),grid_information.getfloat('x2boundary')],
                                               grid_information.getint('resolutionX')) #TODO: Implement different grids

        if numerical_method_information['method'] == 'spatially_adaptive':
            boundaryInterfaces = numerical_method_information['boundaryInterfaces']
            boundaryInterfaces = [float(boundaryInterface) for boundaryInterface in boundaryInterfaces.split(',')]
            orders = numerical_method_information['orders']
            orders = [int(order) for order in orders.split(',')]

            _simulation = simulation.SpatiallyAdaptiveSimulation1D(
                [float(boundaryInterface) for boundaryInterface in numerical_method_information['boundaryInterfaces'].split(',')],
                [int(order) for order in numerical_method_information['orders'].split(',')],
                _pde,
                _mesh,
                numerical_method_information['boundaryCondition'],
                pde_information['initialCondition'],
                pde_information['breakdown_criterion'],
                _spatialDiscretization
            )
        
        elif numerical_method_information['method'] == 'classical':
            
            _simulation = simulation.ClassicalSimulation1D(
                numerical_method_information.getint('order'),
                _pde,
                _mesh,
                numerical_method_information['boundaryCondition'],
                pde_information['initialCondition'],
                _spatialDiscretization)
            
        elif numerical_method_information['method'] == 'micro_macro':
            
            _simulation = simulation.Micro_macro(
                numerical_method_information.getint('order'),
                _pde,
                _mesh,
                numerical_method_information['boundaryCondition'],
                pde_information['initialCondition'],
                _spatialDiscretization)

        start = timeit.default_timer()
        data_array = _simulation.run_simulation(numerical_method_information.getfloat('t_end'))
        stop = timeit.default_timer()
        print('Time: ', stop - start)
        data_frame = pd.DataFrame(data_array)
        #data_frame.to_csv('Data-processing/Results/test_LF.csv', index=False,header=False)

        z = np.linspace(0,1,100)
        if numerical_method_information.getboolean('spatiallyAdaptive'):
            velocity_profile = _pde.compute_vertical_velocity_profile(np.max([int(order) for order in numerical_method_information['orders'].split(',')]),
                                                                      data_array,
                                                                      z)
            number_of_variables = _simulation.max_number_of_variables
        else: 
            velocity_profile = _pde.compute_vertical_velocity_profile(numerical_method_information.getint('order'),
                                                                      data_array,
                                                                      z)
            number_of_variables = _simulation.number_of_variables

        plt.figure()
        plt.subplot(2,3,1)
        plt.plot(velocity_profile[200,:], z)
        plt.title('Velocity profile')

        plt.subplot(2,3,2)
        plt.plot(_mesh.cell_center_positions, data_array[:,1])
        plt.title('Height')

        plt.subplot(2,3,3)
        plt.plot(_mesh.cell_center_positions, data_array[:,2])
        plt.title('Velocity')

        plt.subplot(2,3,4)
        plt.plot(_mesh.cell_center_positions,_pde.compute_breakdown_criterion(data_array[:,1:],number_of_variables,'height_gradient',_mesh.resolution))
        plt.title('Height gradient')

        plt.subplot(2,3,5)
        plt.plot(_mesh.cell_center_positions,_pde.compute_breakdown_criterion(data_array[:,1:],number_of_variables,'momentum_gradient',_mesh.resolution))
        plt.title('Velocity gradient')

        plt.subplot(2,3,6)
        plt.plot(_mesh.cell_center_positions,_pde.compute_breakdown_criterion(data_array[:,1:],number_of_variables,'last_moment',_mesh.resolution))
        plt.title('Absolute value last moment')
        
        plt.show()
        
    else:
        print('2D not implemented yet')
    

if __name__ == '__main__':
    main()