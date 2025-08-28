import Simulation
import PDE
import Mesh
import SpatialDiscretization
import TimeIntegration
import numpy as np
import configparser
import timeit

def main():

    config = configparser.ConfigParser()
    config.read('config.txt')
    pde_information = config['pde_information']
    grid_information = config['grid_information']
    numerical_method_information = config['numerical_method_information']

    if pde_information['pde_type'] == 'SWME1D':
        _pde = PDE.SWME1D(pde_information['initialCondition'],
                          pde_information.getfloat('viscosity'),
                          pde_information.getfloat('slipLength'),
                          hyperbolic=False)
    elif pde_information['pde_type'] == 'HSWME1D':
        _pde = PDE.SWME1D(pde_information['initialCondition'],
                          pde_information.getfloat('viscosity'),
                          pde_information.getfloat('slipLength'),
                          hyperbolic=True)
    elif pde_information['pde_type'] == 'VegetationSWME1D':
        _pde = PDE.VegetationSWME1D(pde_information['initialCondition'],
                                    pde_information.getfloat('viscosity'),
                                    pde_information.getfloat('slipLength'),
                                    False,
                                    1,
                                    1,
                                    1)
    elif pde_information['pde_type'] == 'SGSWME1D' and numerical_method_information.getboolean('stochasticGalerkin') and not numerical_method_information.getboolean('spatiallyAdaptive') and not numerical_method_information.getboolean('monteCarlo'):
        _pde = PDE.SGSWME1D(pde_information['initialCondition'],
                            pde_information['distr'],
                            pde_information.getfloat('mu'),
                            pde_information.getfloat('sigma'),
                            pde_information.getfloat('slipLength'),
                            hyperbolic=False)
    
    elif pde_information['pde_type'] == 'HSGSWME1D' and numerical_method_information.getboolean('stochasticGalerkin') and not numerical_method_information.getboolean('spatiallyAdaptive') and not numerical_method_information.getboolean('monteCarlo'):
        _pde = PDE.SGSWME1D(pde_information['initialCondition'],
                            pde_information['distr'],
                            pde_information.getfloat('mu'),
                            pde_information.getfloat('sigma'),
                            pde_information.getfloat('slipLength'),
                            hyperbolic=True)
    
    elif pde_information['pde_type'] == 'SGSWME1D' or pde_information['pde_type'] == 'HSGSWME1D':
        print("pde_type can only be SGSWME1D if stochasticGalerkin is True and spatiallyAdaptive and monteCarlo are False")
    
    else:
        print('This pde_type is not implemented yet')
    
    ##########################################################################

    if numerical_method_information['fvm_type'] == 'PVM':
        if numerical_method_information['pvm'] == 'PRICE':
            _spatialDiscretization = SpatialDiscretization.PRICE()
        else:
            print('This pvm method is not implemented yet')
    else:
        print('This finite volume type is not implemented yet')


    #########################################################################

    if pde_information.getboolean('1D'):

        _mesh = Mesh.UniformRectangularMesh1D([grid_information.getfloat('x1boundary'),grid_information.getfloat('x2boundary')],
                                               grid_information.getint('resolutionX')) #TODO: Implement different grids
        
        if not numerical_method_information.getboolean('stochasticGalerkin'):
            if numerical_method_information.getboolean('spatiallyAdaptive'):
                boundaryInterfaces = numerical_method_information['boundaryInterfaces']
                boundaryInterfaces = [float(boundaryInterface) for boundaryInterface in boundaryInterfaces.split(',')]
                orders = numerical_method_information['orders']
                orders = [int(order) for order in orders.split(',')]

                _simulation = Simulation.SpatiallyAdaptiveSimulation1D(
                    [float(boundaryInterface) for boundaryInterface in numerical_method_information['boundaryInterfaces'].split(',')],
                    [int(order) for order in numerical_method_information['orders'].split(',')],
                    _pde,
                    _mesh,
                    numerical_method_information['boundaryCondition'],
                    pde_information['initialCondition'],
                    _spatialDiscretization)
            
            else:
                _simulation = Simulation.ClassicalSimulation1D(
                    numerical_method_information.getint('order'),
                    _pde,
                    _mesh,
                    numerical_method_information['boundaryCondition'],
                    pde_information['initialCondition'],
                    _spatialDiscretization)
        
        else:
            _simulation = Simulation.ClassicalGalerkinSimulation1D(
                    numerical_method_information.getint('momOrder'),
                    numerical_method_information.getint('SGOrder'),
                    _pde,
                    _mesh,
                    numerical_method_information['boundaryCondition'],
                    pde_information['initialCondition'],
                    _spatialDiscretization)

        start = timeit.default_timer()
        
        if not numerical_method_information.getboolean('monteCarlo'):
            data_array = _simulation.run_simulation(numerical_method_information.getfloat('t_end'))
        
        elif numerical_method_information.getboolean('monteCarlo') and pde_information['pde_type'] == 'SWME1D' and not numerical_method_information.getboolean('stochasticGalerkin') and not numerical_method_information.getboolean('spatiallyAdaptive'):
            data_array = np.zeros((numerical_method_information.getint('n_MC'), grid_information.getint('resolutionX'), numerical_method_information.getint('order') + 3))
            if pde_information['distr'] == "normal":
                for n in range(numerical_method_information.getint('n_MC')):
                    data_array[n,:,:] = _simulation.run_simulation(numerical_method_information.getfloat('t_end'), viscosity = np.random.normal(pde_information.getfloat('mu'), pde_information.getfloat('sigma')))
            elif pde_information['distr'] == "uniform":
                for n in range(numerical_method_information.getint('n_MC')):
                    data_array[n,:,:] = _simulation.run_simulation(numerical_method_information.getfloat('t_end'), viscosity = np.random.uniform(pde_information.getfloat('mu') - pde_information.getfloat('sigma'), pde_information.getfloat('mu') + pde_information.getfloat('sigma')))
            else:
                print("This distribution is not implemented yet for the Monte Carlo loop")
        
        else:
            print("Monte Carlo loop can only be used on pde_type SWME1D and cannot be used in combination with stochasticGalerkin and/or spatiallyAdaptive = True")
        
        stop = timeit.default_timer()
        print('Time: ', stop - start)

        if numerical_method_information.getboolean('spatiallyAdaptive'):
            np.save("Data\data_{0}_orders={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_information['pde_type'], numerical_method_information['orders'], pde_information.getfloat('viscosity'), pde_information.getfloat('slipLength'), pde_information['initialCondition'], numerical_method_information.getfloat('t_end')), data_array)
        
        elif numerical_method_information.getboolean('monteCarlo'):
            np.save("Data\data_{0}_{1}_order={2}_N={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_information['pde_type'], pde_information['distr'], numerical_method_information.getint('order'), numerical_method_information.getint('n_MC'), pde_information.getfloat('mu'), pde_information.getfloat('sigma'), pde_information.getfloat('slipLength'), pde_information['initialCondition'], numerical_method_information.getfloat('t_end')), data_array)
        
        elif numerical_method_information.getboolean('stochasticGalerkin'):
            np.save("Data\data_{0}_{1}_MO={2},SO={3}_mu={4}_sigma={5}_lambda={6}_IC={7}_T={8}.npy".format(pde_information['pde_type'], pde_information['distr'], numerical_method_information.getint('momOrder'),numerical_method_information.getint('SGOrder'), pde_information.getfloat('mu'), pde_information.getfloat('sigma'), pde_information.getfloat('slipLength'), pde_information['initialCondition'], numerical_method_information.getfloat('t_end')), data_array)
        
        else:
            np.save("Data\data_{0}_order={1}_nu={2}_lambda={3}_IC={4}_T={5}.npy".format(pde_information['pde_type'], numerical_method_information['order'], pde_information.getfloat('viscosity'), pde_information.getfloat('slipLength'), pde_information['initialCondition'], numerical_method_information.getfloat('t_end')), data_array)
    
    else:
        print('2D not implemented yet')

if __name__ == '__main__':
    main()