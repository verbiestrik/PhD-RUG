from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
import pde
import mesh
import simulation

class Plotting(ABC):

    """
    This abstract class represents a plotting object (for the plotting of the simulation results).

    ...

    Attributes
    ----------
    pde_type : PDE
        The partial differential equation that has been simulated
    mesh : RectangularMesh
        The simulation mesh
    simulation : Simulation
        The simulation object
    
    Class methods
    -------------
    def __init__(self,pde_type):
        initializes the plotting object

    Abstract methods
    ---------------
    def plot(self):
        creates a plotting object and plots the simulation results
    """
    def __init__(self,
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 simulation: simulation.Simulation):

        """
        initializes the plotting object

        Parameters
        ------------
        pde_type : PDE
            the PDE model
        mesh : RectangularMesh
            the numerical simulation mesh
        simulation : Simulation
            the simulation object

        Returns
        --------        
        None

        """

        self.pde_type = pde_type 
        self.mesh = mesh
        self.simulation = simulation 

    @abstractmethod
    def plot(self,
             data_array: np.ndarray):
        """
        Creates a plot of the data listed in data_array

        Parameters
        ----------
        data_array : numpy array

        Returns
        -------
        None

        """
        pass

class SWME1DPlotClassical(Plotting):

    """
    This class represents a plotting object for the plotting of numerical results of the 1D SWME 
    of a classical simulation.

    ...

    Attributes
    ----------
    pde_type : SWME1D
        the 1D SWME object
    mesh : RectangularMesh
        The simulation mesh
    simulation : ClassicalSimulation1D
        The classical 1D simulation object

    Implemented methods from abstract parent class 'Plotting'
    ---------------------------------------------------------
    def plot(self):
        creates a plotting object and plots the simulation results

    Methods overriden from abstract parent class 'Plotting
    ------------------------------------------------------
    def __init__(self,pde_type,mesh,simulation):
        initializes the plotting object

    """

    def __init__(self,
                 pde_type: pde.SWME1D,
                 mesh: mesh.RectangularMesh,
                 simulation: simulation.ClassicalSimulation1D):
        """
        initializes the classical SWME1D plotting object

        Parameters
        ------------
        pde_type : SWME1D
            the SWME1D moment model
        mesh : RectangularMesh
            the numerical simulation mesh
        simulation : ClassicalSimulation1D
            the simulation object

        Returns
        --------        
        None

        """

        self.pde_type = pde_type 
        self.mesh = mesh
        self.simulation = simulation

    def plot(self,
             data_array: np.ndarray):

        z = np.linspace(0,1,100)

        velocity_profile = self.pde_type.compute_vertical_velocity_profile(self.simulation.order,
                                                                    data_array,
                                                                    z)
        order = self.simulation.order

        print('total mass = ',np.sum(data_array[:,1]*data_array[:,2]))

        plt.figure()
        plt.subplot(3,3,1)
        plt.plot(velocity_profile[np.floor_divide(self.mesh.resolution,2),:], z)
        plt.title('Velocity profile')

        plt.subplot(3,3,2)
        plt.plot(self.mesh.cell_center_positions, data_array[:,1])
        plt.title('Height')

        plt.subplot(3,3,3)
        plt.plot(self.mesh.cell_center_positions, data_array[:,2])
        plt.title('Velocity')

        k = 4
        for i in range(order):
            plt.subplot(3,3,k)
            plt.plot(self.mesh.cell_center_positions, data_array[:,3+i])
            plt.title('alpha_'+str(i+1))
            k += 1

        plt.show()

class SWME1DPlotAdaptive(Plotting):

    """
    This class represents a plotting object for the plotting of numerical results of the 1D SWME of an adaptive simulation.

    ...

    Attributes
    ----------
    pde_type : SWME1D
        the 1D SWME object
    mesh : RectangularMesh
        The simulation mesh
    simulation : SpatiallyAdaptiveSimulation1D
        The adaptive 1D simulation object

    Implemented methods from abstract parent class 'Plotting'
    ---------------------------------------------------------
    def plot(self,data_array):
        creates a plotting object and plots the simulation results

    Methods overriden from abstract parent class 'Plotting
    ------------------------------------------------------
    def __init__(self,pde_type,mesh,simulation):
        initializes the plotting object

    """

    def __init__(self,
                 pde_type: pde.SWME1D,
                 mesh: mesh.RectangularMesh,
                 simulation: simulation.ModelAdaptiveMomentSimulation1D,
                 type_model_error_estimator: str):
        """
        initializes the adaptive SWME1D plotting object

        Parameters
        ------------
        pde_type : SWME1D
            the SWME1D moment model
        mesh : RectangularMesh
            the numerical simulation mesh
        simulation : SpatiallyAdaptiveSimulation1D
            the adaptive 1D simulation object
        type_model_error_estimator : str
            the type of model error estimator that is used for the domain decomposition

        Returns
        --------        
        None

        """
        self.pde_type = pde_type 
        self.mesh = mesh
        self.simulation = simulation
        self.type_model_error_estimator = type_model_error_estimator

    def plot(self,
             data_array: np.ndarray):
        
        z = np.linspace(0,1,100)

        velocity_profile = self.pde_type.compute_vertical_velocity_profile(self.simulation.max_order,
                                                                    data_array,
                                                                    z)
        order = self.simulation.max_order

        print('total mass = ',np.sum(data_array[:,1]))

        plt.figure()
        plt.subplot(4,4,1)
        plt.plot(velocity_profile[np.floor_divide(self.mesh.resolution,2),:], z)
        plt.title('Velocity profile')

        plt.subplot(4,4,2)
        plt.plot(self.mesh.cell_center_positions, data_array[:,1])
        plt.title('Height')

        plt.subplot(4,4,3)
        plt.plot(self.mesh.cell_center_positions, data_array[:,2])
        plt.title('Velocity')

        k = 4
        for i in range(order):
            plt.subplot(4,4,k)
            plt.plot(self.mesh.cell_center_positions, data_array[:,3+i])
            plt.title('alpha_'+str(i))
            k += 1

        if self.type_model_error_estimator == 'heuristics_plus_discretization':

            plt.subplot(4,4,k)
            plt.plot(self.mesh.cell_center_positions,self.simulation.breakdown_estimators_coarsening[:,1])
            plt.title('Absolute value last moment')

            plt.subplot(4,4,k+1)
            plt.plot(self.mesh.cell_center_positions,self.simulation.breakdown_estimators_coarsening[:,-2])
            plt.title('transport residual')

            plt.subplot(4,4,k+2)
            plt.plot(self.mesh.cell_center_positions,self.simulation.breakdown_estimators_coarsening[:,-1])
            plt.title('source residual')

            plt.subplot(4,4,k+3)
            plt.plot(self.mesh.cell_center_positions,self.simulation.breakdown_estimators_refinement[:,1])
            plt.title('Height gradient')

            plt.subplot(4,4,k+4)
            # plt.plot(self.mesh.cell_center_positions[:-1],momentum_gradient)
            plt.plot(self.mesh.cell_center_positions,self.simulation.breakdown_estimators_refinement[:,2])
            plt.title('Velocity gradient')

            plt.subplot(4,4,k+5)
            plt.plot(self.mesh.cell_center_positions,self.simulation.breakdown_estimators_refinement[:,0])
            plt.title('source term last entry')

            plt.subplot(4,4,k+6)
            plt.plot(self.mesh.cell_center_positions,
                     order*self.simulation.breakdown_estimators_coarsening[:,-2]/max(np.max(self.simulation.breakdown_estimators_coarsening[:,-2]),0.001))
            plt.scatter(self.mesh.cell_center_positions,data_array[:,-1],s=order,color = 'hotpink')
            plt.title('orders vs coars. transp. res.')

            # plt.subplot(4,4,k+6)
            # plt.plot(self.mesh.cell_center_positions,data_array[:,1])
            # plt.scatter(self.mesh.cell_center_positions,
            #             (data_array[:,-1]*np.max(self.simulation.breakdown_estimators_coarsening[:,-2])+(5-data_array[:,-1])*np.min(self.simulation.breakdown_estimators_coarsening[:,-2])),s=5,color = 'hotpink')
            # plt.title('orders vs coars. transp. res.')

        plt.show()

class HME1DPlotClassical(Plotting):

    """
    This class represents a plotting object for the plotting of numerical results of the 1D HME of a classical simulation.

    ...

    Attributes
    ----------
    pde_type : HME1D
        the 1D HME object
    mesh : RectangularMesh
        The simulation mesh
    simulation : ClassicalSimulation1D
        The classical 1D simulation object

    Implemented methods from abstract parent class 'Plotting'
    ---------------------------------------------------------
    def plot(self,data_array):
        creates a plotting object and plots the simulation results

    Methods overriden from abstract parent class 'Plotting
    ------------------------------------------------------
    def __init__(self,pde_type,mesh,simulation):
        initializes the plotting object

    """

    def __init__(self,
                 pde_type: pde.HermiteMomentEquations,
                 mesh: mesh.RectangularMesh,
                 simulation: simulation.ClassicalSimulation1D):
        """
        initializes the classical HME1D plotting object

        Parameters
        ------------
        pde_type : HME1D
            the HME1D moment model
        mesh : RectangularMesh
            the numerical simulation mesh
        simulation : ClassicalSimulation1D
            the classical 1D simulation object

        Returns
        --------        
        None

        """
        self.pde_type = pde_type 
        self.mesh = mesh
        self.simulation = simulation

    def plot(self,
             data_array: np.ndarray):
        
        order = self.simulation.order

        plt.figure()

        plt.subplot(4,4,1)
        plt.plot(self.mesh.cell_center_positions, data_array[:,1])
        plt.title('Density')

        plt.subplot(4,4,2)
        plt.plot(self.mesh.cell_center_positions, data_array[:,2])
        plt.title('Velocity')

        plt.subplot(4,4,3)
        plt.plot(self.mesh.cell_center_positions, data_array[:,3])
        plt.title('Temperature')

        k = 4
        for i in range(3,min(order+1,11)):
            plt.subplot(4,4,k)
            plt.plot(self.mesh.cell_center_positions, data_array[:,i+1])
            plt.title('f_'+str(i))
            k += 1

        plt.show()

class HME1DPlotAdaptive(Plotting):

    """
    This class represents a plotting object for the plotting of numerical results of the 1D HME of an adaptive simulation.

    ...

    Attributes
    ----------
    pde_type : HME1D
        the 1D HME object
    mesh : RectangularMesh
        The simulation mesh
    simulation : SpatiallyAdaptiveSimulation1D
        The adaptive 1D simulation object

    Implemented methods from abstract parent class 'Plotting'
    ---------------------------------------------------------
    def plot(self,data_array):
        creates a plotting object and plots the simulation results

    Methods overriden from abstract parent class 'Plotting
    ------------------------------------------------------
    def __init__(self,pde_type,mesh,simulation):
        initializes the plotting object

    """

    def __init__(self,
                 pde_type: pde.HermiteMomentEquations,
                 mesh: mesh.RectangularMesh,
                 simulation: simulation.ModelAdaptiveMomentSimulation1D):
        """
        initializes the adaptive HME1D plotting object

        Parameters
        ------------
        pde_type : HME1D
            the HME1D moment model
        mesh : RectangularMesh
            the numerical simulation mesh
        simulation : SpatiallyAdaptiveSimulation1D
            the adaptive 1D simulation object

        Returns
        --------        
        None

        """
        self.pde_type = pde_type 
        self.mesh = mesh
        self.simulation = simulation

    def plot(self,
             data_array: np.ndarray):
        
        order = self.simulation.max_order

        plt.figure()

        plt.subplot(4,4,1)
        plt.plot(self.mesh.cell_center_positions, data_array[:,1])
        plt.title('Density')

        plt.subplot(4,4,2)
        plt.plot(self.mesh.cell_center_positions, data_array[:,2])
        plt.title('Velocity')

        plt.subplot(4,4,3)
        plt.plot(self.mesh.cell_center_positions, data_array[:,3])
        plt.title('Temperature')

        k = 4
        for i in range(3,min(order+1,13)):
            plt.subplot(4,4,k)
            plt.plot(self.mesh.cell_center_positions, data_array[:,i+1])
            plt.title('f_'+str(i))
            k += 1

        plt.subplot(4,4,k)
        plt.plot(self.mesh.cell_center_positions, self.simulation.breakdown_estimators_coarsening[:,0])
        plt.title('Decrease estimator')

        plt.subplot(4,4,k+1)
        plt.plot(self.mesh.cell_center_positions, self.simulation.breakdown_estimators_refinement[:,0])
        plt.title('Increase estimator')

        plt.subplot(4,4,k+2)
        plt.plot(self.mesh.cell_center_positions,data_array[:,1])
        plt.scatter(self.mesh.cell_center_positions,(data_array[:,-1]*np.max(data_array[:,1])\
                +(order-data_array[:,-1])*np.min(data_array[:,1]))/order,s=order,color = 'hotpink')
        plt.title('orders vs density')

        plt.show()
