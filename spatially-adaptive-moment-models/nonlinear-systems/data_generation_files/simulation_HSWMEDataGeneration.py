from abc import ABC, abstractmethod
import numpy as np
import pde_HSWMEDataGeneration
import mesh
import spatialDiscretization
import timeIntegration
from scipy.interpolate import BarycentricInterpolator

class Simulation(ABC):

    """
    This interface represents a simulation.

    ...

    Attributes
    ----------
    pde_type : str
        the partial differential equations that is simulated
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    spatial_discretization : spatial_discretization
        the numerical method for the spatial discretization

    
    Abstract methods
    -------
    def __init__(self):
        initializes the simulation object
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    """

    @abstractmethod
    def __init__(self):
        """
        Implemented in the child classes.
        """
        pass

    @abstractmethod
    def run_simulation(self,
                       t_end: float) -> np.ndarray:
        """
        Runs the simulation until the end time t_end and returns the end values of the variables

        Parameters
        ----------
        t_end : float
            end time of the simulation
        
        Returns
        -------
        values: numpy arrays
            data array containing the positions of the grid cells and the values of the variables at the end of the simulation

        """
        pass

    @abstractmethod
    def _get_initial_conditions(self,
                               cell_centers):
        """
        Implemented and documented in the child classes. 
        """
        pass

    @abstractmethod
    def _update_boundary_conditions(self,
                                   values_boundary):
        """
        Implemented and documented in the child classes.
        """
        pass

    @abstractmethod
    def _post_processing(self,
                         end_values):
        """
        Post processes the end simulation data and prepares it for plotting

        Parameters
        ----------
        end_values : numpy array
            end values of the simulation
        
        Returns
        -------
        data_array: numpy arrays
            post processed data array containing values of the variables at the end of the simulation as 
            well as the cell center positions

        """
        pass

class ClassicalSimulation1D(Simulation):

    """
    This class represents a classical (not spatially adaptive) simulation in 1D.

    ...

    Attributes
    ----------
    order: int
        order of the moment model
    pde_type : str
        the partial differential equations that is simulated
    number_of_variables : int
        number of state variables
    mesh : RectangularMesh
        the used mesh
    boundary_condition: str
        the used boundary condition
    initial_condition: str
        the initial condition for the simulation
    spatial_discretization: SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration: TimeIntegration
        the time integration method for the right-hand side source term

    
    Implemented methods from interface Simulation
    -------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    """

    def __init__(self,
                 order: int,
                 pde_type: pde_HSWMEDataGeneration.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 initial_velocity: str,
                 initial_velocity_magnitude: str,
                 spatial_discretization: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration):
        """
        Constructs all the necessary attributes for the ClassicalSimulation1D object.

        Parameters
        ----------
        order: int
            order of the moment model
        pde_type : str
            the partial differential equations that is simulated
        number_of_variables : int
            number of state variables
        mesh : RectangularMesh
            the used mesh
        boundary_condition: str
            the used boundary condition
        initial_condition: str
            the initial condition for the simulation
        spatial_discretization: spatial_discretization
            the numerical method for the spatial discretization
        time_integration: TimeIntegration
            the time integration method for the right-hand side source term

        """
        self.order = order
        self.pde_type = pde_type
        self.number_of_variables = pde_type.compute_number_of_variables(self.order)
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.initial_velocity = initial_velocity
        self.initial_velocity_magnitude = initial_velocity_magnitude
        self.spatial_discretization = spatial_discretization
        self.time_integration = time_integration

        self.max_vals = np.zeros(8)
        self.max_grad_vals = np.zeros(8)

    def run_simulation(self,
                       t_end: float,
                       g = 1) -> tuple[np.ndarray,np.ndarray]:

        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        values = self._get_initial_conditions(self.mesh.cell_center_positions)
        fluctuations_min = np.zeros((self.mesh.resolution+1,self.number_of_variables))
        fluctuations_plus = np.zeros((self.mesh.resolution+1,self.number_of_variables))

        CFL = 0.6 #TODO: put CFL number in config file
        t = 0

        def system_matrix(cell_values):
            return self.pde_type.compute_system_matrix(self.order,cell_values)

        def source_term(cell_values,delta_t):
            return self.pde_type.compute_source_term(self.order,cell_values,delta_t)

        step = 0

        while t < t_end:

            # update boundary conditions
            values[0,:] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:] = self._update_boundary_conditions(values,'right')
            
            max_speed = self.pde_type.compute_max_wavespeed(self.order,
                                                            values)

            delta_t = CFL*delta_x/max_speed 

            for i in range(self.mesh.resolution+1):
                fluctuations_min[i,:],fluctuations_plus[i,:] = self.spatial_discretization.compute_fluctuation(
                    values[i,:],
                    values[i+1,:],
                    system_matrix,
                    delta_t,
                    delta_x)          

            for i in range(1,self.mesh.resolution+1):
                values[i,:] = values[i,:] - delta_t/delta_x*(fluctuations_plus[i-1,:]+fluctuations_min[i,:])
                values[i,:] = self.time_integration.integrate(values[i,:],source_term,delta_t)
            
            grad = (values[2:, :] - values[:-2, :])/(2*delta_x)
            max_grad_vals = np.max(np.abs(grad), axis=0)
            max_vals = np.max(np.abs(values), axis=0)

            for i in range(8):
                if self.order > i-2:
                    if max_vals[i] > self.max_vals[i]:
                        self.max_vals[i] = max_vals[i]

            for i in range(8):
                if self.order > i-2:
                    if max_grad_vals[i] > self.max_grad_vals[i]:
                        self.max_grad_vals[i] = max_grad_vals[i]

            t += delta_t
            step += 1

        simulation_data = self._post_processing(values)
        return simulation_data

    def _get_initial_conditions(self,
                               cell_centers_x: np.ndarray) -> np.ndarray:

        """
        constructs the initial values for the variables

        Parameters
        ----------
        cell_centers_x : numpy 1D array
            the centers of the cells
        
        Returns
        -------
        initial_values: numpy 2D array
            initial values of the variables in each grid cell

        """
        
        initial_values = np.zeros((self.mesh.resolution+2,self.number_of_variables))

        for i in range(0,self.mesh.resolution):
            initial_values[i+1,:] = self.pde_type.get_initial_values(self.order,
                                                                     self.initial_condition,
                                                                     self.initial_velocity,
                                                                     self.initial_velocity_magnitude,
                                                                     cell_centers_x[i])            
        
        return initial_values
    
    def _update_boundary_conditions(self,
                                   values: np.ndarray,
                                   boundary) -> np.ndarray:
        """
        update the boundary conditions

        Parameters
        ----------
        values : numpy 2D array
            values of the variables in each mesh cell
        boundary : str
            the boundary at which we are prescribing a boundary condition
        
        Returns
        -------
        values_ghost: numpy 1D array
            the values of the variables in the ghost cell

        """

        if self.boundary_condition == 'INFLOW_OUTFLOW':
            if boundary == 'left':
                values_ghost = values[1,:]
            else:
                values_ghost = values[-2,:]
        elif self.boundary_condition == 'PERIODIC':
            if boundary == 'left':
                values_ghost = values[-2,:]
            else:
                values_ghost = values[1,:]

        return values_ghost 
    
    def _post_processing(self,
                         values) -> tuple[np.ndarray,np.ndarray]:

        data_array = np.zeros((self.mesh.resolution,self.pde_type.compute_number_of_variables(self.order)+1)) # rewrite this such that it can be generalized to other PDE models
        covariates = np.zeros(19)

        for i in range(self.mesh.resolution):
            data_array[i,0] = self.mesh.cell_center_positions[i]

        data_array[:,1:] = values[1:-1,:]
        data_array[:,1:] = self.pde_type.convert_to_primitive(self.order,data_array[:,1:])

        covariates[:8] = self.max_vals
        covariates[8:16] = self.max_grad_vals
        covariates[18] = self.order

        return data_array,covariates
