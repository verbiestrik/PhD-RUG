from abc import ABC, abstractmethod
import numpy as np
import PDE
import Mesh
import SpatialDiscretization
import TimeIntegration

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
    boundary_condition: str
        the used boundary condition
    initial_condition: str
        the initial condition for the simulation
    spatial_discretization: spatial_discretization
        the numerical method for the spatial discretization

    
    Abstract methods
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

    @abstractmethod
    def __init__(self):
        """
        Implemented in the child classes.
        """
        pass

    @abstractmethod
    def run_simulation(self,
                       t_end: float) -> np.array:
        """
        Runs the simulation until the end time t_end and returns the end values of the variables

        Parameters
        ----------
        t_end : float
            end time of the simulation
        
        Returns
        -------
        values: numpy arrays
            data array containing values of the variables at the end of the simulation

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
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization,
                 time_integration: TimeIntegration.TimeIntegration):
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
        self.spatial_discretization = spatial_discretization
        self.time_integration = time_integration

    def run_simulation(self,
                       t_end: float,
                       **kwargs) -> np.array:

        g = kwargs["g"] if "g" in kwargs else 1
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)

        CFL = 0.7
        t = 0

        def system_matrix(cell_values, **kwargs):
            return self.pde_type.compute_system_matrix(self.order,cell_values, **kwargs)

        def source_term(cell_values, delta_t, **kwargs):
            return self.pde_type.compute_source_term(self.order,cell_values,delta_t,**kwargs)

        step = 0

        while t < t_end:

            # update boundary conditions
            values[0,:] = self._update_boundary_conditions(values, 'left')
            values[self.mesh.resolution+1,:] = self._update_boundary_conditions(values,'right')
            
            wave_speed_sqrt = values[:,0]*g
            for i in range(self.order):
                wave_speed_sqrt += np.divide(values[:,i+2]*values[:,i+2],values[:,0]*values[:,0])
            max_wave_speed_plus = np.max(np.abs(np.divide(values[:,1],values[:,0]) + np.sqrt(wave_speed_sqrt)))
            max_wave_speed_min = np.max(np.abs(np.divide(values[:,1],values[:,0]) - np.sqrt(wave_speed_sqrt)))
            max_speed = max(max_wave_speed_plus,max_wave_speed_min)

            delta_t = CFL*delta_x/max_speed #TODO implement CFL condition'
            #delta_t = 0.0005

            previous_values = np.copy(values)

            for i in range(1,self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    previous_values[i-1,:],
                    previous_values[i,:],
                    system_matrix,
                    'positive',
                    delta_t,
                    delta_x,
                    **kwargs) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    previous_values[i,:],
                    previous_values[i+1,:],
                    system_matrix,
                    'negative',
                    delta_t,
                    delta_x,
                    **kwargs) 
                values[i,:] = previous_values[i,:] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
                values[i,:] = self.time_integration.integrate(values[i,:],source_term,delta_t,**kwargs)
            print(t)
            t += delta_t
            step += 1
        simulation_data = self._post_processing(values)
        return simulation_data

    def _get_initial_conditions(self,
                               cell_centers_x: np.array) -> np.array:

        """
        construct the initial values for the variables

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
            initial_values[i+1,:] = self.pde_type.get_initial_values(self.order,self.initial_condition,cell_centers_x[i])            
        
        return initial_values
    
    def _update_boundary_conditions(self,
                                    values: np.array,
                                    boundary: str) -> np.array:
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
                         values) -> np.array:

        data_array = np.zeros((self.mesh.resolution,self.number_of_variables+1)) # rewrite this such that it can be generalized to other PDE models

        data_array[:,0] = self.mesh.cell_center_positions
        data_array[:,1] = values[1:-1,0]
        data_array[:,2] = np.divide(values[1:-1,1],data_array[:,1])
        for j in range(self.order): #TODO: this is unnecessary routine here
            data_array[:,j+3] = np.divide(values[1:-1,j+2],data_array[:,1])

        return data_array
    

class SpatiallyAdaptiveSimulation1D(Simulation):

    """
    This abstract class represents a spatially adaptive simulation in 1D.

    ...

    Attributes
    ----------
    boundary_interfaces: list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders: list of integers
        list of the order of the moment model in each subdomain
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition: str
        the used boundary condition
    initial_condition: str
        the initial condition for the simulation
    breakdown_criterion: str
        breadown criterion for domain decomposition
    spatial_discretization: spatial_discretization
        the numerical method for the spatial discretization

    Inherited abstract methods from interface Simulation:
    -------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values

    
    Implemented methods from interface Simulation
    -------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting

    
    Abstract methods
    ----------------
    def _resontruct_subdomains(values):
        find boundary interfaces from the cellwise model orders    
    

    Instance methods
    ----------------
    def _calculate_boundary_interfaces(self)
        calculates the indices corresponding to the boundary interfaces
    def _update_domain_decomposition_pointwise(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions in each point of the domain
    """

    def __init__(self,
                 start_order: list,
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization,
                 time_integration: TimeIntegration.TimeIntegration):

        """
        Constructs all the necessary attributes for the SpatiallyAdaptiveSimulation1D object.

        Parameters
        ----------
        start_order: integer
            list of the orders of the moment model in each subdomain
        pde_type : str
            the partial differential equations that is simulated
        numbers_of_variables
            list of the number of state variables in each subdomain
        mesh : RectangularMesh
            the used mesh
        boundary_condition: str
            the used boundary condition
        initial_condition: str
            the initial condition for the simulation
        breakdown_criterion: str
            breadown criterion for domain decomposition
        spatial_discretization: SpatialDiscretization
            the numerical method for the spatial discretization
        time_integration: TimeIntegration
            the time integration method for the integration of the source term

        """

        self.orders = [start_order,start_order]
        self.pde_type = pde_type
        self.numbers_of_variables = [pde_type.compute_number_of_variables(start_order),pde_type.compute_number_of_variables(start_order)]
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.breakdown_criterion = breakdown_criterion
        self.spatial_discretization = spatial_discretization
        self.time_integration = time_integration

        self.boundary_interfaces_discretized = [np.floor_divide(self.mesh.resolution,2)]
        self.max_order = 5
        self.max_number_of_variables = pde_type.compute_number_of_variables(self.max_order)

        self.orders_cellwise = []
        self.numbers_of_variables_cellwise = []
        for i in range(self.mesh.resolution+2):
            self.orders_cellwise.append(start_order)
            self.numbers_of_variables_cellwise.append(self.pde_type.compute_number_of_variables(start_order))

        self.dom_decomp_val_res1 = np.zeros(self.mesh.resolution)
        self.dom_decomp_val_res2 = np.zeros(self.mesh.resolution)

        self.breakdown_estimators = np.zeros((self.mesh.resolution,self.max_order+4))

    # This method is not used anymore in the current version
    def _calculate_boundary_interfaces(self):

        """
        Calculates the interfaces in the grid that correspond to the physical positions of the boundary interfaces

        Parameters
        ----------
        None
        
        
        Returns
        -------
        None
        """

        self.boundary_interfaces_discretized.append(round((self.boundary_interfaces[0] - self.mesh.boundaries[0])/(self.mesh.boundaries[1] - self.mesh.boundaries[0])*self.mesh.resolution))
        for i in range(1,len(self.boundary_interfaces)):
            self.boundary_interfaces_discretized.append(self.boundary_interfaces_discretized[i-1]
                + round((self.boundary_interfaces[i]-self.boundary_interfaces[i-1])/(self.mesh.boundaries[1]-self.mesh.boundaries[0])*self.mesh.resolution))
              
    def _update_boundary_conditions(self,
                                    values: np.array,
                                    boundary: str) -> np.array:

        """
        update the boundary conditions

        Parameters
        ----------
        values_boundary : numpy 2D array 
            the values of the variables in each mesh cell
        boundary : str
            the boundary at which we want to prescribe a boundary condition
        
        Returns
        -------
        values_ghost: numpy 1D array
            the values of the variables in the ghost cell

        """
        if self.boundary_condition == 'INFLOW_OUTFLOW':
            if boundary == 'left':
                values_ghost = values[1,:self.numbers_of_variables[0]]
            else:
                values_ghost = values[-2,:self.numbers_of_variables[-1]]
        elif self.boundary_condition == 'PERIODIC':
            if boundary == 'left':
                values_ghost = values[-2,:self.numbers_of_variables[-1]]
            else:
                values_ghost = values[1,:self.numbers_of_variables[0]]

        return values_ghost 

    def _get_initial_conditions(self,
                                cell_centers_x: np.array) -> np.array:
            
        """
        construct the initial values for the variables

        Parameters
        ----------
        cell_centers_x : numpy 1D array
            the centers of the cells
        
        Returns
        -------
        initial_values: list of numpy arrays
            initial values of the variables in each grid cell

        """

        initial_values = np.zeros((self.mesh.resolution+2,self.max_number_of_variables))      

        right_boundary_subdomain = 0
        for m in range(len(self.boundary_interfaces_discretized)):
            left_boundary_subdomain = right_boundary_subdomain
            if self.orders[m+1] > self.orders[m]: 
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]-2
            else:
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]+2
            for i in range(left_boundary_subdomain,right_boundary_subdomain):
                initial_values[i+1,:2+self.orders[m]] = self.pde_type.get_initial_values(self.orders[m],self.initial_condition,cell_centers_x[i])  
        for i in range(right_boundary_subdomain,self.mesh.resolution):
            initial_values[i+1,:2+self.orders[-1]] = self.pde_type.get_initial_values(self.orders[-1],self.initial_condition,cell_centers_x[i])
            
        return initial_values

    def _update_domain_decomposition_pointwise(self,
                                               values: np.array,
                                               delta_x: float,
                                               delta_t: float):
        
        """
        updates the domain decomposition in each point of the domain

        Parameters
        ----------
        values : numpy 2D array
            the values of the variables in each mesh cell
        tolerance_up : float
            the threshold value for increasing the order
        tolerance_down : float
            the threshold value for decreasing the order
        
        Returns
        -------
        None

        """

        # breakdown_criteria_flags = self.pde_type.compute_all_breakdown_criteria(values,
        #                                                                         self.orders_cellwise,
        #                                                                         self.numbers_of_variables_cellwise,
        #                                                                         self.mesh.resolution,
        #                                                                         delta_x)

        breakdown_criteria_flags = self._compute_breakdown_criteria_full(values,
                                                                         self.mesh.resolution,
                                                                         delta_x,
                                                                         delta_t)

        for i in range(self.mesh.resolution):
            if self.orders_cellwise[i+1] < self.max_order and breakdown_criteria_flags[i] == 1: 
                self.orders_cellwise[i+1] += 1
                self.numbers_of_variables_cellwise[i+1] += 1
            elif self.orders_cellwise[i+1] > 0 and breakdown_criteria_flags[i] == -1:
                self.orders_cellwise[i+1] -= 1
                self.numbers_of_variables_cellwise[i+1] -= 1        
        self.orders_cellwise[0] = self.orders_cellwise[1]
        self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]
        self.orders_cellwise[-1] = self.orders_cellwise[-2]
        self.numbers_of_variables_cellwise[-1] = self.numbers_of_variables_cellwise[-2]

        # breakdown_criteria = self.pde_type.compute_breakdown_criterion(values, self.numbers_of_variables_cellwise,self.breakdown_criterion,self.mesh.resolution)
        # for i in range(self.mesh.resolution):
        #     if self.orders_cellwise[i+1] < self.max_order and breakdown_criteria[i] > tolerance_up: 
        #         self.orders_cellwise[i+1] += 1
        #         self.numbers_of_variables_cellwise[i+1] += 1
        #     elif self.orders_cellwise[i+1] > 0 and breakdown_criteria[i] < tolerance_down:
        #         self.orders_cellwise[i+1] -= 1
        #         self.numbers_of_variables_cellwise[i+1] -= 1

    def _compute_breakdown_criteria_full(self,
                                   values: np.array,
                                   n,
                                   delta_x,
                                   delta_t,
                                   **kwargs) -> np.array:        
        
        tolerance_up_source              = kwargs["tolerance_up_source"]              if "tolerance_up_source"              in kwargs else 0.5
        tolerance_up_height_gradient     = kwargs["tolerance_up_height_gradient"]     if "tolerance_up_height_gradient"     in kwargs else 0.2
        tolerance_up_momentum_gradient   = kwargs["tolerance_up_momentum_gradient"]   if "tolerance_up_momentum_gradient"   in kwargs else 0.2
        tolerance_up_moment_gradient     = kwargs["tolerance_up_moment_gradient"]     if "tolerance_up_moment_gradient"     in kwargs else 0.2
        tolerance_down_height_gradient   = kwargs["tolerance_down_height_gradient"]   if "tolerance_down_height_gradient"   in kwargs else 0.001
        tolerance_down_momentum_gradient = kwargs["tolerance_down_momentum_gradient"] if "tolerance_down_momentum_gradient" in kwargs else 0.001
        tolerance_down_moment_gradient   = kwargs["tolerance_down_moment_gradient"]   if "tolerance_down_moment_gradient"   in kwargs else 0.001
        tolerance_down_last_moment       = kwargs["tolerance_down_last_moment"]       if "tolerance_down_last_moment"       in kwargs else 0.001
        tolerance_down_res1              = kwargs["tolerance_down_res1"]              if "tolerance_down_res1"              in kwargs else 0.005
        tolerance_down_res2              = kwargs["tolerance_down_res2"]              if "tolerance_down_res2"              in kwargs else 0.25


        """
        TODO

        Parameters
        ----------
        values : list of numpy 1D arrays
            the values of the variables in each mesh cell
        orders : list of integers
            the order in each cell
        number_of_variables : list of integers
            the number of variables in each cell
        
        Returns
        -------
        relative_value_last_moment: numpy 2D array
            modelling complexity quantities in each mesh cell

        """
        max_order = max(self.orders)

        breakdown_criterion_flags = np.zeros(n)

        for i in range(n):
            self.breakdown_estimators[i,0] = 1*np.abs(self.pde_type.compute_source_term_lastentry(self.orders_cellwise[i+1],values[i+1,:self.numbers_of_variables_cellwise[i+1]],True,**kwargs))
            self.breakdown_estimators[i,1] = np.abs(values[i+1,self.numbers_of_variables_cellwise[i+1]-1]/values[i+1,0])
            self.breakdown_estimators[i,2] = 1*np.abs((values[i+2,0] - values[i,0]))/(2*delta_x)
            self.breakdown_estimators[i,3] = 1*np.abs((values[i+2,1] - values[i,1]))/(2*delta_x)
            for j in range(self.orders_cellwise[i+1]):
                self.breakdown_estimators[i,4+j] = 1*np.abs((values[i+2,2+j])-values[i,2+j])/(2*delta_x)
        self.breakdown_estimators[0,0] = 1*np.abs(self.pde_type.compute_source_term_lastentry(self.orders_cellwise[1],values[1,:self.numbers_of_variables_cellwise[1]],True,**kwargs))
        self.breakdown_estimators[0,1] = np.abs(values[1,self.numbers_of_variables_cellwise[1]-1])
        self.breakdown_estimators[0,2] = 1*np.abs((values[2,0] - values[1,0]))/delta_x
        self.breakdown_estimators[0,3] = 1*np.abs((values[2,1] - values[1,1]))/delta_x
        for j in range(self.orders_cellwise[i+1]):
            self.breakdown_estimators[0,4+j] = 1*np.abs((values[2,2+j])-values[2,2+j])/delta_x    

        for i in range(n):
            if self.breakdown_estimators[i,0] > tolerance_up_source:
                breakdown_criterion_flags[i] = 1
            else:
                if self.breakdown_estimators[i,2] > tolerance_up_height_gradient or self.breakdown_estimators[i,3] > tolerance_up_momentum_gradient:
                    breakdown_criterion_flags[i] = 1
                else:
                    for j in range(4,self.orders_cellwise[i+1]+4):
                        if self.breakdown_estimators[i,j] > tolerance_up_moment_gradient:
                            breakdown_criterion_flags[i] = 1
                            break
            if (breakdown_criterion_flags[i] !=1 and\
                  self.dom_decomp_val_res1[i] < tolerance_down_res1\
                    and self.dom_decomp_val_res2[i] < tolerance_down_res2\
                        and self.breakdown_estimators[i,1] < tolerance_down_last_moment):
                breakdown_criterion_flags[i] = -1
                # if self.breakdown_estimators[i,2] > tolerance_down_height_gradient or self.breakdown_estimators[i,3] > tolerance_down_momentum_gradient:
                #     breakdown_criterion_flags[i] = 0
                # else:
                #     for j in range(4,self.orders_cellwise[i+1]+4):
                #         if self.breakdown_estimators[i,j] > tolerance_down_moment_gradient:
                #             breakdown_criterion_flags[i] = 0
                #             break
        
        return breakdown_criterion_flags

    def _resontruct_subdomains(self,
                               values: np.array,
                               delta_x: float,
                               delta_t: float) -> np.array:
        
        """
        updates the domain decomposition in each point of the domain

        Parameters
        ----------
        values : numpy 2D array
            the values of the variables in each mesh cell
        tolerance_up : float
            the threshold value for increasing the order
        tolerance_down : float
            the threshold value for decreasing the order
        
        Returns
        -------
        values : numpy 2D array
            the values, updated (moments set to zero) where necessary

        """

        pass
    
    def _post_processing(self,
                         values: list) -> np.array:

        data_array = np.zeros((self.mesh.resolution,self.max_number_of_variables+2)) # rewrite this such that it can be generalized to other PDE models

        for i in range(self.mesh.resolution):
            data_array[i,0] = self.mesh.cell_center_positions[i]
        data_array[:,1] = values[1:-1,0]
        data_array[:,2] = np.divide(values[1:-1,1],data_array[:,1])
        for j in range(self.max_order): 
            data_array[:,j+3] = np.divide(values[1:-1,j+2],data_array[:,1])
        data_array[:,-1] = self.orders_cellwise[1:-1]

        print("Orders at the end of the simulation:",self.orders)
        
        return data_array
    
class NonConservativeAdaptiveSimulation1D(SpatiallyAdaptiveSimulation1D):

    """
    This class represents an adaptive simulation in 1D that can not be written in conservative form.

    ...

    Attributes
    ----------
    boundary_interfaces: list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders: list of integers
        list of the order of the moment model in each subdomain
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition: str
        the used boundary condition
    initial_condition: str
        the initial condition for the simulation
    breakdown_criterion: str
        breadown criterion for domain decomposition
    spatial_discretization: spatial_discretization
        the numerical method for the spatial discretization

    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _calculate_boundary_interfaces(self)
        calculates the indices corresponding to the boundary interfaces
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

    Implemented methods from interface Simulation
    -------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values


    Instance methods
    ----------------
    None
    """

    def __init__(self,
                 start_order: list,
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization,
                 time_integration: TimeIntegration.TimeIntegration):

        super().__init__(start_order,
                         pde_type,
                         mesh,
                         boundary_condition,
                         initial_condition,
                         breakdown_criterion,
                         spatial_discretization,
                         time_integration)

    
    def run_simulation(self,
                       t_end: float,
                       **kwargs) -> np.array:
        
        g = kwargs["g"] if "g" in kwargs else 1
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        print("Orders at the beginning of the simulation:",self.orders)

        values = self._get_initial_conditions(self.mesh.cell_center_positions)

        # # First domain decomposition
        # values = self._resontruct_subdomains(values,delta_x,delta_t)

        CFL = 0.7
        
        step_count = 0
        t = 0

        while t < t_end:
            # update boundary conditions
            values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values,'right')

            wave_speed_sqrt = values[:,0]*g
            for i in range(self.max_order):
                wave_speed_sqrt += np.divide(values[:,i+2]*values[:,i+2],values[:,0]*values[:,0])
            max_wave_speed_plus = np.max(np.abs(np.divide(values[:,1],values[:,0])+np.sqrt(wave_speed_sqrt)))
            max_wave_speed_min = np.max(np.abs(np.divide(values[:,1],values[:,0])-np.sqrt(wave_speed_sqrt)))
            max_speed = max(max_wave_speed_plus,max_wave_speed_min)
            #TODO: add method to PDE class that computes the wave speed (approximately)
      
            delta_t = CFL*delta_x/max_speed 

            # if step_count%10 == 0:
            #     values = self._resontruct_subdomains(values,delta_x)
            values = self._resontruct_subdomains(values,delta_x,delta_t)

            previous_values = np.copy(values)

            right_boundary_subdomain = 0

            for m in range(len(self.boundary_interfaces_discretized)):
                order_left = self.orders[m]
                n_variables_left = self.numbers_of_variables[m]
                order_right = self.orders[m+1]
                n_variables_right = self.numbers_of_variables[m+1]

                def system_matrix_left(cell_values,**kwargs):
                    return self.pde_type.compute_system_matrix(order_left,cell_values,**kwargs)

                def source_term_left(cell_values,delta_t,**kwargs):
                    return self.pde_type.compute_source_term(order_left,cell_values,delta_t,**kwargs)

                def system_matrix_right(cell_values,**kwargs):
                    return self.pde_type.compute_system_matrix(order_right,cell_values,**kwargs)

                def source_term_right(cell_values,delta_t,**kwargs):
                    return self.pde_type.compute_source_term(order_right,cell_values,delta_t,**kwargs)

                left_boundary_subdomain = right_boundary_subdomain+1
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]
                
                if order_right > order_left:
                    previous_values[right_boundary_subdomain-1,n_variables_left:n_variables_right] = \
                        previous_values[right_boundary_subdomain,n_variables_left:n_variables_right] # update boundary interface boundary condition 
                    for i in range(left_boundary_subdomain,right_boundary_subdomain-2):
                        generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            previous_values[i-1,:n_variables_left],
                            previous_values[i,:n_variables_left],
                            system_matrix_left,
                            'positive',
                            delta_t,
                            delta_x,
                            **kwargs)
                        generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            previous_values[i,:n_variables_left],
                            previous_values[i+1,:n_variables_left],
                            system_matrix_left,
                            'negative',
                            delta_t,
                            delta_x,
                            **kwargs)
                        fluctuation_plus = generalized_roe_plus@(previous_values[i,:n_variables_left]-previous_values[i-1,:n_variables_left])
                        fluctuation_minus = generalized_roe_minus@(previous_values[i+1,:n_variables_left]-previous_values[i,:n_variables_left]) 
                        values[i,:n_variables_left] = previous_values[i,:n_variables_left] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus) 
                        values[i,:n_variables_left] = self.time_integration.integrate(values[i,:n_variables_left],source_term_left,delta_t,**kwargs)
                        self.dom_decomp_val_res1[i-1] \
                            = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[i,n_variables_left-1]-previous_values[i-1,n_variables_left-1])\
                            + generalized_roe_minus[:-1,-1]*(previous_values[i+1,n_variables_left-1]-previous_values[i,n_variables_left-1]),np.inf)
                        self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_left-1]-previous_values[i,n_variables_left-1])
                    
                    # Evolution equation for the cell with index right_boundary_subdomain-2
                    generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain-3,:n_variables_left],
                        previous_values[right_boundary_subdomain-2,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain-2,:n_variables_left],
                        previous_values[right_boundary_subdomain-1,:n_variables_left],
                        system_matrix_left,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs)
                    fluctuation_plus = generalized_roe_plus@(previous_values[right_boundary_subdomain-2,:n_variables_left]\
                                                             -previous_values[right_boundary_subdomain-3,:n_variables_left])
                    fluctuation_minus = generalized_roe_minus@(previous_values[right_boundary_subdomain-1,:n_variables_left]\
                                                               -previous_values[right_boundary_subdomain-2,:n_variables_left]) 
                    values[right_boundary_subdomain-2,:n_variables_left] = (previous_values[right_boundary_subdomain-2,:n_variables_left]\
                        -delta_t/delta_x*(fluctuation_plus+fluctuation_minus)) 
                    values[right_boundary_subdomain-2,:n_variables_left] \
                        = self.time_integration.integrate(values[right_boundary_subdomain-2,:n_variables_left],source_term_left,delta_t, **kwargs)
                    self.dom_decomp_val_res1[right_boundary_subdomain-3] \
                        = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[right_boundary_subdomain-2,n_variables_left-1]-previous_values[right_boundary_subdomain-3,n_variables_left-1])\
                        + generalized_roe_minus[:-1,-1]*(previous_values[right_boundary_subdomain-1,n_variables_left-1]-previous_values[right_boundary_subdomain-2,n_variables_left-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain-3] = np.abs(values[right_boundary_subdomain-2,n_variables_left-1]-previous_values[right_boundary_subdomain-2,n_variables_left-1])                   
                    # Evolution equation for the cell with index right_boundary_subdomain-1
                    generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain-2,:n_variables_left],
                        previous_values[right_boundary_subdomain-1,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain-1,:n_variables_left],
                        previous_values[right_boundary_subdomain,:n_variables_left],
                        system_matrix_left,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs)
                    fluctuation_plus = generalized_roe_plus@(previous_values[right_boundary_subdomain-1,:n_variables_left]\
                                                             -previous_values[right_boundary_subdomain-2,:n_variables_left])
                    fluctuation_minus = generalized_roe_minus@(previous_values[right_boundary_subdomain,:n_variables_left]\
                                                               -previous_values[right_boundary_subdomain-1,:n_variables_left]) 
                    values[right_boundary_subdomain-1,:n_variables_left] = (previous_values[right_boundary_subdomain-1,:n_variables_left]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus))
                    values[right_boundary_subdomain-1,:n_variables_left]\
                        = self.time_integration.integrate(values[right_boundary_subdomain-1,:n_variables_left],source_term_left,delta_t,**kwargs)
                    self.dom_decomp_val_res1[right_boundary_subdomain-2] \
                        = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[right_boundary_subdomain-1,n_variables_left-1]-previous_values[right_boundary_subdomain-2,n_variables_left-1])\
                        + generalized_roe_minus[:-1,-1]*(previous_values[right_boundary_subdomain,n_variables_left-1]-previous_values[right_boundary_subdomain-1,n_variables_left-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain-2] = np.abs(values[right_boundary_subdomain-1,n_variables_left-1]-previous_values[right_boundary_subdomain-1,n_variables_left-1]) 


                    # Evolution equation for the cell with index right_boundary_subdomain
                    generalized_roe_plus_Full = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain-1,:n_variables_right],
                        previous_values[right_boundary_subdomain,:n_variables_right],
                        system_matrix_right,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    generalized_roe_plus_Restricted = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain-1,:n_variables_left],
                        previous_values[right_boundary_subdomain,:n_variables_left],
                        system_matrix_left,'positive',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain,:n_variables_right],
                        previous_values[right_boundary_subdomain+1,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs)  
                    fluctuation_plus_Full = generalized_roe_plus_Full@(previous_values[right_boundary_subdomain,:n_variables_right]\
                                                                       -previous_values[right_boundary_subdomain-1,:n_variables_right])
                    fluctuation_plus_Restricted = generalized_roe_plus_Restricted@(previous_values[right_boundary_subdomain,:n_variables_left]\
                                                                                   -previous_values[right_boundary_subdomain-1,:n_variables_left])
                    fluctuation_minus = generalized_roe_minus@(previous_values[right_boundary_subdomain+1,:n_variables_right]\
                                                               -previous_values[right_boundary_subdomain,:n_variables_right])
                    values[right_boundary_subdomain,:n_variables_left] = (previous_values[right_boundary_subdomain,:n_variables_left]
                    -delta_t/delta_x*(fluctuation_plus_Restricted+fluctuation_minus[:n_variables_left])) # solve FVM equations for first moments
                    values[right_boundary_subdomain,n_variables_left:n_variables_right] = (previous_values[right_boundary_subdomain,n_variables_left:n_variables_right]
                    -delta_t/delta_x*(fluctuation_plus_Full[n_variables_left:n_variables_right]+fluctuation_minus[n_variables_left:n_variables_right])) # solve FVM equations for last moment
                    values[right_boundary_subdomain,:n_variables_right]\
                        = self.time_integration.integrate(values[right_boundary_subdomain,:n_variables_right],source_term_right,delta_t,**kwargs)
                    
                    self.dom_decomp_val_res1[right_boundary_subdomain-1] \
                        = np.linalg.norm(generalized_roe_minus[:-1,-1]*(previous_values[right_boundary_subdomain+1,n_variables_right-1]-previous_values[right_boundary_subdomain,n_variables_right-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain-1] = np.abs(values[right_boundary_subdomain,n_variables_right-1]-previous_values[right_boundary_subdomain,n_variables_right-1]) 
                    
            
                else:
                    previous_values[right_boundary_subdomain+2,n_variables_right:n_variables_left] = \
                        previous_values[right_boundary_subdomain+1,n_variables_right:n_variables_left] # update boundary interface boundary condition
                    for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                        generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            previous_values[i-1,:n_variables_left],
                            previous_values[i,:n_variables_left],
                            system_matrix_left,
                            'positive',
                            delta_t,
                            delta_x,
                            **kwargs) 
                        generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            previous_values[i,:n_variables_left],
                            previous_values[i+1,:n_variables_left],
                            system_matrix_left,
                            'negative',
                            delta_t,
                            delta_x,
                            **kwargs)  
                        fluctuation_plus = generalized_roe_plus@(previous_values[i,:n_variables_left]-previous_values[i-1,:n_variables_left])
                        fluctuation_minus = generalized_roe_minus@(previous_values[i+1,:n_variables_left]-previous_values[i,:n_variables_left])
                        values[i,:n_variables_left] = previous_values[i,:n_variables_left]-delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
                        values[i,:n_variables_left] = self.time_integration.integrate(values[i,:n_variables_left],source_term_left,delta_t,**kwargs)

                        self.dom_decomp_val_res1[i-1] \
                            = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[i,n_variables_left-1]-previous_values[i-1,n_variables_left-1])\
                            + generalized_roe_minus[:-1,-1]*(previous_values[i+1,n_variables_left-1]-previous_values[i,n_variables_left-1]),np.inf)
                        self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_left-1]-previous_values[i,n_variables_left-1])
                    
                    # Evolution equation for the cell with index right_boundary_subdomain+1
                    generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain,:n_variables_left],
                        previous_values[right_boundary_subdomain+1,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs)             
                    generalized_roe_minus_Full = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain+1,:n_variables_left],
                        previous_values[right_boundary_subdomain+2,:n_variables_left],
                        system_matrix_left,'negative',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    generalized_roe_minus_Restricted = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain+1,:n_variables_right],
                        previous_values[right_boundary_subdomain+2,:n_variables_right],
                        system_matrix_right,'negative',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    fluctuation_plus = generalized_roe_plus@(previous_values[right_boundary_subdomain+1,:n_variables_left]\
                                                             -previous_values[right_boundary_subdomain,:n_variables_left])
                    fluctuation_minus_Full = generalized_roe_minus_Full@(previous_values[right_boundary_subdomain+2,:n_variables_left]\
                                                             -previous_values[right_boundary_subdomain+1,:n_variables_left])
                    fluctuation_minus_Restricted = generalized_roe_minus_Restricted@(previous_values[right_boundary_subdomain+2,:n_variables_right]\
                                                             -previous_values[right_boundary_subdomain+1,:n_variables_right])
                    

                    values[right_boundary_subdomain+1,:n_variables_right] = (previous_values[right_boundary_subdomain+1,:n_variables_right] 
                    -delta_t/delta_x*(fluctuation_plus[:n_variables_right]+fluctuation_minus_Restricted)) # solve FVM equations for first moments
                    values[right_boundary_subdomain+1,n_variables_right:n_variables_left] = (previous_values[right_boundary_subdomain+1,n_variables_right:n_variables_left] 
                    -delta_t/delta_x*(fluctuation_plus[n_variables_right:n_variables_left]+fluctuation_minus_Full[n_variables_right:n_variables_left])) # solve FVM equations for last moments
                    values[right_boundary_subdomain+1,:n_variables_left] \
                        = self.time_integration.integrate(values[right_boundary_subdomain+1,:n_variables_left],source_term_left,delta_t,**kwargs)
                    
                    self.dom_decomp_val_res1[right_boundary_subdomain] \
                        = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[right_boundary_subdomain+1,n_variables_left-1]-previous_values[right_boundary_subdomain,n_variables_left-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain] = np.abs(values[right_boundary_subdomain+1,n_variables_left-1]-previous_values[right_boundary_subdomain+1,n_variables_left-1])

                    # Evolution equation for the cell with index right_boundary_subdomain+2
                    generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain+1,:n_variables_right],
                        previous_values[right_boundary_subdomain+2,:n_variables_right],
                        system_matrix_right,'positive',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain+2,:n_variables_right],
                        previous_values[right_boundary_subdomain+3,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    fluctuation_plus = generalized_roe_plus@(previous_values[right_boundary_subdomain+2,:n_variables_right]\
                                                             -previous_values[right_boundary_subdomain+1,:n_variables_right])
                    fluctuation_minus = generalized_roe_minus@(previous_values[right_boundary_subdomain+3,:n_variables_right]\
                                                              -previous_values[right_boundary_subdomain+2,:n_variables_right])
                    values[right_boundary_subdomain+2,:n_variables_right] = (previous_values[right_boundary_subdomain+2,:n_variables_right]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus)) # solve FVM equations
                    values[right_boundary_subdomain+2,:n_variables_right] \
                        = self.time_integration.integrate(values[right_boundary_subdomain+2,:n_variables_right],source_term_right,delta_t,**kwargs)
                    
                    self.dom_decomp_val_res1[right_boundary_subdomain+1] \
                        = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[right_boundary_subdomain+2,n_variables_right-1]-previous_values[right_boundary_subdomain+1,n_variables_right-1])\
                        + generalized_roe_minus[:-1,-1]*(previous_values[right_boundary_subdomain+3,n_variables_right-1]-previous_values[right_boundary_subdomain+2,n_variables_right-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain+1] = np.abs(values[right_boundary_subdomain+2,n_variables_right-1]-previous_values[right_boundary_subdomain+2,n_variables_right-1])
                    
                    # Evolution equation for the cell with index right_boundary_subdomain+3
                    generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain+2,:n_variables_right],
                        previous_values[right_boundary_subdomain+3,:n_variables_right],
                        system_matrix_right,'positive',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain+3,:n_variables_right],
                        previous_values[right_boundary_subdomain+4,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    fluctuation_plus = generalized_roe_plus@(previous_values[right_boundary_subdomain+3,:n_variables_right]\
                                                             -previous_values[right_boundary_subdomain+2,:n_variables_right])
                    fluctuation_minus = generalized_roe_minus@(previous_values[right_boundary_subdomain+4,:n_variables_right]\
                                                               -previous_values[right_boundary_subdomain+3,:n_variables_right])
                    values[right_boundary_subdomain+3,:n_variables_right] = (previous_values[right_boundary_subdomain+3,:n_variables_right]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus)) # solve FVM equations
                    values[right_boundary_subdomain+3,:n_variables_right] \
                        = self.time_integration.integrate(values[right_boundary_subdomain+3,:n_variables_right],source_term_right,delta_t,**kwargs)

                    self.dom_decomp_val_res1[right_boundary_subdomain+2] \
                        = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[right_boundary_subdomain+3,n_variables_right-1]-previous_values[right_boundary_subdomain+2,n_variables_right-1])\
                        + generalized_roe_minus[:-1,-1]*(previous_values[right_boundary_subdomain+4,n_variables_right-1]-previous_values[right_boundary_subdomain+3,n_variables_right-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain+2] = np.abs(values[right_boundary_subdomain+3,n_variables_right-1]-previous_values[right_boundary_subdomain+3,n_variables_right-1])

                    right_boundary_subdomain += 3
            
            for i in range(right_boundary_subdomain+1,self.mesh.resolution+1):
                generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                    previous_values[i-1,:n_variables_right],
                    previous_values[i,:n_variables_right],
                    system_matrix_right,
                    'positive',
                    delta_t,
                    delta_x,
                    **kwargs) 
                generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                    previous_values[i,:n_variables_right],
                    previous_values[i+1,:n_variables_right],
                    system_matrix_right,
                    'negative',
                    delta_t,
                    delta_x,
                    **kwargs) 
                fluctuation_plus = generalized_roe_plus@(previous_values[i,:n_variables_right]-previous_values[i-1,:n_variables_right])
                fluctuation_minus = generalized_roe_minus@(previous_values[i+1,:n_variables_right]-previous_values[i,:n_variables_right])
                values[i,:n_variables_right] = previous_values[i,:n_variables_right] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
                values[i,:n_variables_right] = self.time_integration.integrate(values[i,:n_variables_right],source_term_right,delta_t,**kwargs)
            
                self.dom_decomp_val_res1[i-1] \
                    = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[i,n_variables_right-1]-previous_values[i-1,n_variables_right-1])\
                    + generalized_roe_minus[:-1,-1]*(previous_values[i+1,n_variables_right-1]-previous_values[i,n_variables_right-1]),np.inf)
                self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_right-1]-previous_values[i,n_variables_right-1])

            step_count += 1
            print(t)
            t+=delta_t
        
        self.dom_decomp_val_res1 = self.dom_decomp_val_res1/delta_x
        self.dom_decomp_val_res2 = self.dom_decomp_val_res2/delta_t
        simulation_data = self._post_processing(values)
        return simulation_data
    
    def _resontruct_subdomains(self,
                               values: np.array,
                               delta_x: float,
                               delta_t: float) -> np.array:

        boundary_interfaces = []
        orders_merged = []
        number_of_variables_merged = []
        orders_out = []
        number_of_variables_out = []

        super()._update_domain_decomposition_pointwise(values,delta_x,delta_t)

        for i in range(1,self.mesh.resolution-3,4):
            local_order = max(self.orders_cellwise[i:i+4])
            orders_merged.append(local_order)
            number_of_variables_merged.append(self.pde_type.compute_number_of_variables(local_order))

        orders_merged[-1] = max(self.orders_cellwise[i:-1])
        number_of_variables_merged[-1] = self.pde_type.compute_number_of_variables(orders_merged[-1])

        orders_out.append(orders_merged[0])
        number_of_variables_out.append(number_of_variables_merged[0])
        for i in range(len(orders_merged)-1):
            if orders_merged[i] > orders_merged[i+1]:
                orders_out.append(orders_merged[i+1])
                number_of_variables_out.append(number_of_variables_merged[i+1])
                boundary_interfaces.append(4*(i+1)+1)
            elif orders_merged[i] < orders_merged[i+1]:
                orders_out.append(orders_merged[i+1])
                number_of_variables_out.append(number_of_variables_merged[i+1])
                boundary_interfaces.append(4*(i+1)+1)
            # else:
            #     values[4*i+1:4*(i+1)+1,number_of_variables_merged[i]:self.max_number_of_variables] = 0

        self.orders = orders_out
        self.boundary_interfaces_discretized = boundary_interfaces
        self.numbers_of_variables = number_of_variables_out
        if len(self.boundary_interfaces_discretized) == 0:
            self.orders.append(self.orders[0])
            self.numbers_of_variables.append(int(self.numbers_of_variables[0]))
            self.boundary_interfaces_discretized.append(np.floor_divide(self.mesh.resolution,2))

        # Set undefined moments to zero
        left_boundary = 0
        for m in range(len(self.boundary_interfaces_discretized)):
            right_boundary = self.boundary_interfaces_discretized[m]
            if self.numbers_of_variables[m+1] > self.numbers_of_variables[m]:
                values[left_boundary:right_boundary,self.numbers_of_variables[m]:self.max_number_of_variables]=0
                left_boundary = right_boundary
            else:
                values[left_boundary:right_boundary+2,self.numbers_of_variables[m]:self.max_number_of_variables]=0
                left_boundary = right_boundary+2
        values[left_boundary:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables]=0 

        return values

                

class ConservativeAdaptiveSimulation1D(SpatiallyAdaptiveSimulation1D):
    """
    This class represents an adaptive simulation in 1D that can 
    be written in conservative form.

    ...

    Attributes
    ----------
    boundary_interfaces: list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders: list of integers
        list of the order of the moment model in each subdomain
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition: str
        the used boundary condition
    initial_condition: str
        the initial condition for the simulation
    breakdown_criterion: str
        breadown criterion for domain decomposition
    spatial_discretization: spatial_discretization
        the numerical method for the spatial discretization

    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    -------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _calculate_boundary_interfaces(self)
        calculates the indices corresponding to the boundary interfaces
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions


    implemented methods from interface SpatiallyAdaptiveSimulation1D
    -----------------------------------------------------------------
    def _reconstruct_subdomains(self,values)
        finds the boundary interfaces from the cellwise model orders


    Implemented methods from interface Simulation
    ---------------------------------------------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values


    Instance methods
    ----------------
    None
    """

    def __init__(self,
                 start_order: list,
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization,
                 time_integration: TimeIntegration.TimeIntegration):

        super().__init__(start_order,
                         pde_type,
                         mesh,
                         boundary_condition,
                         initial_condition,
                         breakdown_criterion,
                         spatial_discretization,
                         time_integration)


    
    def run_simulation(self,
                       t_end: float,
                       **kwargs) -> np.array:
        
        g = kwargs["g"] if "g" in kwargs else 1
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        #self._calculate_boundary_interfaces()

        print("Orders at the beginning of the simulation:",self.orders)
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)

        # # compute first domain decomposition
        # values = self._resontruct_subdomains(values,delta_x,delta_t)

        CFL = 0.7
        
        step_count = 0
        t = 0

        while t < t_end:

            # update boundary conditions
            values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values,'right')

            wave_speed_sqrt = values[:,0]*g
            for i in range(self.max_order):
                wave_speed_sqrt += np.divide(values[:,i+2]*values[:,i+2],values[:,0]*values[:,0])
            max_wave_speed_plus = np.max(np.abs(np.divide(values[:,1],values[:,0])+np.sqrt(wave_speed_sqrt)))
            max_wave_speed_min = np.max(np.abs(np.divide(values[:,1],values[:,0])-np.sqrt(wave_speed_sqrt)))
            max_speed = max(max_wave_speed_plus,max_wave_speed_min)
            #TODO: add method to PDE class that computes the wave speed (approximately)
      
            delta_t = CFL*delta_x/max_speed 

            # if step_count%10 == 0:
            #     values = self._resontruct_subdomains(values,delta_x)
            values = self._resontruct_subdomains(values,delta_x,delta_t)

            previous_values = np.copy(values)

            right_boundary_subdomain = -1

            for m in range(len(self.boundary_interfaces_discretized)):
                order_left = self.orders[m]
                n_variables_left = self.numbers_of_variables[m]
                order_right = self.orders[m+1]
                n_variables_right = self.numbers_of_variables[m+1]

                def system_matrix_left(cell_values, **kwargs):
                    return self.pde_type.compute_system_matrix(order_left,cell_values,**kwargs)

                def source_term_left(cell_values,delta_t,**kwargs):
                    return self.pde_type.compute_source_term(order_left,cell_values,delta_t,self.pde_type.linear_source,**kwargs)

                def system_matrix_right(cell_values, **kwargs):
                    return self.pde_type.compute_system_matrix(order_right,cell_values,**kwargs)

                def source_term_right(cell_values,delta_t,**kwargs):
                    return self.pde_type.compute_source_term(order_right,cell_values,delta_t,self.pde_type.linear_source,**kwargs)

                left_boundary_subdomain = right_boundary_subdomain+2
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]
                
                if order_right > order_left:
                    previous_values[right_boundary_subdomain,n_variables_left:n_variables_right] = np.zeros(n_variables_right-n_variables_left) 
                    for i in range(left_boundary_subdomain,right_boundary_subdomain):
                        generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            previous_values[i-1,:n_variables_left],
                            previous_values[i,:n_variables_left],
                            system_matrix_left,
                            'positive',
                            delta_t,
                            delta_x,
                            **kwargs)
                        generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            previous_values[i,:n_variables_left],
                            previous_values[i+1,:n_variables_left],
                            system_matrix_left,
                            'negative',
                            delta_t,
                            delta_x,
                            **kwargs)
                        fluctuation_plus = generalized_roe_plus@(previous_values[i,:n_variables_left]-previous_values[i-1,:n_variables_left])
                        fluctuation_minus = generalized_roe_minus@(previous_values[i+1,:n_variables_left]-previous_values[i,:n_variables_left]) 
                        values[i,:n_variables_left] = previous_values[i,:n_variables_left] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus) 
                        values[i,:n_variables_left] = self.time_integration.integrate(values[i,:n_variables_left],source_term_left,delta_t,**kwargs)
                        self.dom_decomp_val_res1[i-1] \
                            = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[i,n_variables_left-1]-previous_values[i-1,n_variables_left-1])\
                            + generalized_roe_minus[:-1,-1]*(previous_values[i+1,n_variables_left-1]-previous_values[i,n_variables_left-1]),np.inf)
                        self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_left-1]-previous_values[i,n_variables_left-1])

                    values_boundary_help = previous_values[right_boundary_subdomain+1,:n_variables_right]
                    values_boundary_help[n_variables_left:n_variables_right] = 0

                    previous_values[right_boundary_subdomain,n_variables_left:n_variables_right] = 0 # set padded value to zero

                    # Evolution equation for the cell with index right_boundary_subdomain
                    generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain-1,:n_variables_left],
                        previous_values[right_boundary_subdomain,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs)
                    generalized_roe_minus1 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain,:n_variables_right],
                        values_boundary_help,
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs)
                    generalized_roe_minus2 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        values_boundary_help,
                        previous_values[right_boundary_subdomain+1,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs)
                    # fluctuation_plus = generalized_roe_plus@(previous_values[right_boundary_subdomain,:n_variables_left]\
                    #                                          -previous_values[right_boundary_subdomain-1,:n_variables_left])
                    # fluctuation_minus = generalized_roe_minus1@(previous_values[right_boundary_subdomain+1,:n_variables_right]\
                    #                                            -previous_values[right_boundary_subdomain,:n_variables_right])
                    fluctuation_plus = generalized_roe_plus@(previous_values[right_boundary_subdomain,:n_variables_left]\
                                                             -previous_values[right_boundary_subdomain-1,:n_variables_left])
                    fluctuation_minus = generalized_roe_minus1@(values_boundary_help\
                                                               -previous_values[right_boundary_subdomain,:n_variables_right])\
                                        + generalized_roe_minus2@(previous_values[right_boundary_subdomain+1,:n_variables_right]\
                                                                  -values_boundary_help)

                    values[right_boundary_subdomain,:n_variables_left] = (previous_values[right_boundary_subdomain,:n_variables_left]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus[:n_variables_left]))
                    values[right_boundary_subdomain,:n_variables_left] \
                        = self.time_integration.integrate(values[right_boundary_subdomain,:n_variables_left],source_term_left,delta_t,**kwargs)

                    self.dom_decomp_val_res1[right_boundary_subdomain-1] \
                        = np.linalg.norm(generalized_roe_minus1[:-1,-1]*(values_boundary_help[n_variables_right-1]\
                                                                        -previous_values[right_boundary_subdomain,n_variables_right-1])\
                                        +generalized_roe_minus2[:-1,-1]*(previous_values[right_boundary_subdomain+1,n_variables_right-1]\
                                                                        -values_boundary_help[n_variables_right-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain-1] \
                        = np.abs(values[right_boundary_subdomain,n_variables_left-1]-previous_values[right_boundary_subdomain,n_variables_left-1])

                    # Evolution equation for the cell with index right_boundary_subdomain+1
                    generalized_roe_plus1 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain,:n_variables_right],
                        values_boundary_help,
                        system_matrix_right,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs)
                    generalized_roe_plus2 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        values_boundary_help,
                        previous_values[right_boundary_subdomain+1,:n_variables_right],
                        system_matrix_right,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs)
                    generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain+1,:n_variables_right],
                        previous_values[right_boundary_subdomain+2,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs)
                    fluctuation_plus = generalized_roe_plus1@(values_boundary_help\
                                                             -previous_values[right_boundary_subdomain,:n_variables_right])\
                                        +generalized_roe_plus2@(previous_values[right_boundary_subdomain+1,:n_variables_right]\
                                                                -values_boundary_help)
                    fluctuation_minus = generalized_roe_minus@(previous_values[right_boundary_subdomain+2,:n_variables_right]\
                                                               -previous_values[right_boundary_subdomain+1,:n_variables_right])                     
                    values[right_boundary_subdomain+1,:n_variables_right] = (previous_values[right_boundary_subdomain+1,:n_variables_right]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus))
                    values[right_boundary_subdomain+1,:n_variables_right] \
                        = self.time_integration.integrate(values[right_boundary_subdomain+1,:n_variables_right],source_term_right,delta_t,**kwargs) 
                    self.dom_decomp_val_res1[right_boundary_subdomain] \
                        = np.linalg.norm(generalized_roe_plus1[:-1,-1]*(values_boundary_help[n_variables_right-1]\
                                                             -previous_values[right_boundary_subdomain,n_variables_right-1]\
                                        +generalized_roe_plus2[:-1,-1]*(previous_values[right_boundary_subdomain+1,n_variables_right-1]\
                                                                -values_boundary_help[n_variables_right-1]))\
                        + generalized_roe_minus[:-1,-1]*(previous_values[right_boundary_subdomain+2,n_variables_right-1]-previous_values[right_boundary_subdomain+1,n_variables_right-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain] \
                        = np.abs(values[right_boundary_subdomain+1,n_variables_right-1]-previous_values[right_boundary_subdomain+1,n_variables_right-1])

                else:

                    values_boundary_help = previous_values[right_boundary_subdomain,:n_variables_left]
                    values_boundary_help[n_variables_right:n_variables_left] = 0

                    previous_values[right_boundary_subdomain+1,n_variables_right:n_variables_left] = 0 # set padded value to zero
                    for i in range(left_boundary_subdomain,right_boundary_subdomain):
                        generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            previous_values[i-1,:n_variables_left],
                            previous_values[i,:n_variables_left],
                            system_matrix_left,
                            'positive',
                            delta_t,
                            delta_x,
                            **kwargs)
                        generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            previous_values[i,:n_variables_left],
                            previous_values[i+1,:n_variables_left],
                            system_matrix_left,
                            'negative',
                            delta_t,
                            delta_x,
                            **kwargs)
                        fluctuation_plus = generalized_roe_plus@(previous_values[i,:n_variables_left]-previous_values[i-1,:n_variables_left])
                        fluctuation_minus = generalized_roe_minus@(previous_values[i+1,:n_variables_left]-previous_values[i,:n_variables_left]) 
                        values[i,:n_variables_left] = previous_values[i,:n_variables_left] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus) 
                        values[i,:n_variables_left] = self.time_integration.integrate(values[i,:n_variables_left],source_term_left,delta_t,**kwargs)
                        self.dom_decomp_val_res1[i-1] \
                            = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[i,n_variables_left-1]-previous_values[i-1,n_variables_left-1])\
                            + generalized_roe_minus[:-1,-1]*(previous_values[i+1,n_variables_left-1]-previous_values[i,n_variables_left-1]),np.inf)
                        self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_left-1]-previous_values[i,n_variables_left-1])

                    # Evolution equation for the cell with index right_boundary_subdomain
                    generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain-1,:n_variables_left],
                        previous_values[right_boundary_subdomain,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs)
                    generalized_roe_minus1 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain,:n_variables_left],
                        values_boundary_help,
                        system_matrix_left,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs)
                    generalized_roe_minus2 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        values_boundary_help,
                        previous_values[right_boundary_subdomain+1,:n_variables_left],
                        system_matrix_left,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs)
                    fluctuation_plus = generalized_roe_plus@(previous_values[right_boundary_subdomain,:n_variables_left]\
                                                             -previous_values[right_boundary_subdomain-1,:n_variables_left])
                    fluctuation_minus = generalized_roe_minus1@(values_boundary_help\
                                                               -previous_values[right_boundary_subdomain,:n_variables_left])\
                                        + generalized_roe_minus2@(previous_values[right_boundary_subdomain+1,:n_variables_left]\
                                                                  -values_boundary_help)                     
                    values[right_boundary_subdomain,:n_variables_left] = (previous_values[right_boundary_subdomain,:n_variables_left]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus)) 
                    values[right_boundary_subdomain,:n_variables_left] \
                        = self.time_integration.integrate(values[right_boundary_subdomain,:n_variables_left],source_term_left,delta_t,**kwargs)

                    self.dom_decomp_val_res1[right_boundary_subdomain-1] \
                        = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[right_boundary_subdomain,n_variables_left-1]-previous_values[right_boundary_subdomain-1,n_variables_left-1])\
                        + generalized_roe_minus1[:-1,-1]*(values_boundary_help[n_variables_left-1]\
                                                               -previous_values[right_boundary_subdomain,n_variables_left-1])\
                        + generalized_roe_minus2[:-1,-1]*(previous_values[right_boundary_subdomain+1,n_variables_left-1]\
                                                                  -values_boundary_help[n_variables_left-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain-1] \
                        = np.abs(values[right_boundary_subdomain,n_variables_left-1]-previous_values[right_boundary_subdomain,n_variables_left-1])


                    # Evolution equation for the cell with index right_boundary_subdomain+1
                    generalized_roe_plus1 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain,:n_variables_left],
                        values_boundary_help,
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs)
                    generalized_roe_plus2 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        values_boundary_help,
                        previous_values[right_boundary_subdomain+1,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x,
                        **kwargs)
                    generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        previous_values[right_boundary_subdomain+1,:n_variables_right],
                        previous_values[right_boundary_subdomain+2,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x,
                        **kwargs) 
                    fluctuation_plus = generalized_roe_plus1@(values_boundary_help\
                                                             -previous_values[right_boundary_subdomain,:n_variables_left])\
                                        + generalized_roe_plus2@(previous_values[right_boundary_subdomain+1,:n_variables_left]\
                                                                 -values_boundary_help)
                    fluctuation_minus = generalized_roe_minus@(previous_values[right_boundary_subdomain+2,:n_variables_right]\
                                                               -previous_values[right_boundary_subdomain+1,:n_variables_right])                        
                    values[right_boundary_subdomain+1,:n_variables_right] = (previous_values[right_boundary_subdomain+1,:n_variables_right]
                    -delta_t/delta_x*(fluctuation_plus[:n_variables_right]+fluctuation_minus))
                    values[right_boundary_subdomain+1,:n_variables_right] \
                        = self.time_integration.integrate(values[right_boundary_subdomain+1,:n_variables_right],source_term_right,delta_t,**kwargs) 
                    self.dom_decomp_val_res1[right_boundary_subdomain] \
                        = np.linalg.norm(generalized_roe_plus1[:-1,-1]*(values_boundary_help[n_variables_left-1]\
                                                             -previous_values[right_boundary_subdomain,n_variables_left-1])\
                                        +generalized_roe_plus2[:-1,-1]*(previous_values[right_boundary_subdomain+1,n_variables_left-1]\
                                                                 -values_boundary_help[n_variables_left-1]),np.inf)
                    self.dom_decomp_val_res2[right_boundary_subdomain] \
                        = np.abs(values[right_boundary_subdomain+1,n_variables_right-1]-previous_values[right_boundary_subdomain+1,n_variables_right-1])
                    
            for i in range(right_boundary_subdomain+1,self.mesh.resolution+1):
                generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                    previous_values[i-1,:n_variables_right],
                    previous_values[i,:n_variables_right],
                    system_matrix_right,
                    'positive',
                    delta_t,
                    delta_x,
                    **kwargs) 
                generalized_roe_minus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                    previous_values[i,:n_variables_right],
                    previous_values[i+1,:n_variables_right],
                    system_matrix_right,
                    'negative',
                    delta_t,
                    delta_x,
                    **kwargs) 
                fluctuation_plus = generalized_roe_plus@(previous_values[i,:n_variables_right]-previous_values[i-1,:n_variables_right])
                fluctuation_minus = generalized_roe_minus@(previous_values[i+1,:n_variables_right]-previous_values[i,:n_variables_right])
                values[i,:n_variables_right] = previous_values[i,:n_variables_right] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
                values[i,:n_variables_right] = self.time_integration.integrate(values[i,:n_variables_right],source_term_right,delta_t,**kwargs)
            
                self.dom_decomp_val_res1[i-1] \
                    = np.linalg.norm(generalized_roe_plus[:-1,-1]*(previous_values[i,n_variables_right-1]-previous_values[i-1,n_variables_right-1])\
                    + generalized_roe_minus[:-1,-1]*(previous_values[i+1,n_variables_right-1]-previous_values[i,n_variables_right-1]),np.inf)
                self.dom_decomp_val_res2[i-1] \
                    = np.abs(values[i,n_variables_right-1]-previous_values[i,n_variables_right-1])
            
            step_count += 1
            print(t)
            t+=delta_t

        self.dom_decomp_val_res1 = self.dom_decomp_val_res1/delta_x
        self.dom_decomp_val_res2 = self.dom_decomp_val_res2/delta_t

        values = simulation_data = self._post_processing(values)
        return simulation_data  

    def _resontruct_subdomains(self,
                               values: np.array,
                               delta_x: float,
                               delta_t: float) -> np.array:
        

        super()._update_domain_decomposition_pointwise(values,delta_x,delta_t)

        orders_out = []
        number_of_variables_out = []
        boundary_interfaces = []

        orders_out.append(int(self.orders_cellwise[0]))
        number_of_variables_out.append(self.pde_type.compute_number_of_variables(int(self.orders_cellwise[0])))
        for i in range(1,self.mesh.resolution+1):
            if self.orders_cellwise[i] != self.orders_cellwise[i+1]:
                orders_out.append(int(self.orders_cellwise[i+1]))
                number_of_variables_out.append(self.pde_type.compute_number_of_variables(int(self.orders_cellwise[i+1])))
                boundary_interfaces.append(int(i))

        # Set undefined moments and padded moments to zero
        right_boundary = -1
        for m in range(len(self.boundary_interfaces_discretized)):
            left_boundary = right_boundary+1
            right_boundary = self.boundary_interfaces_discretized[m]
            values[left_boundary:right_boundary+1,self.numbers_of_variables[m]:self.max_number_of_variables]=0
        values[right_boundary+1:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables]=0 

        self.orders = orders_out
        self.boundary_interfaces_discretized = boundary_interfaces
        self.numbers_of_variables = number_of_variables_out
        if len(self.boundary_interfaces_discretized) == 0:
            self.orders.append(self.orders[0])
            self.numbers_of_variables.append(self.numbers_of_variables[0])
            self.boundary_interfaces_discretized.append(np.floor_divide(self.mesh.resolution,2))
        
        return values

    # def _resontruct_subdomains(self,
    #                            values: np.array,
    #                            delta_x: float,
    #                            delta_t: float) -> np.array:

    #     boundary_interfaces = []
    #     orders_merged = []
    #     number_of_variables_merged = []
    #     orders_out = []
    #     number_of_variables_out = []

    #     super()._update_domain_decomposition_pointwise(values,delta_x,delta_t)

    #     for i in range(1,self.mesh.resolution-3,4):
    #         local_order = max(self.orders_cellwise[i:i+4])
    #         orders_merged.append(local_order)
    #         number_of_variables_merged.append(self.pde_type.compute_number_of_variables(local_order))

    #     orders_merged[-1] = max(self.orders_cellwise[i:-1])
    #     number_of_variables_merged[-1] = self.pde_type.compute_number_of_variables(orders_merged[-1])

    #     orders_out.append(orders_merged[0])
    #     number_of_variables_out.append(number_of_variables_merged[0])
    #     for i in range(len(orders_merged)-1):
    #         if orders_merged[i] > orders_merged[i+1]:
    #             orders_out.append(orders_merged[i+1])
    #             number_of_variables_out.append(number_of_variables_merged[i+1])
    #             boundary_interfaces.append(4*(i+1)+1)
    #         elif orders_merged[i] < orders_merged[i+1]:
    #             orders_out.append(orders_merged[i+1])
    #             number_of_variables_out.append(number_of_variables_merged[i+1])
    #             boundary_interfaces.append(4*(i+1)+1)
    #         # else:
    #         #     values[4*i+1:4*(i+1)+1,number_of_variables_merged[i]:self.max_number_of_variables] = 0

    #     self.orders = orders_out
    #     self.boundary_interfaces_discretized = boundary_interfaces
    #     self.numbers_of_variables = number_of_variables_out
    #     if len(self.boundary_interfaces_discretized) == 0:
    #         self.orders.append(self.orders[0])
    #         self.numbers_of_variables.append(int(self.numbers_of_variables[0]))
    #         self.boundary_interfaces_discretized.append(np.floor_divide(self.mesh.resolution,2))

    #     # Set undefined moments to zero
    #     left_boundary = 0
    #     for m in range(len(self.boundary_interfaces_discretized)):
    #         right_boundary = self.boundary_interfaces_discretized[m]+1
    #         values[left_boundary:right_boundary,self.numbers_of_variables[m]:self.max_number_of_variables]=0
    #         left_boundary = right_boundary+1
    #     values[left_boundary:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables]=0 

    #     return values

class Micro_macro(Simulation):
    """
    This interface represents a micro-macro simulation.

    ...

    Attributes
    ----------
    order: list
        orders of the micro and macro models
    pde_type : str
        the partial differential equations that is simulated
    mesh : RectangularMesh
        the used mesh
    boundary_condition: str
        the used boundary condition
    initial_condition: str
        the initial condition for the simulation
    spatial_discretization: spatial_discretization
        the numerical method for the spatial discretization
    
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
                 orders: list,
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization,
                 time_integration = TimeIntegration.TimeIntegration):
 
        self.micro_order = orders[0]
        self.macro_order = orders[1]
        self.pde_type = pde_type
        self.number_of_variables = pde_type.compute_number_of_variables(self.micro_order)
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.spatial_discretization = spatial_discretization
        self.time_integration = time_integration

    def run_simulation(self,
                       t_end: float,
                       **kwargs) -> np.array:

        g = kwargs["g"] if "g" in kwargs else 1
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution

        micro_moments = self._get_initial_conditions(self.mesh.cell_center_positions)
        macro_moments = np.zeros((self.mesh.resolution+2, self.macro_order+2))

        CFL = 0.7

        def micro_system_matrix(cell_values,**kwargs):
            return self.pde_type.compute_system_matrix(self.micro_order,cell_values,**kwargs)

        def micro_source_term(cell_values,**kwargs):
            return self.pde_type.compute_source_term(self.micro_order,cell_values,**kwargs)
        
        def macro_system_matrix(cell_values,**kwargs):
            return self.pde_type.compute_system_matrix(self.macro_order, cell_values,**kwargs)

        def macro_source_term(cell_values,**kwargs):
            return self.pde_type.compute_source_term(self.macro_order, cell_values,**kwargs)

        t = 0
        step = 0

        while t < t_end:

            # MICRO STEP
            micro_moments[0,:] = self._update_boundary_conditions(micro_moments,'left')
            micro_moments[-1,:] = self._update_boundary_conditions(micro_moments,'right')
            
            # Calculate step size using CFL condition
            wave_speed_sqrt = micro_moments[:,0]*g
            for i in range(self.micro_order):
                wave_speed_sqrt += np.divide(micro_moments[:,i+2]*micro_moments[:,i+2],micro_moments[:,0]*micro_moments[:,0])
            max_speed =  np.max(np.abs(np.divide(micro_moments[:,1],micro_moments[:,0]))+ np.sqrt(wave_speed_sqrt))

            micro_delta_t = CFL*delta_x/max_speed

            previous_values = np.copy(micro_moments)

            # Calculate the space derivative term
            for i in range(1,self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    previous_values[i-1,:],
                    previous_values[i,:],
                    micro_system_matrix,
                    'positive',
                    micro_delta_t,
                    delta_x,
                    **kwargs) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    previous_values[i,:],
                    previous_values[i+1,:],
                    micro_system_matrix,
                    'negative',
                    micro_delta_t,
                    delta_x,
                    **kwargs) 
                micro_moments[i,:] = previous_values[i,:] - micro_delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
                micro_moments[i,:] = self.time_integration.integrate(micro_moments[i,:],micro_source_term,micro_delta_t)

            t += micro_delta_t

            # RESTRICTION
            macro_moments = micro_moments[:, :self.macro_order+2]

            # MACRO STEP
            macro_moments[0,:] = self._update_boundary_conditions(macro_moments,'left')
            macro_moments[-1,:] = self._update_boundary_conditions(macro_moments,'right')

            # Calculate step size using CFL condition
            wave_speed_sqrt = macro_moments[:,0]*g
            for i in range(self.macro_order):
                wave_speed_sqrt += np.divide(macro_moments[:,i+2]*macro_moments[:,i+2],macro_moments[:,0]*macro_moments[:,0])
            max_speed =  np.max(np.abs(np.divide(macro_moments[:,1],macro_moments[:,0]))+np.sqrt(wave_speed_sqrt))

            macro_delta_t = CFL*delta_x/max_speed

            # Calculate the space derivative term
            previous_values = np.copy(macro_moments)

            for i in range(1,self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    previous_values[i-1,:],
                    previous_values[i,:],
                    macro_system_matrix,
                    'positive',
                    macro_delta_t,
                    delta_x,
                    **kwargs) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    previous_values[i,:],
                    previous_values[i+1,:],
                    macro_system_matrix,
                    'negative',
                    macro_delta_t,
                    delta_x,
                    **kwargs) 
                macro_moments[i,:] = previous_values[i,:] - macro_delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
                macro_moments[i,:] = self.time_integration.integrate(macro_moments[i,:],macro_source_term,macro_delta_t)

            t += macro_delta_t

            # MATCHING
            for idx in range(self.mesh.resolution+2):
                micro_moments[idx, self.macro_order+2:] = np.multiply(micro_moments[idx, self.macro_order+2:],(macro_moments[idx, 0]/micro_moments[idx, 0]))
            micro_moments[:, :self.macro_order+2] = macro_moments

            step += 1

        simulation_data = self._post_processing(micro_moments)
        return simulation_data

    def _get_initial_conditions(self,
                               cell_centers_x: np.array) -> np.array:

        """
        construct the initial values for the variables

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
            initial_values[i+1,:] = self.pde_type.get_initial_values(self.micro_order,self.initial_condition,cell_centers_x[i])            
        
        return initial_values
    
    def _update_boundary_conditions(self,
                                   values: np.array,
                                   boundary) -> np.array:
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
                         values) -> np.array:

        data_array = np.zeros((self.mesh.resolution,self.number_of_variables+1)) # rewrite this such that it can be generalized to other PDE models

        for i in range(self.mesh.resolution):
            data_array[i,0] = self.mesh.cell_center_positions[i]
        data_array[:,1] = values[1:-1,0]
        data_array[:,2] = np.divide(values[1:-1,1],data_array[:,1])
        for j in range(self.micro_order): #TODO: this is unnecessary routine here
            data_array[:,j+3] = np.divide(values[1:-1,j+2],data_array[:,1])

        return data_array


class ClassicalGalerkinSimulation1D(Simulation):

    """
    This class represents a classical (not spatially adaptive) simulation in 1D with intrusive uncertainty.

    ...

    Attributes
    ----------
    mom_order: integer
        moment order of the model
    SG_order: integer
        stochastic Galerkin order of the model
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

    Implemented methods from interface Simulation
    -------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values, boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    """

    def __init__(self,
                 mom_order: int,
                 SG_order: int,
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization,
                 time_integration: TimeIntegration.TimeIntegration):
        """
        Constructs all the necessary attributes for the ClassicalGalerkinSimulation1D object.

        Parameters
        ----------
        mom_order: integer
            moment order of the model
        SG_order: integer
            stochastic Galerkin order of the model
        pde_type : str
            the partial differential equations that is simulated
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
        self.mom_order = mom_order
        self.SG_order = SG_order
        self.pde_type = pde_type
        self.number_of_variables = pde_type.compute_number_of_variables(self.mom_order, self.SG_order)
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.spatial_discretization = spatial_discretization
        self.time_integration = time_integration

    def run_simulation(self,
                       t_end: float,
                       **kwargs) -> np.array:

        g = kwargs["g"] if "g" in kwargs else 1
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)

        CFL = 0.7
        t = 0

        def system_matrix(cell_values, **kwargs):
            return self.pde_type.compute_system_matrix(self.mom_order, self.SG_order, cell_values, **kwargs)

        def source_term(cell_values,delta_t, **kwargs):
            return self.pde_type.compute_source_term(self.mom_order, self.SG_order, cell_values, delta_t, **kwargs)

        step = 0
        
        while t < t_end:

            # update boundary conditions
            values[0,:] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:] = self._update_boundary_conditions(values,'right')
            
            denominator = values[:,0]
            if self.SG_order > 0:
                for i in range(1, self.SG_order):
                    denominator += values[:,i]
            wave_speed_part1 = np.divide(values[:,self.SG_order+1], denominator)
            wave_speed_part2 = values[:,0]*g
            if self.SG_order > 0:
                for i in range(1, self.SG_order):
                    wave_speed_part1 += np.divide(values[:,self.SG_order+1+i], denominator)
                    wave_speed_part2 += values[:,i]*g
            if self.mom_order > 0:
                wave_speed_part2 += np.divide(values[:,2*self.SG_order+2]*values[:,2*self.SG_order+2],denominator*denominator)
                if self.SG_order > 0:
                    for i in range(1, self.SG_order):
                        wave_speed_part2 += np.divide(values[:,2*self.SG_order+2+i]*values[:,2*self.SG_order+2+i],denominator*denominator)
            max_wave_speed_plus = np.max(np.abs(wave_speed_part1 + np.sqrt(wave_speed_part2)))
            max_wave_speed_min = np.max(np.abs(wave_speed_part1 - np.sqrt(wave_speed_part2)))
            max_speed = max(max_wave_speed_plus,max_wave_speed_min)

            delta_t = CFL*delta_x/max_speed
            previous_values = np.copy(values)

            for i in range(1, self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    previous_values[i-1,:],
                    previous_values[i,:],
                    system_matrix,
                    'positive',
                    delta_t,
                    delta_x,
                    **kwargs) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    previous_values[i,:],
                    previous_values[i+1,:],
                    system_matrix,
                    'negative',
                    delta_t,
                    delta_x,
                    **kwargs)
                values[i,:] = previous_values[i,:] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
                values[i,:] = self.time_integration.integrate(values[i,:],source_term,delta_t,**kwargs)
            print(t)
            t += delta_t
            step += 1
        simulation_data = self._post_processing(values)
        return simulation_data
            

    def _get_initial_conditions(self,
                               cell_centers_x: np.array) -> np.array:

        """
        construct the initial values for the variables

        Parameters
        ----------
        cell_centers_x : numpy 1D array
            the centers of the cells
        
        Returns
        -------
        initial_values: numpy 2D array
            initial values of the variables in each grid cell

        """
        
        initial_values = np.zeros((self.mesh.resolution + 2, self.number_of_variables))

        for i in range(0, self.mesh.resolution):
            initial_values[i+1,:] = self.pde_type.get_initial_values(self.mom_order, self.SG_order, self.initial_condition, cell_centers_x[i])            
        
        return initial_values
    
    def _update_boundary_conditions(self,
                                   values: np.array,
                                   boundary: str) -> np.array:
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
                         values) -> np.array:

        data_array = np.zeros((self.mesh.resolution, self.number_of_variables + 1)) # rewrite this such that it can be generalized to other PDE models
        data_array[:,0] = self.mesh.cell_center_positions
        data_array[:,1:self.number_of_variables+1] = values[1:-1,:]
        
        return data_array    