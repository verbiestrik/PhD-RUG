from abc import ABC, abstractmethod
import numpy as np
import PDE
import Mesh
import SpatialDiscretization

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
    This interface represents a classical (not spatially adaptive) simulation in 1D.

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
                 order: int,
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization):
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

        """
        self.order = order
        self.pde_type = pde_type
        self.number_of_variables = pde_type.compute_number_of_variables(self.order)
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.spatial_discretization = spatial_discretization

    def run_simulation(self,
                       t_end: float,
                       **kwargs) -> np.array:

        g = kwargs["g"] if "g" in kwargs else 9.81
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)

        CFL = 0.7
        t = 0

        def system_matrix(cell_values, **kwargs):
            return self.pde_type.compute_system_matrix(self.order,cell_values, **kwargs)

        def source_term(cell_values, **kwargs):
            return self.pde_type.compute_source_term(self.order,cell_values, **kwargs)


        while t < t_end:

            # update boundary conditions
            values[0,:] = self._update_boundary_conditions(values[1,:])
            values[self.mesh.resolution+1,:] = self._update_boundary_conditions(values[self.mesh.resolution,:])
            
            wave_speed_sqrt = values[:,0]*g
            for i in range(self.order):
                wave_speed_sqrt += np.divide(values[:,i+2]*values[:,i+2],values[:,0]*values[:,0])
            max_wave_speed_plus = np.max(np.abs(np.divide(values[:,1],values[:,0])+wave_speed_sqrt))
            max_wave_speed_min = np.max(np.abs(np.divide(values[:,1],values[:,0])-wave_speed_sqrt))
            max_speed = max(max_wave_speed_plus,max_wave_speed_min)

            delta_t = CFL*delta_x/max_speed #TODO implement CFL condition

            for i in range(1,self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    values[i-1,:],
                    values[i,:],
                    system_matrix,
                    'positive',
                    delta_t,
                    delta_x) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    values[i,:],
                    values[i+1,:],
                    system_matrix,
                    'negative',
                    delta_t,
                    delta_x) 
                source_term_value = source_term(values[i,:], **kwargs) 
                values[i,:] = values[i,:] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value # solve FVM equations

            t += delta_t
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
                                   values_boundary: np.array) -> np.array:
        """
        update the boundary conditions

        Parameters
        ----------
        values_boundary : numpy 1D array #TODO: implement boundary conditions that include more cells
            the values of the variables in the boundary cell
        
        Returns
        -------
        values_ghost: numpy 1D array
            the values of the variables in the ghost cell

        """

        if self.boundary_condition == 'INFLOW_OUTFLOW':
            values_ghost = values_boundary
        
        else:
            print("This boundary condition is not implemented yet")
        
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
    This interface represents a spatially adaptive simulation in 1D.

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


    Instance methods
    ----------------
    def _calculate_boundary_interfaces(self)
        calculates the indices corresponding to the boundary interfaces
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions
    """

    def __init__(self,
                 boundary_interfaces: list,
                 orders: list,
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization):

        """
        Constructs all the necessary attributes for the SpatiallyAdaptiveSimulation1D object.

        Parameters
        ----------
        boundary_interfaces: list of floats
            list of the physical positions of the boundary interfaces between the different subdomains
        orders: list of integers
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
        spatial_discretization: spatial_discretization
            the numerical method for the spatial discretization

        """

        self.boundary_interfaces = boundary_interfaces
        self.orders = orders
        self.pde_type = pde_type
        self.numbers_of_variables = np.empty(len(self.orders),dtype=int)
        for i in range(len(self.orders)):
            self.numbers_of_variables[i] = pde_type.compute_number_of_variables(self.orders[i])
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.spatial_discretization = spatial_discretization

        self.boundary_interfaces_discretized = []
        self.max_order = 5
        self.max_number_of_variables = pde_type.compute_number_of_variables(self.max_order)

    
    def run_simulation(self,
                       t_end: float,
                       g = 9.81) -> np.array:
        
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        self._calculate_boundary_interfaces()

        print(self.orders)

        values = self._get_initial_conditions(self.mesh.cell_center_positions)

        CFL = 0.7
        
        t = 0

        while t < t_end:

            # update boundary conditions
            values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values[1,:self.numbers_of_variables[0]])
            values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values[self.mesh.resolution,:self.numbers_of_variables[-1]])

            wave_speed_sqrt = values[:,0]*g
            for i in range(self.max_order):
                wave_speed_sqrt += np.divide(values[:,i+2]*values[:,i+2],values[:,0]*values[:,0])
            max_wave_speed_plus = np.max(np.abs(np.divide(values[:,1],values[:,0])+wave_speed_sqrt))
            max_wave_speed_min = np.max(np.abs(np.divide(values[:,1],values[:,0])-wave_speed_sqrt))
            max_speed = max(max_wave_speed_plus,max_wave_speed_min)
            #TODO: add method to PDE class that computes the wave speed (approximately)
      
            delta_t = CFL*delta_x/max_speed 

            right_boundary_subdomain = 0

            for m in range(len(self.boundary_interfaces_discretized)):
                order_left = self.orders[m]
                n_variables_left = self.numbers_of_variables[m]
                order_right = self.orders[m+1]
                n_variables_right = self.numbers_of_variables[m+1]

                def system_matrix_left(cell_values):
                    return self.pde_type.compute_system_matrix(order_left,cell_values)

                def source_term_left(cell_values):
                    return self.pde_type.compute_source_term(order_left,cell_values)

                def system_matrix_right(cell_values):
                    return self.pde_type.compute_system_matrix(order_right,cell_values)

                def source_term_right(cell_values):
                    return self.pde_type.compute_source_term(order_right,cell_values)

                left_boundary_subdomain = right_boundary_subdomain+1
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]
                
                if order_right > order_left:
                    values[right_boundary_subdomain-1,n_variables_left:n_variables_right] = \
                        values[right_boundary_subdomain,n_variables_left:n_variables_right] # update boundary interface boundary condition 
                    for i in range(left_boundary_subdomain,right_boundary_subdomain-2):
                        fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                            values[i-1,:n_variables_left],
                            values[i,:n_variables_left],
                            system_matrix_left,
                            'positive',
                            delta_t,
                            delta_x) 
                        fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                            values[i,:n_variables_left],
                            values[i+1,:n_variables_left],
                            system_matrix_left,
                            'negative',
                            delta_t,
                            delta_x) 
                        source_term_value = source_term_left(values[i,:n_variables_left]) 
                        values[i,:n_variables_left] = values[i,:n_variables_left] \
                            - delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value # solve FVM equations
                    
                    # Evolution equation for the cell with index right_boundary_subdomain-2
                    fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain-3,:n_variables_left],
                        values[right_boundary_subdomain-2,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x) 
                    fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain-2,:n_variables_left],
                        values[right_boundary_subdomain-1,:n_variables_left],
                        system_matrix_left,
                        'negative',
                        delta_t,
                        delta_x) 
                    source_term_value = source_term_left(values[right_boundary_subdomain-2,:n_variables_left]) 
                    values[right_boundary_subdomain-2,:n_variables_left] = (values[right_boundary_subdomain-2,:n_variables_left]\
                        -delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value) # solve FVM equations
                    
                    # Evolution equation for the cell with index right_boundary_subdomain-1
                    fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain-2,:n_variables_left],
                        values[right_boundary_subdomain-1,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x) 
                    fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain-1,:n_variables_left],
                        values[right_boundary_subdomain,:n_variables_left],
                        system_matrix_left,
                        'negative',
                        delta_t,
                        delta_x) 
                    source_term_value = source_term_left(values[right_boundary_subdomain-1,:n_variables_left]) 
                    values[right_boundary_subdomain-1,:n_variables_left] = (values[right_boundary_subdomain-1,:n_variables_left]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value) # solve FVM equations


                    # Evolution equation for the cell with index right_boundary_subdomain
                    fluctuation_plus_Full = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain-1,:n_variables_right],
                        values[right_boundary_subdomain,:n_variables_right],
                        system_matrix_right,
                        'positive',
                        delta_t,
                        delta_x) 
                    fluctuation_plus_Restricted = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain-1,:n_variables_left],
                        values[right_boundary_subdomain,:n_variables_left],
                        system_matrix_left,'positive',
                        delta_t,
                        delta_x) 
                    fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain,:n_variables_right],
                        values[right_boundary_subdomain+1,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x) 
                    source_term_value = source_term_right(values[right_boundary_subdomain,:n_variables_right]) 
                    
                    values[right_boundary_subdomain,:n_variables_left] = (values[right_boundary_subdomain,:n_variables_left]
                    -delta_t/delta_x*(fluctuation_plus_Restricted+fluctuation_minus[:n_variables_left])
                    +delta_t*source_term_value[:n_variables_left]) # solve FVM equations for first moments
                    values[right_boundary_subdomain,n_variables_left:n_variables_right] = (values[right_boundary_subdomain,n_variables_left:n_variables_right]
                    -delta_t/delta_x*(fluctuation_plus_Full[n_variables_left:n_variables_right]+fluctuation_minus[n_variables_left:n_variables_right])
                    +delta_t*source_term_value[n_variables_left:n_variables_right]) # solve FVM equations for last moment
                else:
                    values[right_boundary_subdomain+2,n_variables_right:n_variables_left] = \
                        values[right_boundary_subdomain+1,n_variables_right:n_variables_left] # update boundary interface boundary condition
                    for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                        fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                            values[i-1,:n_variables_left],
                            values[i,:n_variables_left],
                            system_matrix_left,
                            'positive',
                            delta_t,
                            delta_x) 
                        fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                            values[i,:n_variables_left],
                            values[i+1,:n_variables_left],
                            system_matrix_left,
                            'negative',
                            delta_t,
                            delta_x) 
                        source_term_value = source_term_left(values[i,:n_variables_left]) 
                        values[i,:n_variables_left] = values[i,:n_variables_left] - \
                            delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value # solve FVM equations
                    
                    # Evolution equation for the cell with index right_boundary_subdomain+1
                    fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain,:n_variables_left],
                        values[right_boundary_subdomain+1,:n_variables_left],
                        system_matrix_left,
                        'positive',
                        delta_t,
                        delta_x)             
                    fluctuation_minus_Full = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain+1,:n_variables_left],
                        values[right_boundary_subdomain+2,:n_variables_left],
                        system_matrix_left,'negative',
                        delta_t,
                        delta_x) 
                    fluctuation_minus_Restricted = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain+1,:n_variables_right],
                        values[right_boundary_subdomain+2,:n_variables_right],
                        system_matrix_right,'negative',
                        delta_t,
                        delta_x) 
                    source_term_value = source_term_left(values[right_boundary_subdomain+1,:n_variables_left]) 
                    values[right_boundary_subdomain+1,:n_variables_right] = (values[right_boundary_subdomain+1,:n_variables_right] 
                    -delta_t/delta_x*(fluctuation_plus[:n_variables_right]+fluctuation_minus_Restricted)
                    +delta_t*source_term_value[:n_variables_right]) # solve FVM equations for first moments
                    values[right_boundary_subdomain+1,n_variables_right:n_variables_left] = (values[right_boundary_subdomain+1,n_variables_right:n_variables_left] 
                    -delta_t/delta_x*(fluctuation_plus[n_variables_right:n_variables_left]+fluctuation_minus_Full[n_variables_right:n_variables_left])
                    +delta_t*source_term_value[n_variables_right:n_variables_left]) # solve FVM equations for last moments
                    
                    # Evolution equation for the cell with index right_boundary_subdomain+2
                    fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain+1,:n_variables_right],
                        values[right_boundary_subdomain+2,:n_variables_right],
                        system_matrix_right,'positive',
                        delta_t,
                        delta_x) 
                    fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain+2,:n_variables_right],
                        values[right_boundary_subdomain+3,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x) 
                    source_term_value = source_term_right(values[right_boundary_subdomain+2,:n_variables_right]) 
                    values[right_boundary_subdomain+2,:n_variables_right] = (values[right_boundary_subdomain+2,:n_variables_right]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value) # solve FVM equations

                    # Evolution equation for the cell with index right_boundary_subdomain+3
                    fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain+2,:n_variables_right],
                        values[right_boundary_subdomain+3,:n_variables_right],
                        system_matrix_right,'positive',
                        delta_t,
                        delta_x) 
                    fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                        values[right_boundary_subdomain+3,:n_variables_right],
                        values[right_boundary_subdomain+4,:n_variables_right],
                        system_matrix_right,
                        'negative',
                        delta_t,
                        delta_x) 
                    source_term_value = source_term_right(values[right_boundary_subdomain+3,:n_variables_right]) 
                    values[right_boundary_subdomain+3,:n_variables_right] = (values[right_boundary_subdomain+3,:n_variables_right]
                    -delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value) # solve FVM equations

                    right_boundary_subdomain += 3
            
            for i in range(right_boundary_subdomain+1,self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    values[i-1,:n_variables_right],
                    values[i,:n_variables_right],
                    system_matrix_right,
                    'positive',
                    delta_t,
                    delta_x) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    values[i,:n_variables_right],
                    values[i+1,:n_variables_right],
                    system_matrix_right,
                    'negative',
                    delta_t,
                    delta_x) 
                source_term_value = source_term_right(values[i,:n_variables_right]) 
                values[i,:n_variables_right] = values[i,:n_variables_right] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value # solve FVM equations
            t+=delta_t

            self._update_domain_decomposition(values, 0.01, 0.001)

        simulation_data = self._post_processing(values)
        return simulation_data
    
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
                                    values_boundary: np.array) -> np.array:

        """
        update the boundary conditions

        Parameters
        ----------
        values_boundary : numpy 1D array #TODO: implement boundary conditions that include more cells
            the values of the variables in the boundary cell
        
        Returns
        -------
        values_ghost: numpy 1D array
            the values of the variables in the ghost cell

        """

        if self.boundary_condition == 'INFLOW_OUTFLOW':
            values_ghost = values_boundary
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

    def _update_domain_decomposition(self,
                                     values,
                                     tolerance_up = 0.01,
                                     tolerance_down = 0.001):
        
        """
        updates the domain decomposition

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

        breakdown_criteria = self.pde_type.compute_breakdown_criterion(values, 'height_gradient',self.mesh.resolution)

        interface_left = 0
        for i in range(len(self.boundary_interfaces_discretized)):
            interface_right = self.boundary_interfaces_discretized[i]
            if np.max(breakdown_criteria[interface_left:interface_right]) > tolerance_up:
                self.orders[i] = min(self.orders[i]+1, 5)
                self.numbers_of_variables[i] = self.pde_type.compute_number_of_variables(self.orders[i])
            elif np.max(breakdown_criteria[interface_left:interface_right]) < tolerance_down:
                self.orders[i] = max(self.orders[i]-1, 0)
                self.numbers_of_variables[i] = self.pde_type.compute_number_of_variables(self.orders[i])
            interface_left = interface_right + 1
        
        if np.max(breakdown_criteria[interface_left:self.mesh.resolution]) > tolerance_up:
            self.orders[-1] = min(self.orders[-1]+1, 5)
            self.numbers_of_variables[-1] = self.pde_type.compute_number_of_variables(self.orders[-1])
        elif np.max(breakdown_criteria[interface_left:self.mesh.resolution]) < tolerance_down:
            self.orders[-1] = max(self.orders[-1]-1, 0)
            self.numbers_of_variables[-1] = self.pde_type.compute_number_of_variables(self.orders[-1])
    
    def _post_processing(self,
                         values: list) -> np.array:

        data_array = np.zeros((self.mesh.resolution,self.max_number_of_variables+1)) # rewrite this such that it can be generalized to other PDE models

        for i in range(self.mesh.resolution):
            data_array[i,0] = self.mesh.cell_center_positions[i]
        data_array[:,1] = values[1:-1,0]
        data_array[:,2] = np.divide(values[1:-1,1],data_array[:,1])
        for j in range(self.max_order): #TODO: this is unnecessary routine here
            data_array[:,j+3] = np.divide(values[1:-1,j+2],data_array[:,1])
        
        return data_array


class ClassicalGalerkinSimulation1D(Simulation):

    """
    This interface represents a classical (not spatially adaptive) simulation in 1D with intrusive uncertainty.

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
                 mom_order: int,
                 SG_order: int,
                 pde_type: PDE.PDE,
                 mesh: Mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 spatial_discretization: SpatialDiscretization.SpatialDiscretization):
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

        """
        self.mom_order = mom_order
        self.SG_order = SG_order
        self.pde_type = pde_type
        self.number_of_variables = pde_type.compute_number_of_variables(self.mom_order, self.SG_order)
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.spatial_discretization = spatial_discretization

    def run_simulation(self,
                       t_end: float,
                       g = 9.81) -> np.array:

        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)

        CFL = 0.7
        t = 0

        def system_matrix(cell_values):
            return self.pde_type.compute_system_matrix(self.mom_order, self.SG_order, cell_values)

        def source_term(cell_values):
            return self.pde_type.compute_source_term(self.mom_order, self.SG_order, cell_values)


        while t < t_end:

            # update boundary conditions
            values[0,:] = self._update_boundary_conditions(values[1,:])
            values[self.mesh.resolution+1,:] = self._update_boundary_conditions(values[self.mesh.resolution,:])
            
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
            max_wave_speed_min = np.max(np.abs(wave_speed_part2 - np.sqrt(wave_speed_part2)))
            max_speed = max(max_wave_speed_plus,max_wave_speed_min)

            delta_t = CFL*delta_x/max_speed

            for i in range(1, self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    values[i-1,:],
                    values[i,:],
                    system_matrix,
                    'positive',
                    delta_t,
                    delta_x) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    values[i,:],
                    values[i+1,:],
                    system_matrix,
                    'negative',
                    delta_t,
                    delta_x) 
                source_term_value = source_term(values[i,:]) 
                values[i,:] = values[i,:] - delta_t/delta_x*(fluctuation_plus+fluctuation_minus) + delta_t*source_term_value # solve FVM equations

            t += delta_t
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
                                   values_boundary: np.array) -> np.array:
        """
        update the boundary conditions

        Parameters
        ----------
        values_boundary : numpy 1D array #TODO: implement boundary conditions that include more cells
            the values of the variables in the boundary cell
        
        Returns
        -------
        values_ghost: numpy 1D array
            the values of the variables in the ghost cell

        """

        if self.boundary_condition == 'INFLOW_OUTFLOW':
            values_ghost = values_boundary
        
        else:
            print("This boundary condition is not implemented yet")
        
        return values_ghost 
    
    def _post_processing(self,
                         values) -> np.array:

        data_array = np.zeros((self.mesh.resolution, self.number_of_variables + 1)) # rewrite this such that it can be generalized to other PDE models
        data_array[:,0] = self.mesh.cell_center_positions
        data_array[:,1:self.number_of_variables+1] = values[1:-1,:]
        
        return data_array    