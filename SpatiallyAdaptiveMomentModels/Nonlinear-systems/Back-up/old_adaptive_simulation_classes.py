import numpy as np

class SpatiallyAdaptiveSimulation1D(Simulation,ABC):

    """
    This abstract class represents a spatially adaptive simulation in 1D.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell

    
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
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders    
    
    Instance methods
    ----------------
    def _update_domain_decomposition_pointwise(self,values,delta_x,delta_t):
        updates the domain decompositions in each point of the domain
    """

    def __init__(self,
                 start_order: list,
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 spatial_discretization: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration):

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
        self.max_order = start_order
        self.max_number_of_variables = pde_type.compute_number_of_variables(self.max_order)

        self.orders_cellwise = []
        self.numbers_of_variables_cellwise = []
        for i in range(self.mesh.resolution+2):
            self.orders_cellwise.append(start_order)
            self.numbers_of_variables_cellwise.append(self.pde_type.compute_number_of_variables(start_order))

        self.dom_decomp_val_res1 = np.zeros(self.mesh.resolution)
        self.dom_decomp_val_res2 = np.zeros(self.mesh.resolution)

        self.breakdown_estimators = np.zeros((self.mesh.resolution,self.max_order+4))
        self.breakdown_criteria_flags = np.full(shape=self.mesh.resolution,dtype=int,fill_value=0)
              
    def _update_boundary_conditions(self,
                                    values: np.ndarray,
                                    boundary: str) -> np.ndarray:

        """
        update the boundary conditions at the specified boundary

        Parameters
        ----------
        values : numpy 2D array 
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
                                cell_centers_x: np.ndarray) -> np.ndarray:
            
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
                initial_values[i+1,:self.numbers_of_variables[m]] = self.pde_type.get_initial_values(self.orders[m],self.initial_condition,cell_centers_x[i])  
        for i in range(right_boundary_subdomain,self.mesh.resolution):
            initial_values[i+1,:self.numbers_of_variables[-1]] = self.pde_type.get_initial_values(self.orders[-1],self.initial_condition,cell_centers_x[i])
            
        return initial_values
    
    def _post_processing(self,
                         values: list) -> np.ndarray:

        data_array = np.zeros((self.mesh.resolution,self.max_number_of_variables+2))

        for i in range(self.mesh.resolution):
            data_array[i,0] = self.mesh.cell_center_positions[i]

        data_array[:,1:-1] = values[1:-1,:]
        data_array[:,1:-1] = self.pde_type.convert_to_primitive(self.max_order,data_array[:,1:-1])
        data_array[:,-1] = self.orders_cellwise[1:-1]

        print("Orders at the end of the simulation:",self.orders)
        print("numbers of vars at the end of the simulation:",self.numbers_of_variables)
        
        return data_array   

class NonConservativeAdaptiveSimulation1D(SpatiallyAdaptiveSimulation1D):

    """
    This class represents an adaptive simulation in 1D that can not be written in conservative form.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell

    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

        
    Implemented methods from interface Simulation
    ----------------------------------------------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values

        
    Implemented methods from abstract class SpatiallyAdaptiveSimulation1D
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     
    
    Instance methods
    ----------------
    None
    """

    def run_simulation(self,
                       t_end: float,
                       g = 1) -> np.ndarray:
        
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids
        print("Orders at the beginning of the simulation:",self.orders)

        values = self._get_initial_conditions(self.mesh.cell_center_positions)
        fluctuations_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        fluctuations_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        res1_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        res1_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))

        CFL = 0.9
        
        step_count = 0
        t = 0

        while t < t_end:
            # update boundary conditions
            values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values,'right')
            max_speed = self.pde_type.compute_max_wavespeed(self.max_order,values)
      
            delta_t = CFL*delta_x/max_speed 

            # if step_count%10 == 0:
            #     values = self._reconstruct_subdomains(prev_values,delta_x)
            values = self._reconstruct_subdomains(values,delta_x,delta_t)
            prev_values = np.copy(values)

            right_boundary_subdomain = 0

            for m in range(len(self.boundary_interfaces_discretized)):
                order_left = self.orders[m]
                n_variables_left = self.numbers_of_variables[m]
                order_right = self.orders[m+1]
                n_variables_right = self.numbers_of_variables[m+1]

                def system_matrix_left(cell_values):
                    return self.pde_type.compute_system_matrix(order_left,cell_values)

                def system_matrix_right(cell_values):
                    return self.pde_type.compute_system_matrix(order_right,cell_values)

                left_boundary_subdomain = right_boundary_subdomain + 1
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]
                
                if order_right > order_left:
                    prev_values[right_boundary_subdomain-1,n_variables_left:n_variables_right] = \
                        prev_values[right_boundary_subdomain,n_variables_left:n_variables_right] # update boundary interface boundary condition 
                    for i in range(left_boundary_subdomain,right_boundary_subdomain):
                        generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            prev_values[i-1,:n_variables_left],
                            prev_values[i,:n_variables_left],
                            system_matrix_left,
                            delta_t,
                            delta_x)
                        fluctuations_plus[i-1,:n_variables_left] =\
                            generalized_roe_plus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left])
                        fluctuations_min[i-1,:n_variables_left] =\
                            generalized_roe_minus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left]) 
                        res1_min[i-1,:n_variables_left-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])
                        res1_plus[i-1,:n_variables_left-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])

                    # Fluctuation between cell with index right_boundary_subdomain-1 and cell with index right_boundary_subdomain
                    generalized_roe_minus_Full, generalized_roe_plus_Full = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[right_boundary_subdomain-1,:n_variables_right],
                        prev_values[right_boundary_subdomain,:n_variables_right],
                        system_matrix_right,
                        delta_t,
                        delta_x) 
                    generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[right_boundary_subdomain-1,:n_variables_left],
                        prev_values[right_boundary_subdomain,:n_variables_left],
                        system_matrix_left,
                        delta_t,
                        delta_x)  
                    fluctuation_plus_Full = generalized_roe_plus_Full@(prev_values[right_boundary_subdomain,:n_variables_right]\
                                                                       -prev_values[right_boundary_subdomain-1,:n_variables_right])
                    fluctuation_plus_Restricted = generalized_roe_plus@(prev_values[right_boundary_subdomain,:n_variables_left]\
                                                                                   -prev_values[right_boundary_subdomain-1,:n_variables_left])
                    fluctuations_plus[right_boundary_subdomain-1,:n_variables_right] \
                        = np.hstack((fluctuation_plus_Restricted, fluctuation_plus_Full[n_variables_left:n_variables_right]))
                    fluctuations_min[right_boundary_subdomain-1,:n_variables_left] = generalized_roe_minus@(prev_values[right_boundary_subdomain,:n_variables_left]\
                                                               -prev_values[right_boundary_subdomain-1,:n_variables_left])
                    
                    res1_min[right_boundary_subdomain-1,:n_variables_left-1] =\
                        generalized_roe_minus[:-1,-1]*(prev_values[right_boundary_subdomain,n_variables_left-1]-prev_values[right_boundary_subdomain-1,n_variables_left-1])
                    res1_plus[right_boundary_subdomain-1,:n_variables_left-1] =\
                        generalized_roe_plus[:-1,-1]*(prev_values[right_boundary_subdomain,n_variables_left-1]-prev_values[right_boundary_subdomain-1,n_variables_left-1])        

                    # Fluctuation between cell with index right_boundary_subdomain and cell with index right_boundary_subdomain+1
                    generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[right_boundary_subdomain,:n_variables_right],
                        prev_values[right_boundary_subdomain+1,:n_variables_right],
                        system_matrix_right,
                        delta_t,
                        delta_x)  
                    fluctuations_plus[right_boundary_subdomain,:n_variables_right] = generalized_roe_plus@(prev_values[right_boundary_subdomain+1,:n_variables_right]\
                                                                       -prev_values[right_boundary_subdomain,:n_variables_right])
                    fluctuations_min[right_boundary_subdomain,:n_variables_right] = generalized_roe_minus@(prev_values[right_boundary_subdomain+1,:n_variables_right]\
                                                               -prev_values[right_boundary_subdomain,:n_variables_right])
                    
                    res1_min[right_boundary_subdomain,:n_variables_left-1] =\
                        generalized_roe_minus[:n_variables_left-1,n_variables_left-1]\
                            *(prev_values[right_boundary_subdomain+1,n_variables_left-1]-prev_values[right_boundary_subdomain,n_variables_left-1])                    
                    res1_plus[right_boundary_subdomain,:n_variables_right-1] =\
                        generalized_roe_plus[:-1,-1]*(prev_values[right_boundary_subdomain+1,n_variables_right-1]-prev_values[right_boundary_subdomain,n_variables_right-1]) 
                    
                    right_boundary_subdomain += 1

                else:
                    prev_values[right_boundary_subdomain+2,n_variables_right:n_variables_left] = \
                        prev_values[right_boundary_subdomain+1,n_variables_right:n_variables_left] # update boundary interface boundary condition
                    for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                        generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            prev_values[i-1,:n_variables_left],
                            prev_values[i,:n_variables_left],
                            system_matrix_left,
                            delta_t,
                            delta_x)
                        fluctuations_plus[i-1,:n_variables_left] =\
                            generalized_roe_plus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left])
                        fluctuations_min[i-1,:n_variables_left] =\
                            generalized_roe_minus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left]) 
                        res1_min[i-1,:n_variables_left-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])
                        res1_plus[i-1,:n_variables_left-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])
                    
                    # Fluctuation between cell with index right_boundary_subdomain and cell with index right_boundary_subdomain+1
                    generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[right_boundary_subdomain,:n_variables_left],
                        prev_values[right_boundary_subdomain+1,:n_variables_left],
                        system_matrix_left,
                        delta_t,
                        delta_x)  
                    fluctuations_plus[right_boundary_subdomain,:n_variables_left] = generalized_roe_plus@(prev_values[right_boundary_subdomain+1,:n_variables_left]\
                                                                       -prev_values[right_boundary_subdomain,:n_variables_left])
                    fluctuations_min[right_boundary_subdomain,:n_variables_left] = generalized_roe_minus@(prev_values[right_boundary_subdomain+1,:n_variables_left]\
                                                               -prev_values[right_boundary_subdomain,:n_variables_left])
                    
                    res1_min[right_boundary_subdomain,:n_variables_left-1] =\
                        generalized_roe_plus[:-1,-1]*(prev_values[right_boundary_subdomain+1,n_variables_left-1]-prev_values[right_boundary_subdomain,n_variables_left-1]) 
                    res1_plus[right_boundary_subdomain,:n_variables_right-1] =\
                        generalized_roe_minus[:n_variables_right-1,n_variables_right-1]\
                            *(prev_values[right_boundary_subdomain+1,n_variables_right-1]-prev_values[right_boundary_subdomain,n_variables_right-1])

                    # Fluctuation between cell with index right_boundary_subdomain+1 and cell with index right_boundary_subdomain+2
                    generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[right_boundary_subdomain+1,:n_variables_right],
                        prev_values[right_boundary_subdomain+2,:n_variables_right],
                        system_matrix_right,
                        delta_t,
                        delta_x)             
                    generalized_roe_minus_Full, generalized_roe_plus_Full = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[right_boundary_subdomain+1,:n_variables_left],
                        prev_values[right_boundary_subdomain+2,:n_variables_left],
                        system_matrix_left,
                        delta_t,
                        delta_x) 
                    
                    fluctuation_minus_Full = generalized_roe_minus_Full@(prev_values[right_boundary_subdomain+2,:n_variables_left]\
                                                             -prev_values[right_boundary_subdomain+1,:n_variables_left])
                    fluctuation_minus_Restricted = generalized_roe_minus@(prev_values[right_boundary_subdomain+2,:n_variables_right]\
                                                             -prev_values[right_boundary_subdomain+1,:n_variables_right])
                    
                    fluctuations_plus[right_boundary_subdomain+1,:n_variables_right] = generalized_roe_plus@(prev_values[right_boundary_subdomain+2,:n_variables_right]\
                                                             -prev_values[right_boundary_subdomain+1,:n_variables_right])
                    fluctuations_min[right_boundary_subdomain+1,:n_variables_left] = np.hstack((fluctuation_minus_Restricted, fluctuation_minus_Full[n_variables_right:n_variables_left]))
                    
                    res1_min[right_boundary_subdomain+1,:n_variables_right-1] =\
                        generalized_roe_minus[:-1,-1]*(prev_values[right_boundary_subdomain+2,n_variables_right-1]-prev_values[right_boundary_subdomain+1,n_variables_right-1])
                    res1_plus[right_boundary_subdomain+1,:n_variables_right-1] =\
                        generalized_roe_plus[:-1,-1]*(prev_values[right_boundary_subdomain+2,n_variables_right-1]-prev_values[right_boundary_subdomain+1,n_variables_right-1])   

                    right_boundary_subdomain += 2
            
            for i in range(right_boundary_subdomain+1,self.mesh.resolution+2):
                generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                    prev_values[i-1,:n_variables_right],
                    prev_values[i,:n_variables_right],
                    system_matrix_right,
                    delta_t,
                    delta_x)
                fluctuations_plus[i-1,:n_variables_right] =\
                    generalized_roe_plus@(prev_values[i,:n_variables_right]-prev_values[i-1,:n_variables_right])
                fluctuations_min[i-1,:n_variables_right] =\
                    generalized_roe_minus@(prev_values[i,:n_variables_right]-prev_values[i-1,:n_variables_right]) 
                res1_min[i-1,:n_variables_right-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables_right-1]-prev_values[i-1,n_variables_right-1])
                res1_plus[i-1,:n_variables_right-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables_right-1]-prev_values[i-1,n_variables_right-1])

            right_boundary_subdomain = 0
            for m in range(len(self.boundary_interfaces_discretized)):

                order_left = self.orders[m]
                n_variables_left = self.numbers_of_variables[m]
                order_right = self.orders[m+1]
                n_variables_right = self.numbers_of_variables[m+1]

                def source_term_left(cell_values,delta_t):
                    return self.pde_type.compute_source_term(order_left,cell_values,delta_t)

                def source_term_right(cell_values,delta_t):
                    return self.pde_type.compute_source_term(order_right,cell_values,delta_t)

                left_boundary_subdomain = right_boundary_subdomain+1
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]
                
                if order_right > order_left:
                    for i in range(left_boundary_subdomain,right_boundary_subdomain):

                        values[i,:n_variables_left] = prev_values[i,:n_variables_left]\
                            - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables_left]+fluctuations_min[i,:n_variables_left]) 
                        values[i,:n_variables_left] = self.time_integration.integrate(values[i,:n_variables_left],source_term_left,delta_t)

                        self.dom_decomp_val_res1[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables_left-1]+res1_min[i,:n_variables_left-1])
                        self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_left-1]-prev_values[i,n_variables_left-1])

                    # Evolution equation for the cell with index right_boundary_subdomain

                    values[right_boundary_subdomain,:n_variables_right] = (prev_values[right_boundary_subdomain,:n_variables_right]
                    -delta_t/delta_x*(fluctuations_plus[right_boundary_subdomain-1,:n_variables_right] + fluctuations_min[right_boundary_subdomain,:n_variables_right])) 
                    
                    values[right_boundary_subdomain,:n_variables_right] =\
                        self.time_integration.integrate(values[right_boundary_subdomain,:n_variables_right],source_term_right,delta_t)
                    
                    self.dom_decomp_val_res1[right_boundary_subdomain-1] =\
                        np.linalg.norm(res1_plus[right_boundary_subdomain-1,:n_variables_left-1]+res1_min[right_boundary_subdomain,:n_variables_left-1])
                    self.dom_decomp_val_res2[right_boundary_subdomain-1] =\
                        np.abs(values[right_boundary_subdomain,n_variables_left-1]-prev_values[right_boundary_subdomain,n_variables_left-1])                     
            
                else:
                    for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                        values[i,:n_variables_left] = prev_values[i,:n_variables_left]\
                            -delta_t/delta_x*(fluctuations_plus[i-1,:n_variables_left]+fluctuations_min[i,:n_variables_left])
                        values[i,:n_variables_left] = self.time_integration.integrate(values[i,:n_variables_left],source_term_left,delta_t)

                        self.dom_decomp_val_res1[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables_left-1]+res1_min[i,:n_variables_left-1])
                        self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_left-1]-prev_values[i,n_variables_left-1])
                    
                    # Evolution equation for the cell with index right_boundary_subdomain+1
                    values[right_boundary_subdomain+1,:n_variables_left] = (prev_values[right_boundary_subdomain+1,:n_variables_left]
                    -delta_t/delta_x*(fluctuations_plus[right_boundary_subdomain,:n_variables_left] + fluctuations_min[right_boundary_subdomain+1,:n_variables_left])) 
                    
                    values[right_boundary_subdomain+1,:n_variables_left] =\
                        self.time_integration.integrate(values[right_boundary_subdomain+1,:n_variables_left],source_term_left,delta_t)
                    
                    self.dom_decomp_val_res1[right_boundary_subdomain] =\
                        np.linalg.norm(res1_plus[right_boundary_subdomain,:n_variables_right-1]+res1_min[right_boundary_subdomain+1,:n_variables_right-1])
                    self.dom_decomp_val_res2[right_boundary_subdomain] =\
                        np.abs(values[right_boundary_subdomain+1,n_variables_right-1]-prev_values[right_boundary_subdomain+1,n_variables_right-1]) 

                    right_boundary_subdomain += 1
            
            for i in range(right_boundary_subdomain+1,self.mesh.resolution+1):
                values[i,:n_variables_right] = prev_values[i,:n_variables_right]\
                    - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables_right]+fluctuations_min[i,:n_variables_right]) 
                values[i,:n_variables_right] = self.time_integration.integrate(values[i,:n_variables_right],source_term_right,delta_t)

                self.dom_decomp_val_res1[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables_right-1]+res1_min[i,:n_variables_right-1])
                self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_right-1]-prev_values[i,n_variables_right-1])
            
            self.dom_decomp_val_res1 = self.dom_decomp_val_res1/delta_x
            self.dom_decomp_val_res2 = self.dom_decomp_val_res2/delta_t
            step_count += 1
            print(t)
            t+=delta_t

            
        simulation_data = self._post_processing(values)
        return simulation_data
    
    def _reconstruct_subdomains(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float) -> np.ndarray:

        boundary_interfaces = []
        orders_merged = []
        number_of_variables_merged = []
        orders_out = []
        number_of_variables_out = []

        super()._update_domain_decomposition_pointwise(values,delta_x,delta_t)

        smooth_par = 20
        for i in range(1,self.mesh.resolution-(smooth_par-1),smooth_par):
            local_order = max(self.orders_cellwise[i:i+smooth_par])
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            orders_merged.append(local_order)
            number_of_variables_merged.append(local_number_of_variables)
            self.orders_cellwise[i:i+smooth_par] = [local_order]*smooth_par
            self.numbers_of_variables_cellwise[i:i+smooth_par] = [local_number_of_variables]*smooth_par

        self.orders_cellwise[0] = self.orders_cellwise[1]
        self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]

        local_order = max(self.orders_cellwise[i:-1])
        local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
        orders_merged.append(local_order)
        number_of_variables_merged.append(local_number_of_variables)
        self.orders_cellwise[i:-1] = [local_order]*(self.mesh.resolution+1-i)
        self.numbers_of_variables_cellwise[i:-1] = [local_number_of_variables]*(self.mesh.resolution+1-i)

        # orders_merged[-1] = max(self.orders_cellwise[i:-1])
        # number_of_variables_merged[-1] = self.pde_type.compute_number_of_variables(orders_merged[-1])
        if self.boundary_condition == 'PERIODIC':
            local_order = max(orders_merged[-1],orders_merged[0])
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            orders_merged[-1] = local_order
            orders_merged[0] = local_order
            number_of_variables_merged[-1] = local_number_of_variables
            number_of_variables_merged[0] = local_number_of_variables
            self.orders_cellwise[i:-1] = [local_order]*(self.mesh.resolution+1-i)
            self.numbers_of_variables_cellwise[i:-1] = [local_number_of_variables]*(self.mesh.resolution+1-i)
            self.orders_cellwise[:smooth_par+1] = [local_order]*(smooth_par+1)
            self.numbers_of_variables_cellwise[:smooth_par+1] = [local_number_of_variables]*(smooth_par+1)

        orders_out.append(orders_merged[0])
        number_of_variables_out.append(number_of_variables_merged[0])
        for i in range(len(orders_merged)-1):
            if orders_merged[i] != orders_merged[i+1]:
                orders_out.append(orders_merged[i+1])
                number_of_variables_out.append(number_of_variables_merged[i+1])
                boundary_interfaces.append(smooth_par*(i+1))

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
                values[left_boundary:right_boundary,self.numbers_of_variables[m]:self.max_number_of_variables] = 0
                left_boundary = right_boundary
            else:
                values[left_boundary:right_boundary+2,self.numbers_of_variables[m]:self.max_number_of_variables] = 0
                left_boundary = right_boundary+2
        values[left_boundary:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables] = 0 

        return values

class ConservativeAdaptiveSimulation1D(SpatiallyAdaptiveSimulation1D):
    """
    This class represents an adaptive simulation in 1D that can 
    be written in conservative form.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell

    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

        
    Implemented methods from interface Simulation
    ----------------------------------------------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values

        
    Implemented methods from abstract class SpatiallyAdaptiveSimulation1D
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     
    
    Instance methods
    ----------------
    None
    """
    
    def run_simulation(self,
                       t_end: float,
                       g = 1) -> np.ndarray:
        
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        print("Orders at the beginning of the simulation:",self.orders)
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)
        fluctuations_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        fluctuations_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        res1_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        res1_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))

        CFL = 0.5
        
        step_count = 0
        t = 0

        while t < t_end:
            # update boundary conditions
            values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values,'right')
            max_speed = self.pde_type.compute_max_wavespeed(self.max_order,values)
      
            delta_t = CFL*delta_x/max_speed 

            values = self._reconstruct_subdomains(values,delta_x,delta_t)
            prev_values = np.copy(values)

            right_boundary_subdomain = 0

            for m in range(len(self.boundary_interfaces_discretized)):
                order_left = self.orders[m]
                n_variables_left = self.numbers_of_variables[m]
                order_right = self.orders[m+1]
                n_variables_right = self.numbers_of_variables[m+1]

                def system_matrix_left(cell_values):
                    return self.pde_type.compute_system_matrix(order_left,cell_values)

                def system_matrix_right(cell_values):
                    return self.pde_type.compute_system_matrix(order_right,cell_values)

                left_boundary_subdomain = right_boundary_subdomain+1
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]
                
                if order_right > order_left:
                    prev_values[right_boundary_subdomain,n_variables_left:n_variables_right] = 0 
                    values[right_boundary_subdomain,n_variables_left:n_variables_right] = 0 
                    values_boundary_help = prev_values[right_boundary_subdomain+1,:n_variables_right]
                    values_boundary_help[n_variables_left:n_variables_right] = 0
                    for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                        generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            prev_values[i-1,:n_variables_left],
                            prev_values[i,:n_variables_left],
                            system_matrix_left,
                            delta_t,
                            delta_x)
                        fluctuations_plus[i-1,:n_variables_left] =\
                            generalized_roe_plus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left])
                        fluctuations_min[i-1,:n_variables_left] =\
                            generalized_roe_minus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left]) 
                        res1_min[i-1,:n_variables_left-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])
                        res1_plus[i-1,:n_variables_left-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])

                    # Fluctuation between cell with index right_boundary_subdomain and cell with index right_boundary_subdomain+1
                    generalized_roe_minus1,generalized_roe_plus1 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[right_boundary_subdomain,:n_variables_right],
                        values_boundary_help,
                        system_matrix_right,
                        delta_t,
                        delta_x)
                    generalized_roe_minus2,generalized_roe_plus2 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        values_boundary_help,
                        prev_values[right_boundary_subdomain+1,:n_variables_right],
                        system_matrix_right,
                        delta_t,
                        delta_x
                    )

                    fluctuations_plus[right_boundary_subdomain,:n_variables_right] =\
                        generalized_roe_plus1@(values_boundary_help-prev_values[right_boundary_subdomain,:n_variables_right])\
                            +generalized_roe_plus2@(prev_values[right_boundary_subdomain+1,:n_variables_right]-values_boundary_help)
                    fluctuations_min[right_boundary_subdomain,:n_variables_right] =\
                        generalized_roe_minus1@(values_boundary_help-prev_values[right_boundary_subdomain,:n_variables_right])\
                            +generalized_roe_minus2@(prev_values[right_boundary_subdomain+1,:n_variables_right]-values_boundary_help)

                    res1_min[right_boundary_subdomain,:n_variables_left-1] =\
                        generalized_roe_minus1[:n_variables_left-1,n_variables_left-1]\
                            *(values_boundary_help[n_variables_left-1]-prev_values[right_boundary_subdomain,n_variables_left-1])\
                                +generalized_roe_minus2[:n_variables_left-1,n_variables_left-1]*(values[right_boundary_subdomain+1,:n_variables_left-1]-values_boundary_help[n_variables_left-1])                    
                    res1_plus[right_boundary_subdomain,:n_variables_right-1] =\
                        generalized_roe_plus1[:n_variables_right-1,n_variables_right-1]*(values_boundary_help[n_variables_right-1]-prev_values[right_boundary_subdomain,n_variables_right-1])\
                            +generalized_roe_plus2[:n_variables_right-1,n_variables_right-1]*(prev_values[right_boundary_subdomain+1,n_variables_right-1]-values_boundary_help[n_variables_right-1])
                    
                    right_boundary_subdomain += 1

                else:
                    prev_values[right_boundary_subdomain+1,n_variables_right:n_variables_left] = 0
                    values[right_boundary_subdomain+1,n_variables_right:n_variables_left] = 0
                    values_boundary_help = prev_values[right_boundary_subdomain,:n_variables_left]
                    values_boundary_help[n_variables_right:n_variables_left] = 0

                    for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                        generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                            prev_values[i-1,:n_variables_left],
                            prev_values[i,:n_variables_left],
                            system_matrix_left,
                            delta_t,
                            delta_x)
                        fluctuations_plus[i-1,:n_variables_left] =\
                            generalized_roe_plus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left])
                        fluctuations_min[i-1,:n_variables_left] =\
                            generalized_roe_minus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left]) 
                        res1_min[i-1,:n_variables_left-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])
                        res1_plus[i-1,:n_variables_left-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])

                    # Fluctuation between cell with index right_boundary_subdomain and cell with index right_boundary_subdomain+1
                    generalized_roe_minus1,generalized_roe_plus1 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[right_boundary_subdomain,:n_variables_left],
                        values_boundary_help,
                        system_matrix_left,
                        delta_t,
                        delta_x)
                    generalized_roe_minus2,generalized_roe_plus2 = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        values_boundary_help,
                        prev_values[right_boundary_subdomain+1,:n_variables_left],
                        system_matrix_left,
                        delta_t,
                        delta_x
                    )
                    fluctuations_plus[right_boundary_subdomain,:n_variables_left] =\
                        generalized_roe_plus1@(values_boundary_help-prev_values[right_boundary_subdomain,:n_variables_left])\
                            +generalized_roe_plus2@(prev_values[right_boundary_subdomain+1,:n_variables_left]-values_boundary_help)
                    fluctuations_min[right_boundary_subdomain,:n_variables_left] =\
                        generalized_roe_minus1@(values_boundary_help-prev_values[right_boundary_subdomain,:n_variables_left])\
                            +generalized_roe_minus2@(prev_values[right_boundary_subdomain+1,:n_variables_left]-values_boundary_help)

                    res1_min[right_boundary_subdomain,:n_variables_left-1] =\
                        generalized_roe_minus1[:n_variables_left-1,n_variables_left-1]\
                            *(values_boundary_help[n_variables_left-1]-prev_values[right_boundary_subdomain,n_variables_left-1])\
                                +generalized_roe_minus2[:n_variables_left-1,n_variables_left-1]*(values[right_boundary_subdomain+1,:n_variables_left-1]-values_boundary_help[n_variables_left-1])                    
                    res1_plus[right_boundary_subdomain,:n_variables_right-1] =\
                        generalized_roe_plus1[:n_variables_right-1,n_variables_right-1]*(values_boundary_help[n_variables_right-1]-prev_values[right_boundary_subdomain,n_variables_right-1])\
                            +generalized_roe_plus2[:n_variables_right-1,n_variables_right-1]*(prev_values[right_boundary_subdomain+1,n_variables_right-1]-values_boundary_help[n_variables_right-1])
                    
                    right_boundary_subdomain += 1

            for i in range(right_boundary_subdomain+1,self.mesh.resolution+2):
                generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                    prev_values[i-1,:n_variables_right],
                    prev_values[i,:n_variables_right],
                    system_matrix_right,
                    delta_t,
                    delta_x)
                fluctuations_plus[i-1,:n_variables_right] =\
                    generalized_roe_plus@(prev_values[i,:n_variables_right]-prev_values[i-1,:n_variables_right])
                fluctuations_min[i-1,:n_variables_right] =\
                    generalized_roe_minus@(prev_values[i,:n_variables_right]-prev_values[i-1,:n_variables_right]) 
                res1_min[i-1,:n_variables_right-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables_right-1]-prev_values[i-1,n_variables_right-1])
                res1_plus[i-1,:n_variables_right-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables_right-1]-prev_values[i-1,n_variables_right-1])

            right_boundary_subdomain = 0

            for m in range(len(self.boundary_interfaces_discretized)):
                order_left = self.orders[m]
                n_variables_left = self.numbers_of_variables[m]
                order_right = self.orders[m+1]
                n_variables_right = self.numbers_of_variables[m+1]

                def source_term_left(cell_values,delta_t):
                    return self.pde_type.compute_source_term(order_left,cell_values,delta_t)

                def source_term_right(cell_values,delta_t):
                    return self.pde_type.compute_source_term(order_right,cell_values,delta_t)

                left_boundary_subdomain = right_boundary_subdomain+1
                right_boundary_subdomain = self.boundary_interfaces_discretized[m]
                
                for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                    values[i,:n_variables_left] = prev_values[i,:n_variables_left]\
                        - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables_left]+fluctuations_min[i,:n_variables_left]) 
                    values[i,:n_variables_left] = self.time_integration.integrate(values[i,:n_variables_left],source_term_left,delta_t)

                    self.dom_decomp_val_res1[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables_left-1]+res1_min[i,:n_variables_left-1])
                    self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_left-1]-prev_values[i,n_variables_left-1])
                    
            for i in range(right_boundary_subdomain+1,self.mesh.resolution+1):
                values[i,:n_variables_right] = prev_values[i,:n_variables_right]\
                    - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables_right]+fluctuations_min[i,:n_variables_right]) 
                values[i,:n_variables_right] = self.time_integration.integrate(values[i,:n_variables_right],source_term_right,delta_t)

                self.dom_decomp_val_res1[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables_right-1]+res1_min[i,:n_variables_right-1])
                self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables_right-1]-prev_values[i,n_variables_right-1])
            
            step_count += 1
            print(t)
            self.dom_decomp_val_res1 = self.dom_decomp_val_res1/delta_x
            self.dom_decomp_val_res2 = self.dom_decomp_val_res2/delta_t
            t+=delta_t

        values = simulation_data = self._post_processing(values)
        return simulation_data  

    def _reconstruct_subdomains(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float) -> np.ndarray:
        

        super()._update_domain_decomposition_pointwise(values,delta_x,delta_t)

        if self.boundary_condition == 'PERIODIC':
            self.orders_cellwise[1] = max(self.orders_cellwise[1],self.orders_cellwise[-2])
            self.numbers_of_variables_cellwise[1] = self.pde_type.compute_number_of_variables(self.orders_cellwise[1])
            self.orders_cellwise[-2] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[-2] = self.numbers_of_variables_cellwise[1]
            self.orders_cellwise[0] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]
            self.orders_cellwise[-1] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[-1] = self.numbers_of_variables_cellwise[1]

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

class InterpolatedAdaptiveSimulation1D(SpatiallyAdaptiveSimulation1D):
    """
    This class represents an adaptive simulation in 1D that uses a path-conservative numerical scheme 
    (and spatial coupling), and an interpolation technique to fill in the moment values when the order is increased.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell
    smooth_par : float
        the size of the subregions, relative to the size of the entire domain
    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

    Implemented methods from abstract class SpatiallyAdaptiveSimulation1D
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     

    Overriden methods from abstract class SpatiallAdaptiveSimulation1D
    def __init__(self,start_order,pde_type,mesh.RectangularMesh,boundary_condition,initial_condition,
                    breakdown_criterion,spatial_discretization,time_integration):
        initializes the SmoothedSubdomainReconstruction object
    
    Instance methods
    ----------------
    None
    """

    def __init__(self,
                 start_order: list,
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 spatial_discretization: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration):

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

        self.pde_type = pde_type
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.breakdown_criterion = breakdown_criterion
        self.spatial_discretization = spatial_discretization
        self.time_integration = time_integration

        self.max_order = start_order
        self.max_number_of_variables = pde_type.compute_number_of_variables(self.max_order)

        self.orders_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        self.numbers_of_variables_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        for i in range(self.mesh.resolution+2):
            self.orders_cellwise[i] = start_order
            self.numbers_of_variables_cellwise[i] = self.pde_type.compute_number_of_variables(start_order)

        self.dom_decomp_val_res1 = np.zeros(self.mesh.resolution)
        self.dom_decomp_val_res2 = np.zeros(self.mesh.resolution)

        self.breakdown_estimators = np.zeros((self.mesh.resolution,self.max_order+4))

        self.smooth_par = 8
        self.orders = np.full(shape=self.smooth_par+1,fill_value=start_order,dtype=int)
        self.breakdown_criteria_flags_subdomains = np.full(shape=self.smooth_par+1,fill_value=0,dtype=int)
        self.breakdown_criteria_flags = np.full(shape=self.mesh.resolution,dtype=int,fill_value=0)
        self.numbers_of_variables = np.full(shape=self.smooth_par+1,fill_value=self.pde_type.compute_number_of_variables(start_order),dtype=int)
        self.boundary_interfaces = np.zeros(self.smooth_par,dtype=int)

        self.n_cells_subdomain = int(np.floor(self.mesh.resolution/self.smooth_par))
        self.subdomain_start = int(np.ceil(self.n_cells_subdomain/2))+1
        self.boundary_interfaces[0] = self.subdomain_start - 1
        for i in range(1,self.smooth_par):
            self.boundary_interfaces[i] = self.subdomain_start + i*self.n_cells_subdomain - 1

        if self.boundary_condition != 'PERIODIC':
            self.orders = np.full(shape=self.smooth_par,fill_value=start_order,dtype=int)
            self.breakdown_criteria_flags_subdomains = np.full(shape=self.smooth_par,fill_value=0,dtype=int)
            self.numbers_of_variables = np.full(shape=self.smooth_par,fill_value=self.pde_type.compute_number_of_variables(start_order),dtype=int)
            self.boundary_interfaces = np.zeros(self.smooth_par-1,dtype=int)
            for i in range(self.smooth_par-1):
                self.boundary_interfaces[i] = (i+1)*self.n_cells_subdomain - 1

    def _get_initial_conditions(self,
                                cell_centers_x: np.ndarray) -> np.ndarray:

        initial_values = np.zeros((self.mesh.resolution+2,self.max_number_of_variables))      

        l_bound_subdom = 1
        for m in range(len(self.orders)-1):
            r_bound_subdom = self.boundary_interfaces[m]
            for i in range(l_bound_subdom,r_bound_subdom+1):
                initial_values[i,:self.numbers_of_variables[m]] =\
                    self.pde_type.get_initial_values(self.orders[m],self.initial_condition,cell_centers_x[i-1])
            l_bound_subdom = r_bound_subdom + 1
        for i in range(l_bound_subdom,self.mesh.resolution+1):
            initial_values[i,:self.numbers_of_variables[-1]] =\
                self.pde_type.get_initial_values(self.orders[-1],self.initial_condition,cell_centers_x[i-1])
            
        return initial_values

    def run_simulation(self,
                       t_end: float,
                       g = 1) -> np.ndarray:
    
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        print("Orders at the beginning of the simulation:",self.orders)
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)
        fluctuations_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        fluctuations_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        res1_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        res1_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))

        reconstruct_subdomains = self._reconstruct_subdomains()
        interpolate_added_moments = self._interpolate_subdomains()

        CFL = 0.8
        
        step_count = 0
        t = 0

        while t < t_end:
            # update boundary conditions
            values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values,'right')
            max_speed = self.pde_type.compute_max_wavespeed(self.max_order,values)
      
            delta_t = CFL*delta_x/max_speed 

            prev_values = np.copy(values)

            l_bound_subdom = 1

            for m in range(len(self.orders)-1):

                order_left = self.orders[m]
                n_variables_left = self.numbers_of_variables[m]
                order_right = self.orders[m+1]
                n_variables_right = self.numbers_of_variables[m+1]

                r_bound_subdom = self.boundary_interfaces[m]

                # if order_left < order_right:
                #     prev_values[r_bound_subdom,n_variables_left:n_variables_right]=\
                #     prev_values[r_bound_subdom+1,n_variables_left:n_variables_right]                
                # elif order_left > order_right:
                #     prev_values[r_bound_subdom+1,n_variables_right:n_variables_left]=\
                #     prev_values[r_bound_subdom,n_variables_right:n_variables_left]             

                def system_matrix(cell_values):
                    return self.pde_type.compute_system_matrix(order_left,cell_values)

                for i in range(l_bound_subdom, r_bound_subdom + 1):
                    generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                        prev_values[i-1,:n_variables_left],
                        prev_values[i,:n_variables_left],
                        system_matrix,
                        delta_t,
                        delta_x)
                    fluctuations_plus[i-1,:n_variables_left] =\
                        generalized_roe_plus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left])
                    fluctuations_min[i-1,:n_variables_left] =\
                        generalized_roe_minus@(prev_values[i,:n_variables_left]-prev_values[i-1,:n_variables_left]) 
                    res1_min[i-1,:n_variables_left-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])
                    res1_plus[i-1,:n_variables_left-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables_left-1]-prev_values[i-1,n_variables_left-1])                       
                l_bound_subdom = r_bound_subdom + 1

            def system_matrix(cell_values):
                return self.pde_type.compute_system_matrix(order_right,cell_values)
            for i in range(l_bound_subdom, self.mesh.resolution+2):
                generalized_roe_minus, generalized_roe_plus = self.spatial_discretization.compute_generalized_roe_and_viscosity(
                    prev_values[i-1,:n_variables_right],
                    prev_values[i,:n_variables_right],
                    system_matrix,
                    delta_t,
                    delta_x)
                fluctuations_plus[i-1,:n_variables_right] =\
                    generalized_roe_plus@(prev_values[i,:n_variables_right]-prev_values[i-1,:n_variables_right])
                fluctuations_min[i-1,:n_variables_right] =\
                    generalized_roe_minus@(prev_values[i,:n_variables_right]-prev_values[i-1,:n_variables_right]) 
                res1_min[i-1,:n_variables_right-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables_right-1]-prev_values[i-1,n_variables_right-1])
                res1_plus[i-1,:n_variables_right-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables_right-1]-prev_values[i-1,n_variables_right-1])                    

            l_bound_subdom = 1

            for m in range(len(self.orders)-1):
                def source_term(cell_values,delta_t):
                    return self.pde_type.compute_source_term(self.orders[m],cell_values,delta_t)

                r_bound_subdom = self.boundary_interfaces[m]

                n_variables = self.numbers_of_variables[m]

                for i in range(l_bound_subdom,r_bound_subdom + 1):
                    values[i,:n_variables] = values[i,:n_variables]\
                        - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables]+fluctuations_min[i,:n_variables]) 
                    values[i,:n_variables] = self.time_integration.integrate(values[i,:n_variables],source_term,delta_t)

                    self.dom_decomp_val_res1[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables-1]+res1_min[i,:n_variables-1])
                    self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables-1]-prev_values[i,n_variables-1])

                l_bound_subdom = r_bound_subdom + 1
                    
            for i in range(l_bound_subdom,self.mesh.resolution+1):
                
                def source_term(cell_values,delta_t):
                    return self.pde_type.compute_source_term(self.orders[-1],cell_values,delta_t)
                n_variables = self.numbers_of_variables[-1]

                values[i,:n_variables] = values[i,:n_variables]\
                    - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables]+fluctuations_min[i,:n_variables]) 
                values[i,:n_variables] = self.time_integration.integrate(values[i,:n_variables],source_term,delta_t)

                self.dom_decomp_val_res1[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables-1]+res1_min[i,:n_variables-1])
                self.dom_decomp_val_res2[i-1] = np.abs(values[i,n_variables-1]-prev_values[i,n_variables-1])

            for m in range(len(self.boundary_interfaces)):
                n_variables_left = self.numbers_of_variables[m]
                n_variables_right = self.numbers_of_variables[m+1]

                if n_variables_left < n_variables_right:
                    self.dom_decomp_val_res1[self.boundary_interfaces[m]-1] =\
                        2*np.linalg.norm(res1_plus[self.boundary_interfaces[m]-1,:n_variables_left-1])
                else: 
                    self.dom_decomp_val_res1[self.boundary_interfaces[m]] =\
                        2*np.linalg.norm(res1_min[self.boundary_interfaces[m],:n_variables_right-1])
            
            step_count += 1
            print(t)
            # if min(self.orders) == 8:
            #     break
            values = reconstruct_subdomains(values,delta_x,delta_t)
            # values = interpolate_added_moments(values)
            self.dom_decomp_val_res1 = self.dom_decomp_val_res1/delta_x
            self.dom_decomp_val_res2 = self.dom_decomp_val_res2/delta_t
            t+=delta_t
        values = simulation_data = self._post_processing(values)
        return simulation_data  

    def _reconstruct_subdomains(self) -> np.ndarray:
        
        _reconstruct_subdomains_fun = self._reconstruct_subdomains_periodicBoundary if self.boundary_condition == 'PERIODIC'\
            else self._reconstruct_subdomains_nonPeriodicBoundary

        return _reconstruct_subdomains_fun

    def _reconstruct_subdomains_interior(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float) -> np.ndarray:

        self._update_domain_decomposition_pointwise(values,delta_x,delta_t)

        for i in range(1,len(self.boundary_interfaces)):
            self.breakdown_criteria_flags_subdomains[i] =\
                np.max(self.breakdown_criteria_flags[self.boundary_interfaces[i-1]:self.boundary_interfaces[i]])
            local_order = np.max(self.orders_cellwise[self.boundary_interfaces[i-1]+1:self.boundary_interfaces[i]+1])
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            self.orders[i] = local_order
            self.numbers_of_variables[i] = local_number_of_variables
            self.orders_cellwise[self.boundary_interfaces[i-1]+1:self.boundary_interfaces[i]+1] = local_order
            self.numbers_of_variables_cellwise[self.boundary_interfaces[i-1]+1:self.boundary_interfaces[i]+1] = local_number_of_variables

    def _reconstruct_subdomains_nonPeriodicBoundary(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float) -> np.ndarray:

        self._reconstruct_subdomains_interior(values,delta_x,delta_t)

        self.orders[0] = np.max(self.orders_cellwise[1:self.subdomain_start])
        self.breakdown_criteria_flags_subdomains[0] = np.max(self.breakdown_criteria_flags[1:self.subdomain_start])
        self.orders[-1] = np.max(self.orders_cellwise[self.boundary_interfaces[-1]:-1])
        self.breakdown_criteria_flags_subdomains[-1] = np.max(self.breakdown_criteria_flags[self.boundary_interfaces[-1]:-1])
        self.numbers_of_variables[0] = self.pde_type.compute_number_of_variables(self.orders[0])
        self.numbers_of_variables[-1] = self.pde_type.compute_number_of_variables(self.orders[-1])

        values = self._process_domain_decomposition(values)

        return values

    def _reconstruct_subdomains_periodicBoundary(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float) -> np.ndarray:

        self._reconstruct_subdomains_interior(values,delta_x,delta_t)

        local_order = max(np.max(self.orders_cellwise[1:self.subdomain_start]),\
                          np.max(self.orders_cellwise[self.boundary_interfaces[-1]+1:-1]))
        self.breakdown_criteria_flags_subdomains[0] = max(np.max(self.breakdown_criteria_flags[1:self.subdomain_start]),\
            np.max(self.breakdown_criteria_flags[self.boundary_interfaces[-1]:-1]))
        self.breakdown_criteria_flags_subdomains[-1] = self.breakdown_criteria_flags_subdomains[0]
        local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
        self.orders[0] = local_order
        self.orders[-1] = local_order
        self.numbers_of_variables[0] = local_number_of_variables
        self.numbers_of_variables[-1] = local_number_of_variables
        self.orders_cellwise[:self.subdomain_start] = local_order
        self.orders_cellwise[self.boundary_interfaces[-1]+1:] = local_order
        self.numbers_of_variables_cellwise[:self.subdomain_start] = local_number_of_variables
        self.numbers_of_variables_cellwise[self.boundary_interfaces[-1]+1:] = local_number_of_variables

        values = self._process_domain_decomposition(values)

        return values

    def _process_domain_decomposition(self,values):

        # Set undefined moments to zero
        values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[0]:] = 0
        for i in range(1,len(self.boundary_interfaces)):
            values[self.boundary_interfaces[i-1]+1:self.boundary_interfaces[i]+1,self.numbers_of_variables[i]:] = 0
        values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[-1]:] = 0

        return values

    def _interpolate_subdomains(self) -> np.ndarray:
        
        _interpolate_subdomains_fun = self._linear_interpolate_subdomains_periodicBoundary if self.boundary_condition == 'PERIODIC' \
            else self._linear_interpolate_subdomains_nonPeriodicBoundary

        return _interpolate_subdomains_fun   

    def _interpolate_subdomains_interior(self,values):
        
        for i in range(2,len(self.boundary_interfaces)-1):
            if self.breakdown_criteria_flags_subdomains[i] > 0:
                values[self.boundary_interfaces[i-1]+1:self.boundary_interfaces[i]+1,\
                       self.numbers_of_variables[i]-self.breakdown_criteria_flags_subdomains[i]\
                            :self.numbers_of_variables[i]] =\
                    self._interpolate(values[self.boundary_interfaces[i-1]+1:self.boundary_interfaces[i]+1,
                                             self.numbers_of_variables[i]-self.breakdown_criteria_flags_subdomains[i]\
                                                :self.numbers_of_variables[i]],
                                      values[self.boundary_interfaces[i-2]+1:self.boundary_interfaces[i-1]+1,
                                             self.numbers_of_variables[i]-self.breakdown_criteria_flags_subdomains[i]\
                                                :self.numbers_of_variables[i]],
                                      values[self.boundary_interfaces[i]+1:self.boundary_interfaces[i+1]+1,
                                             self.numbers_of_variables[i]-self.breakdown_criteria_flags_subdomains[i]\
                                                :self.numbers_of_variables[i]])

        return values

    def _interpolate_subdomains_periodicBoundary(self,values):

        values = self._interpolate_subdomains_interior(values)

        if self.breakdown_criteria_flags_subdomains[0] > 0:
            interpol_values =\
                self._interpolate(np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[0]\
                                            -self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]],\
                                        values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[0]\
                                            -self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]]), axis=0),\
                        values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,\
                        self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:\
                            self.numbers_of_variables[0]],\
                        values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                        self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:\
                            self.numbers_of_variables[0]])
            values[self.boundary_interfaces[-1]+1:,\
                   self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]] =\
                   interpol_values[:self.mesh.resolution+1-self.boundary_interfaces[-1],:]
            values[:self.boundary_interfaces[0]+1,\
                   self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]] =\
                   interpol_values[self.mesh.resolution+1-self.boundary_interfaces[-1]:,:]            

        if self.smooth_par > 2:
            if self.breakdown_criteria_flags_subdomains[1] > 0:
                values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                    self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]] =\
                    self._interpolate(values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]],
                                      np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]],\
                                                     values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]]),axis=0),
                                      values[self.boundary_interfaces[1]+1:self.boundary_interfaces[2]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]])
                
            if self.breakdown_criteria_flags[-2] > 0:
                values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,\
                    self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:self.numbers_of_variables[-2]] =\
                    self._interpolate(values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]],
                                      values[self.boundary_interfaces[-3]+1:self.boundary_interfaces[-2]+1,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]],
                                      np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[-2]\
                                                        -self.breakdown_criteria_flags_subdomains[-2]:self.numbers_of_variables[-2]],\
                                                     values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[-2]\
                                                        -self.breakdown_criteria_flags_subdomains[-2]:self.numbers_of_variables[-2]]),axis=0))
        else:
            if self.breakdown_criteria_flags_subdomains[1] > 0:
                values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                    self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]] =\
                    self._interpolate(values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]],
                                      np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]],\
                                                     values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]]),axis=0),
                                      np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]],\
                                                     values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]]),axis=0))

        return values            

    def _interpolate_subdomains_nonPeriodicBoundary(self,values):

        values = self._interpolate_subdomains_interior(values)

        if self.breakdown_criteria_flags_subdomains[0] > 0:
            values[:self.boundary_interfaces[0]+1,\
                   self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]] =\
                        self._interpolate_from_right_data(values[:self.boundary_interfaces[0]+1,\
                                                            self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:\
                                                            self.numbers_of_variables[0]],
                                                        values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                                                            self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:\
                                                            self.numbers_of_variables[0]])
        if self.breakdown_criteria_flags_subdomains[-1] > 0:
            values[self.boundary_interfaces[-1]+1:,\
                   self.numbers_of_variables[-1]-self.breakdown_criteria_flags_subdomains[-1]:self.numbers_of_variables[-1]] =\
                        self._interpolate_from_left_data(values[self.boundary_interfaces[-1]+1:,\
                                                            self.numbers_of_variables[-1]-self.breakdown_criteria_flags_subdomains[-1]:\
                                                            self.numbers_of_variables[-1]],
                                                        values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,\
                                                            self.numbers_of_variables[-1]-self.breakdown_criteria_flags_subdomains[-1]:\
                                                            self.numbers_of_variables[-1]])

        if self.smooth_par > 2:
            if self.breakdown_criteria_flags_subdomains[1] > 0:
                values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                    self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]] =\
                    self._interpolate(values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]],
                                        values[:self.boundary_interfaces[0]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]],
                                        values[self.boundary_interfaces[1]+1:self.boundary_interfaces[2]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]])
        if self.smooth_par > 3:
            if self.breakdown_criteria_flags_subdomains[-2] > 0:
                values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,\
                    self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:self.numbers_of_variables[-2]] =\
                    self._interpolate(values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]],
                                        values[self.boundary_interfaces[-3]+1:self.boundary_interfaces[-2]+1,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]],
                                        values[self.boundary_interfaces[-1]+1:,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]])

        return values    

    def _interpolate(self,values,interpolation_data_left,interpolation_data_right):

        # data_points_x = np.concatenate((np.arange(interpolation_data_left.shape[0]),\
        #                                np.arange(interpolation_data_left.shape[0]+values.shape[0],
        #                                          interpolation_data_left.shape[0]+values.shape[0]+interpolation_data_right.shape[0],1)))
        # data_points_y = np.concatenate((interpolation_data_left,interpolation_data_right),axis=0)

        # interpolators = [BarycentricInterpolator(data_points_x, data_points_y[:, j]) for j in range(data_points_y.shape[1])]
        # interpolated_values = np.column_stack([interp(
        #     np.arange(interpolation_data_left.shape[0],interpolation_data_left.shape[0]+values.shape[0],1)
        # ) for interp in interpolators])

        # return interpolated_values
        return self._linear_interpolate(values,interpolation_data_left,interpolation_data_right)

    def _interpolate_from_left_data(self,values,interpolation_data):

        # data_points_x = np.arange(interpolation_data.shape[0])

        # interpolators = [BarycentricInterpolator(data_points_x, interpolation_data[:, j]) for j in range(interpolation_data.shape[1])]
        # interpolated_values = np.column_stack([interp(
        #     np.arange(interpolation_data.shape[0],interpolation_data.shape[0]+values.shape[0],1)
        # ) for interp in interpolators])

        # return interpolated_values   
        return self._linear_interpolate_from_left_data(values,interpolation_data)

    def _interpolate_from_right_data(self,values,interpolation_data):
        # data_points_x = np.arange(values.shape[0],values.shape[0]+interpolation_data.shape[0],1)

        # interpolators = [BarycentricInterpolator(data_points_x, interpolation_data[:, j]) for j in range(interpolation_data.shape[1])]
        # interpolated_values = np.column_stack([interp(
        #     np.arange(values.shape[0])
        # ) for interp in interpolators])

        # return interpolated_values  
        return self._linear_interpolate_from_right_data(values,interpolation_data)

    def _linear_interpolate_subdomains_interior(self,values):
        
        for i in range(1,len(self.boundary_interfaces)):
            if self.breakdown_criteria_flags_subdomains[i] > 0:
                values[self.boundary_interfaces[i-1]+1:self.boundary_interfaces[i]+1,\
                       self.numbers_of_variables[i]-self.breakdown_criteria_flags_subdomains[i]\
                            :self.numbers_of_variables[i]] =\
                    self._linear_interpolate(
                        values[self.boundary_interfaces[i-1],\
                            self.numbers_of_variables[i]-self.breakdown_criteria_flags_subdomains[i]:self.numbers_of_variables[i]],
                        values[self.boundary_interfaces[i]+1,\
                            self.numbers_of_variables[i]-self.breakdown_criteria_flags_subdomains[i]:self.numbers_of_variables[i]],
                        1+self.n_cells_subdomain
                        )

        return values

    def _linear_interpolate_subdomains_periodicBoundary(self,values):

        values = self._linear_interpolate_subdomains_interior(values)
        if self.breakdown_criteria_flags_subdomains[0] > 0:
            interpol_values = self._linear_interpolate(
                values[self.boundary_interfaces[-1],\
                    self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]],
                values[self.boundary_interfaces[0]+1,
                    self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]],
                self.mesh.resolution-self.boundary_interfaces[-1]+self.boundary_interfaces[0]+3
            )
            
            values[self.boundary_interfaces[-1]+1:,\
                   self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]] =\
                   interpol_values[:self.mesh.resolution+1-self.boundary_interfaces[-1],:]
            values[:self.boundary_interfaces[0]+1,\
                   self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]] =\
                   interpol_values[self.mesh.resolution+1-self.boundary_interfaces[-1]:,:] 

        return values                     

    def _linear_interpolate_subdomains_nonPeriodicBoundary(self,values):

        values = self._interpolate_subdomains_interior(values)

        if self.breakdown_criteria_flags_subdomains[0] > 0:
            values[:self.boundary_interfaces[0]+1,\
                   self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]] =\
                        self._interpolate_from_right_data(values[:self.boundary_interfaces[0]+1,\
                                                            self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:\
                                                            self.numbers_of_variables[0]],
                                                        values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                                                            self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:\
                                                            self.numbers_of_variables[0]])
        if self.breakdown_criteria_flags_subdomains[-1] > 0:
            values[self.boundary_interfaces[-1]+1:,\
                   self.numbers_of_variables[-1]-self.breakdown_criteria_flags_subdomains[-1]:self.numbers_of_variables[-1]] =\
                        self._interpolate_from_left_data(values[self.boundary_interfaces[-1]+1:,\
                                                            self.numbers_of_variables[-1]-self.breakdown_criteria_flags_subdomains[-1]:\
                                                            self.numbers_of_variables[-1]],
                                                        values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,\
                                                            self.numbers_of_variables[-1]-self.breakdown_criteria_flags_subdomains[-1]:\
                                                            self.numbers_of_variables[-1]])

        if self.smooth_par > 2:
            if self.breakdown_criteria_flags_subdomains[1] > 0:
                values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                    self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]] =\
                    self._interpolate(values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]],
                                        values[:self.boundary_interfaces[0]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]],
                                        values[self.boundary_interfaces[1]+1:self.boundary_interfaces[2]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]])
        if self.smooth_par > 3:
            if self.breakdown_criteria_flags_subdomains[-2] > 0:
                values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,\
                    self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:self.numbers_of_variables[-2]] =\
                    self._interpolate(values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]],
                                        values[self.boundary_interfaces[-3]+1:self.boundary_interfaces[-2]+1,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]],
                                        values[self.boundary_interfaces[-1]+1:,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]])

        return values 

    def _linear_interpolate(self,value_left,value_right,delta_i):

        slope = (value_right - value_left)/(delta_i)
        # interpolated_values = value_right + np.outer(np.arange(1,delta_i), slope)
        interpolated_values = np.zeros((delta_i-1,1))

        return interpolated_values

    def _linear_interpolate_from_left_data(self,values,interpolation_data):

        slope = (interpolation_data[-1,:] - interpolation_data[0,:])/(1+interpolation_data.shape[0])
        interpolated_values = interpolation_data[0,:] + np.outer(np.arange(1,values.shape[0]+1), slope)

        return interpolated_values 

    def _linear_interpolate_from_right_data(self,values,interpolation_data):

        slope = (interpolation_data[-1,:] - interpolation_data[0,:])/(1+interpolation_data.shape[0])
        interpolated_values = interpolation_data[0,:] + np.outer(np.arange(-values.shape[0],0), slope)

        return interpolated_values  

class SmoothedSubdomainReconstruction(SpatiallyAdaptiveSimulation1D):

    """
    This class represents a smoothed adaptive simulation in 1D. 
    It smooths the domain decomposition by grouping cells together and giving them the same model order.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell
    smooth_par : float
        the size of the subregions, relative to the size of the entire domain
    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

    Implemented methods from abstract class SpatiallyAdaptiveSimulation1D
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     

    Overriden methods from abstract class SpatiallAdaptiveSimulation1D
    def __init__(self,start_order,pde_type,mesh.RectangularMesh,boundary_condition,initial_condition,
                    breakdown_criterion,spatial_discretization,time_integration):
        initializes the SmoothedSubdomainReconstruction object
    
    Instance methods
    ----------------
    None
    """

    def __init__(self,
                 start_order: list,
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 spatial_discretization: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration):

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

        self.orders = [start_order,start_order,start_order]
        self.pde_type = pde_type
        self.numbers_of_variables = [pde_type.compute_number_of_variables(start_order),
                                     pde_type.compute_number_of_variables(start_order),
                                     pde_type.compute_number_of_variables(start_order)]
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.breakdown_criterion = breakdown_criterion
        self.spatial_discretization = spatial_discretization
        self.time_integration = time_integration

        self.boundary_interfaces_discretized = [np.floor_divide(2*self.mesh.resolution,3),np.floor_divide(self.mesh.resolution,3)]
        self.max_order = start_order
        self.max_number_of_variables = pde_type.compute_number_of_variables(self.max_order)

        self.orders_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        self.numbers_of_variables_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        for i in range(self.mesh.resolution+2):
            self.orders_cellwise[i] = start_order
            self.numbers_of_variables_cellwise[i] = self.pde_type.compute_number_of_variables(start_order)

        self.dom_decomp_val_res1 = np.zeros(self.mesh.resolution)
        self.dom_decomp_val_res2 = np.zeros(self.mesh.resolution)

        self.breakdown_estimators = np.zeros((self.mesh.resolution,self.max_order+4))
        self.breakdown_criteria_flags = np.full(shape=self.mesh.resolution,dtype=int,fill_value=0)

        self.smooth_par = 20
        self.orders = np.full(shape=self.smooth_par,fill_value=start_order,dtype=int)
        self.numbers_of_variables = np.full(shape=self.smooth_par,fill_value=self.pde_type.compute_number_of_variables(start_order),dtype=int)

        self.n_cells_subdomain = int(np.floor(self.mesh.resolution/self.smooth_par))
        self.subdomain_start = int(np.ceil(self.n_cells_subdomain/2))

    def _reconstruct_subdomains(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float) -> np.ndarray:

        boundary_interfaces = []
        orders_out = []
        number_of_variables_out = []

        # self._update_domain_decomposition_pointwise(values,delta_x,delta_t)

        k = 1
        for i in range(self.subdomain_start,self.mesh.resolution-self.n_cells_subdomain,self.n_cells_subdomain):
            local_order = np.max(self.orders_cellwise[i:i+self.n_cells_subdomain])
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            self.orders[k] = local_order
            self.numbers_of_variables[k] = local_number_of_variables
            self.orders_cellwise[i:i+self.n_cells_subdomain] = local_order
            self.numbers_of_variables_cellwise[i:i+self.n_cells_subdomain] = local_number_of_variables
            k += 1

        if self.boundary_condition == 'PERIODIC':
            local_order = max(np.max(self.orders_cellwise[1:self.subdomain_start]),np.max(self.orders_cellwise[i+self.smooth_par:-1]))
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            self.orders[0] = local_order
            self.numbers_of_variables[0] = local_number_of_variables
            self.orders_cellwise[:self.subdomain_start] = local_order
            self.orders_cellwise[i+self.smooth_par:] = local_order
            self.numbers_of_variables_cellwise[:self.subdomain_start] = local_number_of_variables
            self.numbers_of_variables_cellwise[i+self.smooth_par:] = local_number_of_variables

            orders_out.append(int(self.orders[0]))
            number_of_variables_out.append(int(self.numbers_of_variables[0]))
            for i in range(self.smooth_par-1):
                if self.orders[i] != self.orders[i+1]:
                    orders_out.append(int(self.orders[i+1]))
                    number_of_variables_out.append(int(self.numbers_of_variables[i+1]))
                    boundary_interfaces.append(self.subdomain_start+self.n_cells_subdomain*(i+1)-1)
            if self.orders[-1] != self.orders[0]:
                boundary_interfaces.append(self.subdomain_start+self.n_cells_subdomain*(i+1)-1)
            orders_out.append(int(self.orders[0]))
            number_of_variables_out.append(int(self.numbers_of_variables[0]))
        else:
            self.orders_cellwise[:self.subdomain_start] = self.orders[1]
            self.orders_cellwise[i+self.smooth_par:] = self.orders[-1]
            self.numbers_of_variables_cellwise[:self.subdomain_start] = self.numbers_of_variables[1]
            self.numbers_of_variables_cellwise[i+self.smooth_par:] = self.numbers_of_variables[-1]

            orders_out.append(int(self.orders[1]))
            number_of_variables_out.append(int(self.numbers_of_variables[1]))
            for i in range(1,self.smooth_par-1):
                if self.orders[i] != self.orders[i+1]:
                    orders_out.append(int(self.orders[i+1]))
                    number_of_variables_out.append(int(self.numbers_of_variables[i+1]))
                    boundary_interfaces.append(self.subdomain_start+self.n_cells_subdomain*(i+1)-1)

        self.orders = orders_out
        self.boundary_interfaces_discretized = boundary_interfaces
        self.numbers_of_variables = number_of_variables_out
        if len(self.boundary_interfaces_discretized) == 0:
            self.orders.append(self.orders[0])
            self.orders.append(self.orders[0])
            self.numbers_of_variables.append(int(self.numbers_of_variables[0]))
            self.numbers_of_variables.append(int(self.numbers_of_variables[0]))
            self.boundary_interfaces_discretized.append(np.floor_divide(self.mesh.resolution,3))
            self.boundary_interfaces_discretized.append(np.floor_divide(2*self.mesh.resolution,3))

        # Set undefined moments to zero
        left_boundary = 0
        for m in range(len(self.boundary_interfaces_discretized)):
            right_boundary = self.boundary_interfaces_discretized[m]
            if self.numbers_of_variables[m+1] > self.numbers_of_variables[m]:
                values[left_boundary:right_boundary,self.numbers_of_variables[m]:self.max_number_of_variables] = 0
                left_boundary = right_boundary
            else:
                values[left_boundary:right_boundary+2,self.numbers_of_variables[m]:self.max_number_of_variables] = 0
                left_boundary = right_boundary+2
        values[left_boundary:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables] = 0 

        return values

class SmoothedConsAdaptiveSimulation1D(SmoothedSubdomainReconstruction,ConservativeAdaptiveSimulation1D):

    """
    This class represents a smoothed adaptive simulation in 1D that uses the conservative interface flux coupling.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell
    smooth_par : float
        the size of the subregions, relative to the size of the entire domain
    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

        
    Inherited methods from parent class ConservativeAdaptiveSimulation
    -------------------------------------------------------------------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values

    Inherited methods from parent class SmoothedSubdomainReconstruction
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     
    def __init__(self,start_order,pde_type,mesh.RectangularMesh,boundary_condition,initial_condition,
                    breakdown_criterion,spatial_discretization,time_integration):
        initializes the SmoothedConsAdaptiveSimulation1D object
    
    Instance methods
    ----------------
    None
    """
    pass

class SmoothedNonConsAdaptiveSimulation1D(SmoothedSubdomainReconstruction,NonConservativeAdaptiveSimulation1D):

    """
    This class represents a smoothed adaptive simulation in 1D that uses the nonconservative padded buffer cell coupling.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell
    smooth_par : float
        the size of the subregions, relative to the size of the entire domain
    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

        
    Inherited methods from parent class ConservativeAdaptiveSimulation
    -------------------------------------------------------------------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values

    Inherited methods from parent class SmoothedSubdomainReconstruction
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     
    def __init__(self,start_order,pde_type,mesh.RectangularMesh,boundary_condition,initial_condition,
                    breakdown_criterion,spatial_discretization,time_integration):
        initializes the SmoothedNonConsAdaptiveSimulation1D object
    
    Instance methods
    ----------------
    None
    """
    pass

class ModelAdaptiveSimulation1D(SpatiallyAdaptiveSimulation1D):
    """
    This class represents a model-adaptive simulation that can be written in path-conservative form.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell

    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

        
    Implemented methods from interface Simulation
    ----------------------------------------------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values

        
    Implemented methods from abstract class SpatiallyAdaptiveSimulation1D
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     
    
    Instance methods
    ----------------
    None
    """
    
    """
    This abstract class represents a spatially adaptive simulation in 1D.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell

    
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
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders    
    
    Instance methods
    ----------------
    def _update_domain_decomposition_pointwise(self,values,delta_x,delta_t):
        updates the domain decompositions in each point of the domain
    """

    def __init__(self,
                 start_order: list,
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 spatial_discretization: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration):

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
        self.spatial_discretization_predictor = spatialDiscretization.LF(1)
        self.spatial_discretization_interface = self.spatial_discretization
        self.path_conservative_interface_coupling = True
        self.time_integration = time_integration

        self.boundary_interfaces_discretized = [np.floor_divide(self.mesh.resolution,2)]
        self.max_order = start_order
        self.max_number_of_variables = pde_type.compute_number_of_variables(self.max_order)

        # self.orders_cellwise = []
        # self.numbers_of_variables_cellwise = []
        # for i in range(self.mesh.resolution+2):
        #     self.orders_cellwise.append(start_order)
        #     self.numbers_of_variables_cellwise.append(self.pde_type.compute_number_of_variables(start_order))

        self.orders_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        self.numbers_of_variables_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        for i in range(self.mesh.resolution+2):
            self.orders_cellwise[i] = start_order
            self.numbers_of_variables_cellwise[i] = self.pde_type.compute_number_of_variables(start_order)

        self.breakdown_estimators = np.zeros((self.mesh.resolution,self.max_order+4))
        self.breakdown_criteria_flags = np.full(shape=self.mesh.resolution,dtype=int,fill_value=0)

    def run_simulation(self,
                       t_end: float) -> np.ndarray:
        
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        print("Orders at the beginning of the simulation:",self.orders)
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)

        CFL = 0.7
        
        step_count = 0
        t = 0

        while t < t_end:
            # update boundary conditions
            values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values,'right')
            # max_speed = self.pde_type.compute_max_wavespeed(self.max_order,values)
            max_speed = self.pde_type.compute_max_wavespeed(max(self.orders),values)
      
            delta_t = CFL*delta_x/max_speed 

            prev_values = np.copy(values)
            predicted_values = np.copy(values)

            padded_vectors_left,padded_vectors_right = self._construct_padded_vectors(prev_values)
            
            if t > 0:
                predicted_values_transport = self._transport_step_augmented(values,
                                                        delta_t,
                                                        delta_x,
                                                        self.spatial_discretization_predictor)
                predicted_values = self._source_step_augmented(predicted_values_transport,
                                                    delta_t,
                                                    delta_x)

                # padded_vectors_predicted_left,padded_vectors_predicted_right = self._construct_padded_vectors(predicted_values)

                # estimators_increase_prev,flags_increase_prev = self.pde_type.compute_breakdown_criteria_increase(
                #                                                             prev_values,
                #                                                             # padded_vectors_left,
                #                                                             # padded_vectors_right,
                #                                                             self.mesh.resolution,
                #                                                             delta_x,
                #                                                             delta_t,
                #                                                             self.max_order,
                #                                                             self.orders_cellwise,
                #                                                             self.orders,
                #                                                             self.numbers_of_variables,
                #                                                             self.boundary_interfaces_discretized)
                # estimators_increase_next,flags_increase_next = self.pde_type.compute_breakdown_criteria_increase(
                #                                                             predicted_values,
                #                                                             # padded_vectors_predicted_left,
                #                                                             # padded_vectors_predicted_right,
                #                                                             self.mesh.resolution,
                #                                                             delta_x,
                #                                                             delta_t,
                #                                                             self.max_order,
                #                                                             self.orders_cellwise,
                #                                                             self.orders,
                #                                                             self.numbers_of_variables,
                #                                                             self.boundary_interfaces_discretized)

                estimators_increase_prev,flags_increase_prev = self.pde_type.compute_refinement_criterion(
                                                                    prev_values,
                                                                    self.orders,
                                                                    self.max_order,
                                                                    self.numbers_of_variables,
                                                                    self.mesh.resolution,
                                                                    self.boundary_interfaces_discretized,
                                                                    delta_x)

                estimators_increase_next,flags_increase_next = self.pde_type.compute_refinement_criterion(
                                                                    predicted_values,
                                                                    self.orders,
                                                                    self.max_order,
                                                                    self.numbers_of_variables,
                                                                    self.mesh.resolution,
                                                                    self.boundary_interfaces_discretized,
                                                                    delta_x)

                flags_increase_max = np.maximum(flags_increase_next,flags_increase_prev)
                estimators_increase_max = np.maximum(estimators_increase_next,estimators_increase_prev)

                # estimators_decrease_prev,flags_decrease_prev = self.pde_type.compute_breakdown_criteria_decrease(
                #                                                             prev_values,
                #                                                             # padded_vectors_left,
                #                                                             # padded_vectors_right,
                #                                                             self.mesh.resolution,
                #                                                             delta_x,
                #                                                             delta_t,
                #                                                             self.max_order,
                #                                                             self.orders_cellwise,
                #                                                             self.orders,
                #                                                             self.numbers_of_variables,
                #                                                             self.boundary_interfaces_discretized,
                #                                                             flags_increase_max)

                # estimators_decrease_next,flags_decrease_next = self.pde_type.compute_breakdown_criteria_decrease(
                #                                                             predicted_values,
                #                                                             # padded_vectors_predicted_left,
                #                                                             # padded_vectors_predicted_right,
                #                                                             self.mesh.resolution,
                #                                                             delta_x,
                #                                                             delta_t,
                #                                                             self.max_order,
                #                                                             self.orders_cellwise,
                #                                                             self.orders,
                #                                                             self.numbers_of_variables,
                #                                                             self.boundary_interfaces_discretized,
                #                                                             flags_increase_max)

                estimators_decrease_prev,flags_decrease_prev = self.pde_type.compute_coarsening_criterion(
                                                                    prev_values,
                                                                    self.orders,
                                                                    self.max_order,
                                                                    self.numbers_of_variables,
                                                                    self.mesh.resolution,
                                                                    self.boundary_interfaces_discretized,
                                                                    delta_x,
                                                                    flags_increase_max)

                estimators_decrease_next,flags_decrease_next = self.pde_type.compute_coarsening_criterion(
                                                                    predicted_values,
                                                                    self.orders,
                                                                    self.max_order,
                                                                    self.numbers_of_variables,
                                                                    self.mesh.resolution,
                                                                    self.boundary_interfaces_discretized,
                                                                    delta_x,
                                                                    flags_increase_max)

                flags_decrease_min = np.minimum(flags_decrease_next,flags_decrease_prev)
                estimators_decrease_max = np.maximum(estimators_decrease_next,estimators_decrease_prev)

                self.breakdown_estimators = np.stack((estimators_decrease_max, estimators_increase_max),axis=1)
                self.breakdown_criteria_flags = self.pde_type.decompose_domain(self.mesh.resolution,
                                                                                self.max_order,
                                                                                self.orders_cellwise,
                                                                                flags_decrease_min,
                                                                                flags_increase_max)
                
                values = self._reconstruct_subdomains(values,delta_x,delta_t,not self.path_conservative_interface_coupling)
                padded_vectors_left,padded_vectors_right = self._construct_padded_vectors(values)

            values = self._transport_step(values,
                                          padded_vectors_left,
                                          padded_vectors_right,
                                          delta_t,
                                          delta_x,
                                          self.spatial_discretization,
                                          not self.path_conservative_interface_coupling)
            values = self._source_step(values,
                                       padded_vectors_left,
                                       padded_vectors_right,
                                       delta_t,
                                       delta_x,
                                       not self.path_conservative_interface_coupling)
            
            step_count += 1

            print()
            print('time: '+str(t))
            print('step size: '+str(delta_t))
            print()

            t+=delta_t
        print(self.boundary_interfaces_discretized)
        values = simulation_data = self._post_processing(values)
        return simulation_data  

    def _transport_step(self,
                        values,
                        padded_vectors_left,
                        padded_vectors_right,
                        delta_t,
                        delta_x,
                        spatial_discretization,
                        padding) -> np.ndarray:
        
        fluctuations_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        fluctuations_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))

        prev_values = np.copy(values)

        right_boundary_subdomain = 0

        for m in range(len(self.boundary_interfaces_discretized)):
            order = self.orders[m]
            n_variables = self.numbers_of_variables[m]
            def system_matrix(cell_values):
                return self.pde_type.compute_system_matrix(order,cell_values)               

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = self.boundary_interfaces_discretized[m]                 

            if padding:
                # if order > self.orders[m+1]:
                #     prev_values[right_boundary_subdomain+1,:] = padded_vectors_right[m]
                #     prev_values[right_boundary_subdomain,:] = values[right_boundary_subdomain,:]
                # else:
                #     prev_values[right_boundary_subdomain,:] = padded_vectors_left[m]
                #     prev_values[right_boundary_subdomain+1,:] = values[right_boundary_subdomain+1,:]
                if order < self.orders[m+1]:
                    prev_values[right_boundary_subdomain,:] = padded_vectors_left[m,:]
                else:
                    prev_values[right_boundary_subdomain+1,:] = padded_vectors_right[m,:]    

            for i in range(left_boundary_subdomain,right_boundary_subdomain):
                fluctuations_min[i-1,:n_variables], fluctuations_plus[i-1,:n_variables] =\
                    spatial_discretization.compute_fluctuation(
                                                        prev_values[i-1,:n_variables],
                                                        prev_values[i,:n_variables],
                                                        system_matrix,
                                                        delta_t,
                                                        delta_x)   
            i = right_boundary_subdomain
            fluctuations_min[i-1,:n_variables], fluctuations_plus[i-1,:n_variables] =\
                self.spatial_discretization_interface.compute_fluctuation(
                                                    prev_values[i-1,:n_variables],
                                                    prev_values[i,:n_variables],
                                                    system_matrix,
                                                    delta_t,
                                                    delta_x) 
        order = self.orders[-1]
        n_variables = self.numbers_of_variables[-1]

        # if padding:
        #     if order < self.orders[-2]:
        #         prev_values[left_boundary_subdomain,:] = padded_vectors_right[-1]
        #         prev_values[left_boundary_subdomain-1,:] = values[right_boundary_subdomain-1,:]
        #     else:
        #         prev_values[left_boundary_subdomain-1,:] = padded_vectors_left[-1]
        #         prev_values[left_boundary_subdomain,:] = values[left_boundary_subdomain,:]

        def system_matrix(cell_values):
            return self.pde_type.compute_system_matrix(order,cell_values)
        
        for i in range(right_boundary_subdomain+1,self.mesh.resolution+2):
                fluctuations_min[i-1,:n_variables], fluctuations_plus[i-1,:n_variables] =\
                    spatial_discretization.compute_fluctuation(
                                                        prev_values[i-1,:n_variables],
                                                        prev_values[i,:n_variables],
                                                        system_matrix,
                                                        delta_t,
                                                        delta_x) 

        right_boundary_subdomain = 0

        for m in range(len(self.boundary_interfaces_discretized)):
            n_variables = self.numbers_of_variables[m]

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = self.boundary_interfaces_discretized[m]
            
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                values[i,:n_variables] = prev_values[i,:n_variables]\
                    - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables]+fluctuations_min[i,:n_variables]) 

        n_variables = self.numbers_of_variables[-1]

        for i in range(right_boundary_subdomain+1,self.mesh.resolution+1):
            values[i,:n_variables] = prev_values[i,:n_variables]\
                - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables]+fluctuations_min[i,:n_variables]) 

        return values      

    def _source_step(self,
                    values,
                    padded_vectors_left,
                    padded_vectors_right,
                    delta_t,
                    delta_x,
                    padding) -> np.ndarray:

        right_boundary_subdomain = 0

        for m in range(len(self.boundary_interfaces_discretized)):
            order = self.orders[m]
            n_variables = self.numbers_of_variables[m]

            def source_term(cell_values,delta_t):
                return self.pde_type.compute_source_term(order,cell_values,delta_t)

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = self.boundary_interfaces_discretized[m]
            
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                values[i,:n_variables] = self.time_integration.integrate(values[i,:n_variables],source_term,delta_t)
            
            if padding:
                if order < self.orders[m+1]:
                    values[right_boundary_subdomain,n_variables:] = \
                        values[right_boundary_subdomain+1,n_variables:] # update boundary interface boundary condition 
                else:
                    values[right_boundary_subdomain+1,self.numbers_of_variables[m+1]:] = \
                        values[right_boundary_subdomain,self.numbers_of_variables[m+1]:] # update boundary interface boundary condition 

        order = self.orders[-1]
        n_variables = self.numbers_of_variables[-1]

        def source_term(cell_values,delta_t):
            return self.pde_type.compute_source_term(order,cell_values,delta_t)         
                    
        for i in range(right_boundary_subdomain+1,self.mesh.resolution+1): 
            values[i,:n_variables] = self.time_integration.integrate(values[i,:n_variables],source_term,delta_t)

        if padding:
            if order < self.orders[m-1]:
                values[right_boundary_subdomain+1,n_variables:] = \
                    values[right_boundary_subdomain,n_variables:] # update boundary interface boundary condition 
            else:
                values[right_boundary_subdomain,self.numbers_of_variables[m-1]:] = \
                    values[right_boundary_subdomain+1,self.numbers_of_variables[m-1]:] # update boundary interface boundary condition 

        return values 

    def _transport_step_augmented(self,
                                values,
                                delta_t,
                                delta_x,
                                spatial_discretization) -> np.ndarray:
        
        # nr_augmented_variables = 1 #Generalize this such that it includes the general case

        # fluctuations_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        # fluctuations_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))

        # prev_values = np.copy(values)

        # right_boundary_subdomain = 0

        # for m in range(len(self.boundary_interfaces_discretized)):
        #     order = self.orders[m]
        #     n_variables = self.numbers_of_variables[m]

        #     left_boundary_subdomain = right_boundary_subdomain+1
        #     right_boundary_subdomain = self.boundary_interfaces_discretized[m]                 

        #     if order < self.max_order:
        #         nr_augmented_variables = 1
        #         def system_matrix(cell_values):
        #             # return self.pde_type.compute_system_matrix_augmented(order,cell_values)
        #             return self.pde_type.compute_system_matrix(order+1,cell_values)
        #     else:
        #         def system_matrix(cell_values):
        #             return self.pde_type.compute_system_matrix(order,cell_values)
        #         nr_augmented_variables = 0                                             

        #     if padding:
        #         # if m > 0 and order < self.orders[m-1]:
        #         #     prev_values[left_boundary_subdomain,:] = padded_vectors_right[m-1]
        #         #     prev_values[left_boundary_subdomain-1,:] = values[right_boundary_subdomain-1,:]
        #         # else:
        #         #     prev_values[left_boundary_subdomain-1,:] = padded_vectors_left[m-1]
        #         #     prev_values[left_boundary_subdomain,:] = values[left_boundary_subdomain,:]
        #         # if order < self.orders[m+1]:
        #         #     prev_values[right_boundary_subdomain,:] = padded_vectors_left[m,:]
        #         # else:
        #         #     prev_values[right_boundary_subdomain+1,:] = padded_vectors_right[m,:]    

        #         if order < self.orders[m+1]:
        #             prev_values[right_boundary_subdomain,:] =\
        #                 padded_vectors_left[m,:]
        #         else:
        #             prev_values[right_boundary_subdomain+1,:] =\
        #                 padded_vectors_right[m,:]

        #     for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
        #         fluctuations_min[i-1,:n_variables+nr_augmented_variables],\
        #             fluctuations_plus[i-1,:n_variables+nr_augmented_variables] =\
        #             spatial_discretization.compute_fluctuation(
        #                                                 prev_values[i-1,:n_variables+nr_augmented_variables],
        #                                                 prev_values[i,:n_variables+nr_augmented_variables],
        #                                                 system_matrix,
        #                                                 delta_t,
        #                                                 delta_x)   

        # order = self.orders[-1]
        # n_variables = self.numbers_of_variables[-1]

        # # if padding:
        # #     if order < self.orders[-2]:
        # #         prev_values[left_boundary_subdomain,:] = padded_vectors_right[-1]
        # #         prev_values[left_boundary_subdomain-1,:] = values[right_boundary_subdomain-1,:]
        # #     else:
        # #         prev_values[left_boundary_subdomain-1,:] = padded_vectors_left[-1]
        # #         prev_values[left_boundary_subdomain,:] = values[left_boundary_subdomain,:]

        # if order < self.max_order:
        #     nr_augmented_variables = 1
        #     def system_matrix(cell_values):
        #         return self.pde_type.compute_system_matrix_augmented(order,cell_values)
        # else:
        #     def system_matrix(cell_values):
        #         return self.pde_type.compute_system_matrix(order,cell_values)
        #     nr_augmented_variables = 0   
        
        # for i in range(right_boundary_subdomain+1,self.mesh.resolution+2):
        #         fluctuations_min[i-1,:n_variables+nr_augmented_variables],\
        #             fluctuations_plus[i-1,:n_variables+nr_augmented_variables] =\
        #             spatial_discretization.compute_fluctuation(
        #                                                 prev_values[i-1,:n_variables+nr_augmented_variables],
        #                                                 prev_values[i,:n_variables+nr_augmented_variables],
        #                                                 system_matrix,
        #                                                 delta_t,
        #                                                 delta_x) 

        # right_boundary_subdomain = 0

        # for m in range(len(self.boundary_interfaces_discretized)):
        #     n_variables = self.numbers_of_variables[m]

        #     if n_variables == self.max_number_of_variables:
        #         nr_augmented_variables = 0   
        #     else:
        #         nr_augmented_variables = 1

        #     left_boundary_subdomain = right_boundary_subdomain+1
        #     right_boundary_subdomain = self.boundary_interfaces_discretized[m]
            
        #     for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
        #         values[i,:n_variables+nr_augmented_variables] = values[i,:n_variables+nr_augmented_variables]\
        #             - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables+nr_augmented_variables]+\
        #                                fluctuations_min[i,:n_variables+nr_augmented_variables]) 

        # n_variables = self.numbers_of_variables[-1]
        # if n_variables == self.max_number_of_variables:
        #     nr_augmented_variables = 0   
        # else:
        #     nr_augmented_variables = 1

        # for i in range(right_boundary_subdomain+1,self.mesh.resolution+1):
        #     values[i,:n_variables+nr_augmented_variables] = values[i,:n_variables+nr_augmented_variables]\
        #         - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables+nr_augmented_variables]+\
        #                            fluctuations_min[i,:n_variables+nr_augmented_variables]) 

        # return values      

        nr_augmented_variables = 1 #Generalize this such that it includes the general case

        fluctuations_min = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))
        fluctuations_plus = np.zeros((self.mesh.resolution+1,self.max_number_of_variables))

        prev_values = np.copy(values)
        predicted_values_transport = np.copy(values)

        right_boundary_subdomain = 0
        for m in range(len(self.boundary_interfaces_discretized)):
            order = self.orders[m]
            n_variables = self.numbers_of_variables[m]

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = self.boundary_interfaces_discretized[m]                 

            if order < self.max_order:
                nr_augmented_variables = 1
                def system_matrix(cell_values):
                    return self.pde_type.compute_system_matrix_augmented(order,cell_values)
                    # return self.pde_type.compute_system_matrix(order+1,cell_values)
            else:
                def system_matrix(cell_values):
                    return self.pde_type.compute_system_matrix(order,cell_values)
                nr_augmented_variables = 0                                             

            left_boundary_value = np.zeros(self.max_number_of_variables)
            right_boundary_value = np.zeros(self.max_number_of_variables)
            left_boundary_value[:n_variables] = prev_values[left_boundary_subdomain-1,:n_variables]
            right_boundary_value[:n_variables] = prev_values[right_boundary_subdomain+1,:n_variables]
            # fluctuations_plus[left_boundary_subdomain-1,:] = np.zeros(self.max_number_of_variables)
            # # if m == 0 and self.orders[m] < self.max_order:
            # #     print(np.max(values[left_boundary_subdomain:right_boundary_subdomain+1,self.numbers_of_variables[m]-1]))
            # # if m == 0 and self.orders[m] < self.max_order:
            # #     print(values[left_boundary_subdomain:right_boundary_subdomain+1,self.numbers_of_variables[m]])
            y,fluctuations_plus[left_boundary_subdomain-1,:n_variables+nr_augmented_variables] =\
                spatial_discretization.compute_fluctuation(
                                                    left_boundary_value[:n_variables+nr_augmented_variables],
                                                    values[left_boundary_subdomain,:n_variables+nr_augmented_variables],
                                                    system_matrix,
                                                    delta_t,
                                                    delta_x)  
            for i in range(left_boundary_subdomain+1,right_boundary_subdomain+1):
                fluctuations_min[i-1,:n_variables+nr_augmented_variables],\
                    fluctuations_plus[i-1,:n_variables+nr_augmented_variables] =\
                    spatial_discretization.compute_fluctuation(
                                                        values[i-1,:n_variables+nr_augmented_variables],
                                                        values[i,:n_variables+nr_augmented_variables],
                                                        system_matrix,
                                                        delta_t,
                                                        delta_x)   
                predicted_values_transport[i-1,:n_variables+nr_augmented_variables] =\
                    values[i-1,:n_variables+nr_augmented_variables]\
                    - delta_t/delta_x*(fluctuations_plus[i-2,:n_variables+nr_augmented_variables]+\
                                       fluctuations_min[i-1,:n_variables+nr_augmented_variables]) 
            i = right_boundary_subdomain
            fluctuations_min[right_boundary_subdomain,:n_variables+nr_augmented_variables],y =\
                spatial_discretization.compute_fluctuation(
                                                    values[i,:n_variables+nr_augmented_variables],
                                                    right_boundary_value[:n_variables+nr_augmented_variables],
                                                    system_matrix,
                                                    delta_t,
                                                    delta_x)              
            predicted_values_transport[i,:n_variables+nr_augmented_variables] =\
                values[i,:n_variables+nr_augmented_variables]\
                - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables+nr_augmented_variables]+\
                                   fluctuations_min[i,:n_variables+nr_augmented_variables]) 
            # if m == 0 and self.orders[m] < self.max_order:
            #     print(values[left_boundary_subdomain:right_boundary_subdomain+1,self.numbers_of_variables[m]])

        order = self.orders[-1]
        n_variables = self.numbers_of_variables[-1]

        if order < self.max_order:
            nr_augmented_variables = 1
            def system_matrix(cell_values):
                return self.pde_type.compute_system_matrix_augmented(order,cell_values)
        else:
            def system_matrix(cell_values):
                return self.pde_type.compute_system_matrix(order,cell_values)
            nr_augmented_variables = 0   

        left_boundary_subdomain = right_boundary_subdomain+1 
        left_boundary_value = np.zeros(self.max_number_of_variables)
        left_boundary_value[:n_variables] = prev_values[left_boundary_subdomain-1,:n_variables]
        y,fluctuations_plus[left_boundary_subdomain-1,:n_variables+nr_augmented_variables] =\
            spatial_discretization.compute_fluctuation(
                                                left_boundary_value[:n_variables+nr_augmented_variables],
                                                values[left_boundary_subdomain,:n_variables+nr_augmented_variables],
                                                system_matrix,
                                                delta_t,
                                                delta_x)  
        # fluctuations_plus[left_boundary_subdomain-1,:] = np.zeros(self.max_number_of_variables)
        for i in range(left_boundary_subdomain+1,self.mesh.resolution+2):
            fluctuations_min[i-1,:n_variables+nr_augmented_variables],\
                fluctuations_plus[i-1,:n_variables+nr_augmented_variables] =\
                spatial_discretization.compute_fluctuation(
                                                    values[i-1,:n_variables+nr_augmented_variables],
                                                    values[i,:n_variables+nr_augmented_variables],
                                                    system_matrix,
                                                    delta_t,
                                                    delta_x)   
            predicted_values_transport[i-1,:n_variables+nr_augmented_variables] =\
                values[i-1,:n_variables+nr_augmented_variables]\
                - delta_t/delta_x*(fluctuations_plus[i-2,:n_variables+nr_augmented_variables]+\
                                    fluctuations_min[i-1,:n_variables+nr_augmented_variables]) 

        return predicted_values_transport 

    def _source_step_augmented(self,
                    values,
                    delta_t,
                    delta_x) -> np.ndarray:

        predicted_values = np.copy(values)

        nr_augmented_variables = 1 #Generalize this such that it includes the general case

        right_boundary_subdomain = 0

        for m in range(len(self.boundary_interfaces_discretized)):
            order = self.orders[m]
            n_variables = self.numbers_of_variables[m]

            if order < self.max_order:
                nr_augmented_variables = 1
                def source_term(cell_values,delta_t):
                    return self.pde_type.compute_source_term_augmented(order,cell_values,delta_t)
            else:
                def source_term(cell_values,delta_t):
                    return self.pde_type.compute_source_term(order,cell_values,delta_t)        
                nr_augmented_variables = 0       

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = self.boundary_interfaces_discretized[m]

            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                predicted_values[i,:n_variables+nr_augmented_variables] \
                    = self.time_integration.integrate(values[i,:n_variables+nr_augmented_variables],source_term,delta_t)
            
            # if padding:
            #     if order < self.orders[m+1]:
            #         values[right_boundary_subdomain,n_variables:] = \
            #             values[right_boundary_subdomain+1,n_variables:] # update boundary interface boundary condition 
            #     else:
            #         values[right_boundary_subdomain+1,self.numbers_of_variables[m+1]:] = \
            #             values[right_boundary_subdomain,self.numbers_of_variables[m+1]:] # update boundary interface boundary condition 

        order = self.orders[-1]
        n_variables = self.numbers_of_variables[-1]

        if order < self.max_order:
            nr_augmented_variables = 1
            def source_term(cell_values,delta_t):
                return self.pde_type.compute_source_term_augmented(order,cell_values,delta_t)
        else:
            def source_term(cell_values,delta_t):
                return self.pde_type.compute_source_term(order,cell_values,delta_t)        
            nr_augmented_variables = 0         
                    
        for i in range(right_boundary_subdomain+1,self.mesh.resolution+1): 
            predicted_values[i,:n_variables+nr_augmented_variables] \
                = self.time_integration.integrate(values[i,:n_variables+nr_augmented_variables],source_term,delta_t)

        # for m in range(len(self.boundary_interfaces_discretized)):
        #     if self.orders[m] < self.orders[m+1]:
        #         values[self.boundary_interfaces_discretized[m],self.numbers_of_variables[m]:] =\
        #             values[self.boundary_interfaces_discretized[m]-1,self.numbers_of_variables[m]:]
        #     else:
        #         values[self.boundary_interfaces_discretized[m]+1,self.numbers_of_variables[m+1]:] =\
        #             values[self.boundary_interfaces_discretized[m]+2,self.numbers_of_variables[m+1]:]
            # if self.orders[m+1] > self.orders[m]:
            #     values[self.boundary_interfaces_discretized[m],self.numbers_of_variables[m]:] =\
            #         values[self.boundary_interfaces_discretized[m]+1,self.numbers_of_variables[m]:]
            # else:
            #     values[self.boundary_interfaces_discretized[m]+1,self.numbers_of_variables[m+1]:] =\
            #         values[self.boundary_interfaces_discretized[m],self.numbers_of_variables[m+1]:]
        return predicted_values 

    def _reconstruct_subdomains(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float,
                               padding) -> np.ndarray:
        

        for i in range(self.mesh.resolution):
            self.orders_cellwise[i+1] += int(self.breakdown_criteria_flags[i])
            self.numbers_of_variables_cellwise[i+1] += int(self.breakdown_criteria_flags[i])     
        self.orders_cellwise[0] = self.orders_cellwise[1]
        self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]
        self.orders_cellwise[-1] = self.orders_cellwise[-2]
        self.numbers_of_variables_cellwise[-1] = self.numbers_of_variables_cellwise[-2]

        if self.boundary_condition == 'PERIODIC':
            self.orders_cellwise[1] = max(self.orders_cellwise[1],self.orders_cellwise[-2])
            self.numbers_of_variables_cellwise[1] = self.pde_type.compute_number_of_variables(self.orders_cellwise[1])
            self.orders_cellwise[-2] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[-2] = self.numbers_of_variables_cellwise[1]
            self.orders_cellwise[0] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]
            self.orders_cellwise[-1] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[-1] = self.numbers_of_variables_cellwise[1]

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

        self.orders = orders_out
        self.boundary_interfaces_discretized = boundary_interfaces
        self.numbers_of_variables = number_of_variables_out

        # Set undefined moments and padded moments to zero
        right_boundary = -1
        for m in range(len(self.boundary_interfaces_discretized)):
            left_boundary = right_boundary+1
            right_boundary = self.boundary_interfaces_discretized[m]
            values[left_boundary:right_boundary+1,self.numbers_of_variables[m]:self.max_number_of_variables]=0
            if padding:
                if self.orders[m] < self.orders[m+1]:
                    values[right_boundary,self.numbers_of_variables[m]:] = \
                        values[right_boundary+1,self.numbers_of_variables[m]:] # update boundary interface boundary condition 
                else:
                    values[right_boundary+1,self.numbers_of_variables[m+1]:] = \
                        values[right_boundary,self.numbers_of_variables[m+1]:] # update boundary interface boundary condition 
        values[right_boundary+1:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables]=0 
                
        if len(self.boundary_interfaces_discretized) == 0:
            self.orders.append(self.orders[0])
            self.numbers_of_variables.append(self.numbers_of_variables[0])
            self.boundary_interfaces_discretized.append(np.floor_divide(self.mesh.resolution,2))
        
        if padding and self.orders[-2] > self.orders[-1]:
                values[right_boundary+1,self.numbers_of_variables[-1]:] = \
                    values[right_boundary,self.numbers_of_variables[-1]:] # update boundary interface boundary condition 

        return values

    def _construct_padded_vectors(self,values):

        padded_vectors_left = np.zeros((len(self.boundary_interfaces_discretized),self.max_number_of_variables))
        padded_vectors_right = np.zeros((len(self.boundary_interfaces_discretized),self.max_number_of_variables))        

        for i in range(len(self.boundary_interfaces_discretized)):
            if self.orders[i+1] > self.orders[i]:
                padded_vectors_left[i,:] = values[self.boundary_interfaces_discretized[i]+1,:]
                padded_vectors_left[i,:self.numbers_of_variables[i]] =\
                    values[self.boundary_interfaces_discretized[i],:self.numbers_of_variables[i]]
            else:
                padded_vectors_right[i,:] = values[self.boundary_interfaces_discretized[i],:]
                padded_vectors_right[i,:self.numbers_of_variables[i+1]] =\
                    values[self.boundary_interfaces_discretized[i]+1,:self.numbers_of_variables[i+1]]

        return padded_vectors_left, padded_vectors_right

    def get_breakdown_estimators(self):
        
        return self.breakdown_estimators[:,0], self.breakdown_estimators[:,1]

class SmoothedSubdomainReconstructionNew(SpatiallyAdaptiveSimulation1D):

    """
    This class represents a smoothed adaptive simulation in 1D. 
    It smooths the domain decomposition by grouping cells together and giving them the same model order.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell
    smooth_par : float
        the size of the subregions, relative to the size of the entire domain
    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

    Implemented methods from abstract class SpatiallyAdaptiveSimulation1D
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     

    Overriden methods from abstract class SpatiallAdaptiveSimulation1D
    def __init__(self,start_order,pde_type,mesh.RectangularMesh,boundary_condition,initial_condition,
                    breakdown_criterion,spatial_discretization,time_integration):
        initializes the SmoothedSubdomainReconstruction object
    
    Instance methods
    ----------------
    None
    """

    def __init__(self,
                 start_order: list,
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 smooth_par: int,
                 spatial_discretization_interior: spatialDiscretization.SpatialDiscretization,
                 spatial_discretization_interface: spatialDiscretization.SpatialDiscretization,
                 spatial_discretization_predictor: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration):

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

        self.orders = [start_order,start_order,start_order]
        self.pde_type = pde_type
        self.numbers_of_variables = [pde_type.compute_number_of_variables(start_order),
                                     pde_type.compute_number_of_variables(start_order),
                                     pde_type.compute_number_of_variables(start_order)]
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.breakdown_criterion = breakdown_criterion
        self.spatial_discretization = spatial_discretization_interior
        self.spatial_discretization_predictor = spatial_discretization_predictor
        self.spatial_discretization_interface = spatial_discretization_interface
        self.path_conservative_interface_coupling = True
        self.time_integration = time_integration

        self.boundary_interfaces_discretized = [np.floor_divide(2*self.mesh.resolution,3),np.floor_divide(self.mesh.resolution,3)]
        self.max_order = start_order
        self.max_number_of_variables = pde_type.compute_number_of_variables(self.max_order)

        self.orders_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        self.numbers_of_variables_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        for i in range(self.mesh.resolution+2):
            self.orders_cellwise[i] = start_order
            self.numbers_of_variables_cellwise[i] = self.pde_type.compute_number_of_variables(start_order)

        self.dom_decomp_val_res1 = np.zeros(self.mesh.resolution)
        self.dom_decomp_val_res2 = np.zeros(self.mesh.resolution)

        self.breakdown_estimators = np.zeros((self.mesh.resolution,self.max_order+4))
        self.breakdown_criteria_flags = np.full(shape=self.mesh.resolution,dtype=int,fill_value=0)

        self.smooth_par = smooth_par
        self.orders_subdomains = np.full(shape=self.smooth_par,fill_value=start_order,dtype=int)
        self.numbers_of_variables_subdomains \
            = np.full(shape=self.smooth_par,fill_value=self.pde_type.compute_number_of_variables(start_order),dtype=int)

        self.n_cells_subdomain = int(np.floor(self.mesh.resolution/self.smooth_par))
        self.subdomain_start = int(np.ceil(self.n_cells_subdomain/2))

    def _reconstruct_subdomains(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float,
                               padding) -> np.ndarray:

        for i in range(self.mesh.resolution):
            self.orders_cellwise[i+1] += int(self.breakdown_criteria_flags[i])
            self.numbers_of_variables_cellwise[i+1] += int(self.breakdown_criteria_flags[i])   
        self.orders_cellwise[0] = self.orders_cellwise[1]
        self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]
        self.orders_cellwise[-1] = self.orders_cellwise[-2]
        self.numbers_of_variables_cellwise[-1] = self.numbers_of_variables_cellwise[-2]

        if self.boundary_condition == 'PERIODIC':
            self.orders_cellwise[1] = max(self.orders_cellwise[1],self.orders_cellwise[-2])
            self.numbers_of_variables_cellwise[1] = self.pde_type.compute_number_of_variables(self.orders_cellwise[1])
            self.orders_cellwise[-2] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[-2] = self.numbers_of_variables_cellwise[1]
            self.orders_cellwise[0] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]
            self.orders_cellwise[-1] = self.orders_cellwise[1]
            self.numbers_of_variables_cellwise[-1] = self.numbers_of_variables_cellwise[1]

        boundary_interfaces = []
        orders_out = []
        number_of_variables_out = []

        values_copy = np.copy(values)

        # self._update_domain_decomposition_pointwise(values,delta_x,delta_t)

        k = 1
        for i in range(self.subdomain_start,self.mesh.resolution-self.n_cells_subdomain,self.n_cells_subdomain):
            local_order = np.max(self.orders_cellwise[i:i+self.n_cells_subdomain])
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            self.orders_subdomains[k] = local_order
            self.numbers_of_variables_subdomains[k] = local_number_of_variables
            self.orders_cellwise[i:i+self.n_cells_subdomain] = local_order
            self.numbers_of_variables_cellwise[i:i+self.n_cells_subdomain] = local_number_of_variables
            k += 1

        if self.boundary_condition == 'PERIODIC':
            local_order = max(np.max(self.orders_cellwise[1:self.subdomain_start]),np.max(self.orders_cellwise[i+self.smooth_par:-1]))
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            self.orders_subdomains[0] = local_order
            self.numbers_of_variables_subdomains[0] = local_number_of_variables
            self.orders_cellwise[:self.subdomain_start] = local_order
            self.orders_cellwise[i+self.smooth_par:] = local_order
            self.numbers_of_variables_cellwise[:self.subdomain_start] = local_number_of_variables
            self.numbers_of_variables_cellwise[i+self.smooth_par:] = local_number_of_variables

            orders_out.append(int(self.orders_subdomains[0]))
            number_of_variables_out.append(int(self.numbers_of_variables_subdomains[0]))
            for i in range(self.smooth_par-1):
                if self.orders_subdomains[i] != self.orders_subdomains[i+1]:
                    orders_out.append(int(self.orders_subdomains[i+1]))
                    number_of_variables_out.append(int(self.numbers_of_variables_subdomains[i+1]))
                    boundary_interfaces.append(self.subdomain_start+self.n_cells_subdomain*(i+1)-1)
            if self.orders_subdomains[-1] != self.orders_subdomains[0]:
                boundary_interfaces.append(self.subdomain_start+self.n_cells_subdomain*(i+1)-1)
            orders_out.append(int(self.orders_subdomains[0]))
            number_of_variables_out.append(int(self.numbers_of_variables_subdomains[0]))
        else:
            print(self.orders_subdomains)
            self.orders_cellwise[:self.subdomain_start] = self.orders_subdomains[1]
            self.orders_cellwise[i+self.smooth_par:] = self.orders_subdomains[-1]
            self.numbers_of_variables_cellwise[:self.subdomain_start] = self.numbers_of_variables_subdomains[1]
            self.numbers_of_variables_cellwise[i+self.smooth_par:] = self.numbers_of_variables_subdomains[-1]

            orders_out.append(int(self.orders_subdomains[1]))
            number_of_variables_out.append(int(self.numbers_of_variables_subdomains[1]))
            for i in range(1,self.smooth_par-1):
                if self.orders_subdomains[i] != self.orders_subdomains[i+1]:
                    orders_out.append(int(self.orders_subdomains[i+1]))
                    number_of_variables_out.append(int(self.numbers_of_variables_subdomains[i+1]))
                    boundary_interfaces.append(self.subdomain_start+self.n_cells_subdomain*i-1)

        self.orders = orders_out
        self.boundary_interfaces_discretized = boundary_interfaces
        self.numbers_of_variables = number_of_variables_out
        if len(self.boundary_interfaces_discretized) == 0:
            self.orders.append(self.orders_subdomains[1])
            self.orders.append(self.orders_subdomains[1])
            self.numbers_of_variables.append(int(self.numbers_of_variables_subdomains[1]))
            self.numbers_of_variables.append(int(self.numbers_of_variables_subdomains[1]))
            self.boundary_interfaces_discretized.append(np.floor_divide(self.mesh.resolution,3))
            self.boundary_interfaces_discretized.append(np.floor_divide(2*self.mesh.resolution,3))

        # Set undefined moments and padded moments to zero
        right_boundary = -1
        for m in range(len(self.boundary_interfaces_discretized)):
            left_boundary = right_boundary+1
            right_boundary = self.boundary_interfaces_discretized[m]
            values_copy[left_boundary:right_boundary+1,self.numbers_of_variables[m]:self.max_number_of_variables]=0
            if padding:
                if self.orders[m] < self.orders[m+1]:
                    values_copy[right_boundary,self.numbers_of_variables[m]:] = \
                        values[right_boundary+1,self.numbers_of_variables[m]:] # update boundary interface boundary condition 
                else:
                    values_copy[right_boundary+1,self.numbers_of_variables[m+1]:] = \
                        values[right_boundary,self.numbers_of_variables[m+1]:] # update boundary interface boundary condition 
        values_copy[right_boundary+1:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables]=0 

        return values_copy

class SmoothedSubdomainReconstructionWithInterpolation(SpatiallyAdaptiveSimulation1D):
    """
    This class represents an adaptive simulation in 1D that uses a path-conservative numerical scheme 
    (and spatial coupling), and an interpolation technique to fill in the moment values when the order is increased.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell
    smooth_par : float
        the size of the subregions, relative to the size of the entire domain
    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

    Implemented methods from abstract class SpatiallyAdaptiveSimulation1D
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     

    Overriden methods from abstract class SpatiallAdaptiveSimulation1D
    def __init__(self,start_order,pde_type,mesh.RectangularMesh,boundary_condition,initial_condition,
                    breakdown_criterion,spatial_discretization,time_integration):
        initializes the SmoothedSubdomainReconstruction object
    
    Instance methods
    ----------------
    None
    """

    def __init__(self,
                 start_order: list,
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 breakdown_criterion: str,
                 smooth_par: int,
                 spatial_discretization_interior: spatialDiscretization.SpatialDiscretization,
                 spatial_discretization_interface: spatialDiscretization.SpatialDiscretization,
                 spatial_discretization_predictor: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration):

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

        self.orders = [start_order,start_order,start_order]
        self.pde_type = pde_type
        self.numbers_of_variables = [pde_type.compute_number_of_variables(start_order),
                                     pde_type.compute_number_of_variables(start_order),
                                     pde_type.compute_number_of_variables(start_order)]
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.breakdown_criterion = breakdown_criterion
        self.spatial_discretization = spatial_discretization_interior
        self.spatial_discretization_predictor = spatial_discretization_predictor
        self.spatial_discretization_interface = spatial_discretization_interface
        self.path_conservative_interface_coupling = True
        self.time_integration = time_integration
        self.boundary_interfaces_discretized = [np.floor_divide(2*self.mesh.resolution,3),np.floor_divide(self.mesh.resolution,3)]

        self.max_order = start_order
        self.max_number_of_variables = pde_type.compute_number_of_variables(self.max_order)

        self.orders_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        self.numbers_of_variables_cellwise = np.zeros(self.mesh.resolution+2,dtype=int)
        for i in range(self.mesh.resolution+2):
            self.orders_cellwise[i] = start_order
            self.numbers_of_variables_cellwise[i] = self.pde_type.compute_number_of_variables(start_order)

        self.dom_decomp_val_res1 = np.zeros(self.mesh.resolution)
        self.dom_decomp_val_res2 = np.zeros(self.mesh.resolution)

        self.breakdown_estimators = np.zeros((self.mesh.resolution,self.max_order+4))

        self.smooth_par = smooth_par
        self.orders_subdomains = np.full(shape=self.smooth_par+1,fill_value=start_order,dtype=int)
        self.breakdown_criteria_flags_subdomains = np.full(shape=self.smooth_par+1,fill_value=0,dtype=int)
        self.breakdown_criteria_flags = np.full(shape=self.mesh.resolution,dtype=int,fill_value=0)
        self.numbers_of_variables_subdomains = np.full(shape=self.smooth_par+1,fill_value=self.pde_type.compute_number_of_variables(start_order),dtype=int)
        self.boundary_interfaces_subdomains = np.zeros(self.smooth_par,dtype=int)  

        self.n_cells_subdomain = int(np.floor(self.mesh.resolution/self.smooth_par))
        self.subdomain_start = int(np.ceil(self.n_cells_subdomain/2))+1
        self.boundary_interfaces_subdomains[0] = self.subdomain_start - 1
        for i in range(1,self.smooth_par):
            self.boundary_interfaces_subdomains[i] = self.subdomain_start + i*self.n_cells_subdomain - 1

        if self.boundary_condition != 'PERIODIC':
            self.increase_flags_subdomains = np.full(shape=self.smooth_par,fill_value=0,dtype=int)  
            self.orders_subdomains = np.full(shape=self.smooth_par,fill_value=start_order,dtype=int)
            self.increase_flags_subdomains = np.full(shape=self.smooth_par,fill_value=0,dtype=int)
            self.numbers_of_variables_subdomains = np.full(shape=self.smooth_par,fill_value=self.pde_type.compute_number_of_variables(start_order),dtype=int)
            self.boundary_interfaces_subdomains = np.zeros(self.smooth_par-1,dtype=int)
            for i in range(self.smooth_par-1):
                self.boundary_interfaces_subdomains[i] = (i+1)*self.n_cells_subdomain - 1 

        self._reconstruct_subdomains = self._reconstruct_subdomains()
        self._interpolate_subdomains = self._interpolate_subdomains()

    def _reconstruct_subdomains(self) -> np.ndarray:
        
        _reconstruct_subdomains_fun = self._reconstruct_subdomains_periodicBoundary if self.boundary_condition == 'PERIODIC'\
            else self._reconstruct_subdomains_nonPeriodicBoundary

        return _reconstruct_subdomains_fun

    def _reconstruct_subdomains_nonPeriodicBoundary(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float,
                               padding) -> np.ndarray:

        for i in range(self.mesh.resolution):
            self.orders_cellwise[i+1] += int(self.breakdown_criteria_flags[i])
            self.numbers_of_variables_cellwise[i+1] += int(self.breakdown_criteria_flags[i])   
        self.orders_cellwise[0] = self.orders_cellwise[1]
        self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]
        self.orders_cellwise[-1] = self.orders_cellwise[-2]
        self.numbers_of_variables_cellwise[-1] = self.numbers_of_variables_cellwise[-2]

        boundary_interfaces = []
        orders_out = []
        number_of_variables_out = []

        values_copy = np.copy(values)

        for i in range(self.smooth_par-1):
            local_order = np.max(self.orders_cellwise[i*self.n_cells_subdomain:(i+1)*self.n_cells_subdomain])
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            max_increase_flags = int(np.max(
                self.breakdown_criteria_flags[max(i*self.n_cells_subdomain-1,0):(i+1)*self.n_cells_subdomain-1]))
            self.increase_flags_subdomains[i] = max_increase_flags
            self.orders_subdomains[i] = local_order
            self.numbers_of_variables_subdomains[i] = local_number_of_variables
            self.orders_cellwise[i*self.n_cells_subdomain:(i+1)*self.n_cells_subdomain] = local_order
            self.numbers_of_variables_cellwise[i*self.n_cells_subdomain:(i+1)*self.n_cells_subdomain] = local_number_of_variables
        local_order = np.max(self.orders_cellwise[(self.smooth_par-1)*self.n_cells_subdomain:])
        self.orders_subdomains[-1] = local_order
        self.orders_cellwise[(self.smooth_par-1)*self.n_cells_subdomain:] = local_order
        self.numbers_of_variables_subdomains[-1] = self.pde_type.compute_number_of_variables(local_order)
        self.numbers_of_variables_cellwise[(self.smooth_par-1)*self.n_cells_subdomain:] = self.numbers_of_variables_subdomains[-1]
        max_increase_flags = int(np.max(
            self.breakdown_criteria_flags[(self.smooth_par-1)*self.n_cells_subdomain-1:]))
        self.increase_flags_subdomains[-1] = max_increase_flags

        orders_out.append(int(self.orders_subdomains[0]))
        number_of_variables_out.append(int(self.numbers_of_variables_subdomains[0]))
        for i in range(self.smooth_par-1):
            if self.orders_subdomains[i] != self.orders_subdomains[i+1]:
                orders_out.append(int(self.orders_subdomains[i+1]))
                number_of_variables_out.append(int(self.numbers_of_variables_subdomains[i+1]))
                boundary_interfaces.append(self.n_cells_subdomain*(i+1)-1)

        self.orders = orders_out
        self.boundary_interfaces_discretized = boundary_interfaces
        self.numbers_of_variables = number_of_variables_out
        if len(self.boundary_interfaces_discretized) == 0:
            self.orders.append(self.orders_subdomains[1])
            self.orders.append(self.orders_subdomains[1])
            self.numbers_of_variables.append(int(self.numbers_of_variables_subdomains[1]))
            self.numbers_of_variables.append(int(self.numbers_of_variables_subdomains[1]))
            self.boundary_interfaces_discretized.append(np.floor_divide(self.mesh.resolution,3))
            self.boundary_interfaces_discretized.append(np.floor_divide(2*self.mesh.resolution,3))

        # Set undefined moments and padded moments to zero
        right_boundary = -1
        for m in range(len(self.boundary_interfaces_discretized)):
            left_boundary = right_boundary+1
            right_boundary = self.boundary_interfaces_discretized[m]
            values_copy[left_boundary:right_boundary+1,self.numbers_of_variables[m]:self.max_number_of_variables]=0
            if padding:
                if self.orders[m] < self.orders[m+1]:
                    values_copy[right_boundary,self.numbers_of_variables[m]:] = \
                        values[right_boundary+1,self.numbers_of_variables[m]:] # update boundary interface boundary condition 
                else:
                    values_copy[right_boundary+1,self.numbers_of_variables[m+1]:] = \
                        values[right_boundary,self.numbers_of_variables[m+1]:] # update boundary interface boundary condition 
        values_copy[right_boundary+1:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables]=0 

        values_copy_interpolated = self._interpolate_subdomains(values_copy)
        print("Orders: "+str(self.orders_subdomains))
        print("increase flags: "+str(self.increase_flags_subdomains))        

        return values_copy_interpolated
        # return values_copy

    def _reconstruct_subdomains_periodicBoundary(self,
                               values: np.ndarray,
                               delta_x: float,
                               delta_t: float,
                               padding) -> np.ndarray:

        for i in range(self.mesh.resolution):
            self.orders_cellwise[i+1] += int(self.breakdown_criteria_flags[i])
            self.numbers_of_variables_cellwise[i+1] += int(self.breakdown_criteria_flags[i])   
        self.orders_cellwise[0] = self.orders_cellwise[1]
        self.numbers_of_variables_cellwise[0] = self.numbers_of_variables_cellwise[1]
        self.orders_cellwise[-1] = self.orders_cellwise[-2]
        self.numbers_of_variables_cellwise[-1] = self.numbers_of_variables_cellwise[-2]

        boundary_interfaces = []
        orders_out = []
        number_of_variables_out = []

        values_copy = np.copy(values)

        local_order_start = np.max(self.orders_cellwise[:self.subdomain_start])
        for i in range(self.smooth_par-1):
            local_order = np.max(self.orders_cellwise[self.subdomain_start+i*self.n_cells_subdomain:\
                                                      self.subdomain_start+(i+1)*self.n_cells_subdomain])
            local_number_of_variables = self.pde_type.compute_number_of_variables(local_order)
            self.orders_subdomains[i+1] = local_order
            self.numbers_of_variables_subdomains[i] = local_number_of_variables
            self.orders_cellwise[self.subdomain_start+i*self.n_cells_subdomain:\
                                    self.subdomain_start+(i+1)*self.n_cells_subdomain] = local_order
            self.numbers_of_variables_cellwise[self.subdomain_start+i*self.n_cells_subdomain:\
                                                    self.subdomain_start+(i+1)*self.n_cells_subdomain] = local_number_of_variables
        local_order_end = np.max(self.orders[self.subdomain_start+(self.smooth_par-1)*self.n_cells_subdomain:])
        max_boundary_order = max(local_order_start,local_order_end)
        self.orders_subdomains[0] = max_boundary_order
        self.orders_subdomains[-1] = max_boundary_order
        self.numbers_of_variables_subdomains[0] = self.pde_type.compute_number_of_variables(max_boundary_order)
        self.numbers_of_variables_subdomains[-1] = self.pde_type.compute_number_of_variables(max_boundary_order)

        orders_out.append(int(self.orders_subdomains[0]))
        number_of_variables_out.append(int(self.numbers_of_variables_subdomains[0]))
        for i in range(self.smooth_par):
            if self.orders_subdomains[i] != self.orders_subdomains[i+1]:
                orders_out.append(int(self.orders_subdomains[i+1]))
                number_of_variables_out.append(int(self.numbers_of_variables_subdomains[i+1]))
                boundary_interfaces.append(self.subdomain_start+self.n_cells_subdomain*i-1)

        self.orders = orders_out
        self.boundary_interfaces_discretized = boundary_interfaces
        self.numbers_of_variables = number_of_variables_out
        if len(self.boundary_interfaces_discretized) == 0:
            self.orders.append(self.orders_subdomains[1])
            self.orders.append(self.orders_subdomains[1])
            self.numbers_of_variables.append(int(self.numbers_of_variables_subdomains[1]))
            self.numbers_of_variables.append(int(self.numbers_of_variables_subdomains[1]))
            self.boundary_interfaces_discretized.append(np.floor_divide(self.mesh.resolution,3))
            self.boundary_interfaces_discretized.append(np.floor_divide(2*self.mesh.resolution,3))

        # Set undefined moments and padded moments to zero
        right_boundary = -1
        for m in range(len(self.boundary_interfaces_discretized)):
            left_boundary = right_boundary+1
            right_boundary = self.boundary_interfaces_discretized[m]
            values_copy[left_boundary:right_boundary+1,self.numbers_of_variables[m]:self.max_number_of_variables]=0
            if padding:
                if self.orders[m] < self.orders[m+1]:
                    values_copy[right_boundary,self.numbers_of_variables[m]:] = \
                        values[right_boundary+1,self.numbers_of_variables[m]:] # update boundary interface boundary condition 
                else:
                    values_copy[right_boundary+1,self.numbers_of_variables[m+1]:] = \
                        values[right_boundary,self.numbers_of_variables[m+1]:] # update boundary interface boundary condition 
        values_copy[right_boundary+1:self.mesh.resolution+2,self.numbers_of_variables[-1]:self.max_number_of_variables]=0 

        return values_copy

    def _process_domain_decomposition(self,values):

        # Set undefined moments to zero
        values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[0]:] = 0
        for i in range(1,len(self.boundary_interfaces)):
            values[self.boundary_interfaces[i-1]+1:self.boundary_interfaces[i]+1,self.numbers_of_variables[i]:] = 0
        values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[-1]:] = 0

        return values

    def _interpolate_subdomains(self) -> np.ndarray:
        
        _interpolate_subdomains_fun = self._interpolate_subdomains_periodicBoundary if self.boundary_condition == 'PERIODIC' \
            else self._interpolate_subdomains_nonPeriodicBoundary

        return _interpolate_subdomains_fun   

    def _interpolate_subdomains_interior(self,values):
        
        interpolated_values_interior = np.copy(values)
        for i in range(1,self.smooth_par-1):
            increase_flag = int(self.increase_flags_subdomains[i])
            if increase_flag > 0:
                boundary_left = self.boundary_interfaces_subdomains[i-1]+1
                boundary_right = self.boundary_interfaces_subdomains[i]+1
                nr_variables = self.numbers_of_variables_subdomains[i]
                interpolated_values_interior[boundary_left:boundary_right,nr_variables-increase_flag:nr_variables] =\
                    self._linear_interpolate(
                        values[boundary_left-1,nr_variables-increase_flag:nr_variables],
                        values[boundary_right,nr_variables-increase_flag:nr_variables],
                        1+self.n_cells_subdomain
                        )

        return interpolated_values_interior

    def _interpolate_subdomains_periodicBoundary(self,values):

        values = self._interpolate_subdomains_interior(values)

        if self.breakdown_criteria_flags_subdomains[0] > 0:
            interpol_values =\
                self._interpolate(np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[0]\
                                            -self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]],\
                                        values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[0]\
                                            -self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]]), axis=0),\
                        values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,\
                        self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:\
                            self.numbers_of_variables[0]],\
                        values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                        self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:\
                            self.numbers_of_variables[0]])
            values[self.boundary_interfaces[-1]+1:,\
                   self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]] =\
                   interpol_values[:self.mesh.resolution+1-self.boundary_interfaces[-1],:]
            values[:self.boundary_interfaces[0]+1,\
                   self.numbers_of_variables[0]-self.breakdown_criteria_flags_subdomains[0]:self.numbers_of_variables[0]] =\
                   interpol_values[self.mesh.resolution+1-self.boundary_interfaces[-1]:,:]            

        if self.smooth_par > 2:
            if self.breakdown_criteria_flags_subdomains[1] > 0:
                values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                    self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]] =\
                    self._interpolate(values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]],
                                      np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]],\
                                                     values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]]),axis=0),
                                      values[self.boundary_interfaces[1]+1:self.boundary_interfaces[2]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]])
                
            if self.breakdown_criteria_flags[-2] > 0:
                values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,\
                    self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:self.numbers_of_variables[-2]] =\
                    self._interpolate(values[self.boundary_interfaces[-2]+1:self.boundary_interfaces[-1]+1,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]],
                                      values[self.boundary_interfaces[-3]+1:self.boundary_interfaces[-2]+1,
                                                self.numbers_of_variables[-2]-self.breakdown_criteria_flags_subdomains[-2]:\
                                                    self.numbers_of_variables[-2]],
                                      np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[-2]\
                                                        -self.breakdown_criteria_flags_subdomains[-2]:self.numbers_of_variables[-2]],\
                                                     values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[-2]\
                                                        -self.breakdown_criteria_flags_subdomains[-2]:self.numbers_of_variables[-2]]),axis=0))
        else:
            if self.breakdown_criteria_flags_subdomains[1] > 0:
                values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,\
                    self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]] =\
                    self._interpolate(values[self.boundary_interfaces[0]+1:self.boundary_interfaces[1]+1,
                                                self.numbers_of_variables[1]-self.breakdown_criteria_flags_subdomains[1]:\
                                                    self.numbers_of_variables[1]],
                                      np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]],\
                                                     values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]]),axis=0),
                                      np.concatenate((values[self.boundary_interfaces[-1]+1:,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]],\
                                                     values[:self.boundary_interfaces[0]+1,self.numbers_of_variables[1]\
                                                        -self.breakdown_criteria_flags_subdomains[1]:self.numbers_of_variables[1]]),axis=0))

        return values            

    def _interpolate_subdomains_nonPeriodicBoundary(self,values):

        interpolated_values = self._interpolate_subdomains_interior(values)

        increase_flag = self.increase_flags_subdomains[0]
        boundary_right = self.boundary_interfaces_subdomains[0]+1
        nr_variables = self.numbers_of_variables_subdomains[0]
        if increase_flag > 0:
            interpolated_values[:boundary_right,nr_variables-increase_flag:nr_variables] =\
                    self._linear_interpolate(np.zeros(increase_flag),
                                             values[boundary_right,nr_variables-increase_flag:nr_variables],
                                             self.n_cells_subdomain+1)
        increase_flag = self.increase_flags_subdomains[-1]
        boundary_left = self.boundary_interfaces_subdomains[-1]+1
        nr_variables = self.numbers_of_variables_subdomains[-1]                                          
        if increase_flag > 0:
            interpolated_values[boundary_left:,nr_variables-increase_flag:nr_variables] =\
                    self._linear_interpolate(values[boundary_left-1,nr_variables-increase_flag:nr_variables],
                                             np.zeros(increase_flag),
                                             self.mesh.resolution-boundary_left+3)

        return interpolated_values                     

    def _linear_interpolate(self,value_left,value_right,delta_i): 

        slope = (value_right - value_left)/(delta_i)
        interpolated_values = value_left + np.outer(np.arange(1,delta_i), slope)
        # interpolated_values = np.zeros((delta_i-1,1))

        return interpolated_values

class SmoothedModelAdaptiveSimulation1D(SmoothedSubdomainReconstructionNew,ModelAdaptiveSimulation1D):

    """
    This class represents a smoothed adaptive simulation in 1D that uses the conservative interface flux coupling.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell
    smooth_par : float
        the size of the subregions, relative to the size of the entire domain
    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

        
    Inherited methods from parent class ConservativeAdaptiveSimulation
    -------------------------------------------------------------------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values

    Inherited methods from parent class SmoothedSubdomainReconstruction
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     
    def __init__(self,start_order,pde_type,mesh.RectangularMesh,boundary_condition,initial_condition,
                    breakdown_criterion,spatial_discretization,time_integration):
        initializes the SmoothedConsAdaptiveSimulation1D object
    
    Instance methods
    ----------------
    None
    """
    pass

class SmoothedModelAdaptiveSimulationWithInterpolation1D(SmoothedSubdomainReconstructionWithInterpolation,
                                                         ModelAdaptiveSimulation1D):

    """
    This class represents a smoothed adaptive simulation in 1D that uses the conservative interface flux coupling.

    ...

    Attributes
    ----------
    boundary_interfaces : list of floats
        list of the physical positions of the interfaces that separate the domain into subdomains
    orders : list of integers
        list of the order of the moment model in each subdomain
    max_order : int
        the maximum order in the simulation
    max_number_of_variables : int
        the maximum number of variables in the simulation
    orders_cellwise : list
        list conting the order in each cell
    number_of_variables_cellwise : list
        list containing the number of variables in each cell
    pde_type : str
        the partial differential equations that is simulated
    numbers_of_variables : int
        list of the number of state variables in each subdomain
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    breakdown_criterion : str
        breadown criterion for domain decomposition
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right hand side source term
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the values of the residual res1 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the values of the residual res2 (see definition paper) in each cell
        this value is only used for hierarchical moment equations
    breakdown_estimators : np.ndarray
        numpy array containg the values of each breakdown estimator in each grid cell
    smooth_par : float
        the size of the subregions, relative to the size of the entire domain
    
    Inherited methods from abstract class SpatiallyAdaptiveSimulation1D
    --------------------------------------------------------------------
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processed the end data of the simulation and prepares it for plotting
    def _update_domain_decomposition(self,values,tolerance_up,tolerance_down):
        updates the domain decompositions

        
    Inherited methods from parent class ConservativeAdaptiveSimulation
    -------------------------------------------------------------------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values

    Inherited methods from parent class SmoothedSubdomainReconstruction
    ----------------------------------------------------------------------
    def _reconstruct_subdomains(self,values,delta_x,delta_t):
        find boundary interfaces from the cellwise model orders     
    def __init__(self,start_order,pde_type,mesh.RectangularMesh,boundary_condition,initial_condition,
                    breakdown_criterion,spatial_discretization,time_integration):
        initializes the SmoothedConsAdaptiveSimulation1D object
    
    Instance methods
    ----------------
    None
    """
    pass