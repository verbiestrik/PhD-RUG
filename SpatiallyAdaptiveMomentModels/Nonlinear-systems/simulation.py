from abc import ABC, abstractmethod
import numpy as np
import pde
import mesh
import spatialDiscretization
import timeIntegration
from scipy.interpolate import BarycentricInterpolator

class Simulation(ABC):

    """
    This abstract class represents a simulation.

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
    time_integration : TimeIntegration
        the numerical method for the time integration of the source term
    CFL : float
        the CFL constant that is used to compute the time step size

    Instance methods
    ----------------    
    def __init__(self,pde_type,mesh,boundary_condition,initial_condition,time_integration,CFL_number):
        initializes the attributes that all simulation methods have in common
    
    Abstract methods
    ----------------
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
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 time_integration: timeIntegration.TimeIntegration,
                 CFL_number:float):
        """
        Constructs all the attributes that all simulation classes have in common.

        Parameters
        ----------
        pde_type : str
            the partial differential equations that is simulated
        mesh : RectangularMesh
            the used mesh
        boundary_condition: str
            the used boundary condition
        initial_condition: str
            the initial condition for the simulation
        time_integration: TimeIntegration
            the time integration method for the right-hand side source term

        """

        self.pde_type = pde_type
        self.mesh = mesh
        self.boundary_condition = boundary_condition
        self.initial_condition = initial_condition
        self.time_integration = time_integration
        self.CFL = CFL_number

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
        pass

    @abstractmethod
    def _update_boundary_conditions(self,
                                   values : np.ndarray,
                                   boundary : str):
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
    pde_type : PDE
        the partial differential equations that are simulated
    number_of_variables : int
        number of state variables
    mesh : RectangularMesh
        the used mesh
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    spatial_discretization : SpatialDiscretization
        the numerical method for the spatial discretization
    time_integration : TimeIntegration
        the time integration method for the right-hand side source term
    CFL : float
        the CFL constant that is used to compute the time step size

    Instance methods
    ----------------
    def __init__(self,order,pde_type,mesh,boundary_condition,initial_condition,spatial_discretization,time_integration,CFL_number):
        initializes the ClassicalSimulation1D object and initializes the attributes
    
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
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 spatial_discretization: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration,
                 CFL_number: float):
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
        
        super().__init__(pde_type,mesh,boundary_condition,initial_condition,time_integration,CFL_number)
        
        self.order = order
        self.number_of_variables = pde_type.compute_number_of_variables(self.order)
        self.spatial_discretization = spatial_discretization

    def run_simulation(self,
                       t_end: float,
                       g = 1) -> np.ndarray:

        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        values = self._get_initial_conditions(self.mesh.cell_center_positions)
        fluctuations_min = np.zeros((self.mesh.resolution+1,self.number_of_variables))
        fluctuations_plus = np.zeros((self.mesh.resolution+1,self.number_of_variables))

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

            delta_t = self.CFL*delta_x/max_speed 

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
            print()
            print('time: '+str(t))
            print('step size: '+str(delta_t))
            print()
            t += delta_t
            step += 1

        simulation_data = self._post_processing(values)
        return simulation_data

    def _get_initial_conditions(self,
                               cell_centers_x: np.ndarray) -> np.ndarray:
        
        initial_values = np.zeros((self.mesh.resolution+2,self.number_of_variables))

        for i in range(0,self.mesh.resolution):
            initial_values[i+1,:] = self.pde_type.get_initial_values(self.order,self.initial_condition,cell_centers_x[i])            
        
        return initial_values
    
    def _update_boundary_conditions(self,
                                   values: np.ndarray,
                                   boundary) -> np.ndarray:

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
                         values) -> np.ndarray:

        data_array = np.zeros((self.mesh.resolution,self.pde_type.compute_number_of_variables(self.order)+1)) # rewrite this such that it can be generalized to other PDE models

        for i in range(self.mesh.resolution):
            data_array[i,0] = self.mesh.cell_center_positions[i]

        data_array[:,1:] = values[1:-1,:]

        data_array[:,1:] = self.pde_type.convert_to_primitive(self.order,data_array[:,1:])

        return data_array 

class ModelAdaptiveMomentSimulation1D(Simulation):
    """
    This class represents a model-adaptive moment model simulation in 1D.

    ...

    Attributes
    ----------
    pde_type : PDE
        the partial differential equations that are simulated
    mesh : Mesh
        the discrete mesh on which the partial differential equations are numerically solved
    CFL_number : float
        the CFL constant that is used to compute the time step size
    boundary_condition : str
        the used boundary condition
    initial_condition : str
        the initial condition for the simulation
    time_integration : TimeIntegration
        the time integration method for the right-hand side source term
    spatial_discretization_interior : SpatialDiscretization
        the spatial discretization scheme for the interior of the subdomains
    spatial_discretization_predictor : SpatialDiscretization
        the spatial discretization scheme for the predictions that may be used in the domain decomposition criteria
    spatial_discretization_interface : SpatialDiscretization
        the spatial discretization scheme that is used at the boundary interfaces
    path_conservative_interface_coupling : bool
        true if a path-conservative interface coupling method is used, false if otherwise
    boundary_interfaces_discretized : list of floats
        list containing the cell indices of cells left of a boundary interface
    max_order : int
        maximum order of the moment model during the simulation
    min_order : int
        minimum order of the moment model during the simulation
    max_number_of_variables : int
        maximum number of variables in the moment model during the simulation
    min_number_of_variables : int
        minimum number of variables in the moment model during the simulation
    orders : list of integers
        list containing the orders in each subdomain
    numbers_of_variables : list of integers
        list containing the numbers of variables in each subdomain
    orders_cellwise : np.ndarray of integers
        numpy array containing the orders in each cell of the mesh
    numbers_of_variables_cellwise : np.ndarray of integers
        numpy array containing the numbers_of_variables in each cell of the mesh
    dom_decomp_val_res1 : np.ndarray of floats
        if the domain decomposition criteria are based on the higher-order numerical fluctuations,
        then these values in each mesh cell are stored in this numpy array
    dom_decomp_val_res2 : np.ndarray of floats
        if the domain decomposition criteria are based on the assumption that the highest-order moment equals zero,
        then the (approximate) time-evolutions of the highest-order moment in each mesh cell are stored 
        in this numpy array
    subdomainReconstruction : SubdomainReconstruction1D
        the object that is used to compute the new subdomains, boundary interface positions, and subdomain orders, 
        given the previously computed domain decomposition indicators
    type_model_error_estimator : str
        the type of model error estimator that is used for computing domain decomposition criteria
    domainDecomposition : DomainDecomposition1D
        the object that is used to compute the domain decomposition criteria and the indicators
    AdaptiveFVMStep1D : AdaptiveFVMStep1D
        the object that is used to advance the discrete state variables on the mesh in time
    number_of_breakdown_estimators_coarsening : int
        the number of breakdown estimators for model coarsening that is computed
    number_of_breakdown_estimators_refinement : int
        the number of breakdown estimators for model refinement that is computed
    breakdown_estimators_coarsening : np.ndarray of floats
        numpy 2D array that contains the values of all breakdown estimators for model coarsening in each mesh cell
    breakdown_estimators_refinement : np.ndarray of floats
        numpy 2D array that contains the values of all breakdown estimators for model refinement in each mesh cell
    tols_coarsening = list of floats
        list containing the tolerance values for model coarsening
        Note: if type_model_error_estimator == 'heuristics_plus_discretization', then
        first entry of the list is used for dom_decomp_val_res1, the second is used for dom_decomp_val_res2,
        and the third is used for the heuristics
    tols_refinement = list of floats
        list containing the tolerance values for model refinement
        Note: if type_model_error_estimator == 'heuristics_plus_discretization', then
        first entry of the list is used for the last entry of the source term, and the remaining entries
        are used for the remaining heuristic values (the tolerance for all remaining heuristic values is one value) 
    order_diff = int
        the number of moments that is added or subtracted when model refinement or model coarsening, respectively,
        is indicated

    Overriden methods from parent class Simulation
    ----------------------------------------------
    def __init__(self,max_order,min_order,pde_type,mesh,CFL_number,boundary_condition,initial_condition,smoothing,
                 smooth_par,interpolation,path_conservation,spatial_discretization_interior,spatial_discretization_interface,
                 spatial_discretization_predictor,time_integration,two_step_domdecomp_approx,hierarchical,
                 type_model_error_estimator,time_step_splitting,order_diff,tols_coarsening,tols_refinement):
        initializes the ModelAdaptiveMomentSimulation1D object and initializes the attributes


    Instance methods
    ----------------
    decompose_domain(self,flags_decrease,flags_increase)
        given the flags for model coarsening and model refinement, update the orders in each cell

    
    Implemented methods from interface Simulation
    -------
    def run_simulation(self,t_end):
        runs the simulation and outputs the end values
    def _get_initial_conditions(self,cell_centers_x):
        constructs the initial values in each grid cell
    def _update_boundary_conditions(self,values_boundary):
        updates the boundary conditions
    def _post_processing(self,values):
        post processes the end data of the simulation and prepares it for plotting
    """
    def __init__(self, 
                 max_order: int,
                 min_order: int,
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 CFL_number : float,
                 boundary_condition: str,
                 initial_condition: str,
                 smoothing: bool,
                 smooth_par: int,
                 interpolation: bool,
                 path_conservation: bool,
                 spatial_discretization_interior: spatialDiscretization.SpatialDiscretization,
                 spatial_discretization_interface: spatialDiscretization.SpatialDiscretization,
                 spatial_discretization_predictor: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration,
                 two_step_domdecomp_approx: bool,
                 hierarchical: bool,
                 type_model_error_estimator: str,
                 time_step_splitting: bool,
                 order_diff: int,
                 tols_coarsening: list,
                 tols_refinement: list):

        """
        Initializes the ModelAdaptiveMomentSimulations1D object

        Parameters
        ----------
        max_order : int
            maximum order of the moment model during the simulation
        min_order : int
            minimum order of the moment model during the simulation
        pde_type : PDE
            the partial differential equations that are simulated
        mesh : Mesh
            the discrete mesh on which the partial differential equations are numerically solved
        CFL_number : float
            the CFL constant that is used to compute the time step size
        boundary_condition : str
            the used boundary condition
        initial_condition : str
            the initial condition for the simulation
        smoothing : bool
            true if the domain decomposition is smoothed (grouping of cells and giving them the same order),
            false if not
        smooth_par : int
            the maximum number of subdomains (of equal size), if the domain decomposition is smoothed
        interpolation : bool
            true if added moments (after model refinement) are initialized using interpolation,
            otherwise they are set to zero
        path_conservative_interface_coupling : bool
            true if a path-conservative interface coupling method is used, false if otherwise
        spatial_discretization_interior : SpatialDiscretization
            the spatial discretization scheme for the interior of the subdomains
        spatial_discretization_predictor : SpatialDiscretization
            the spatial discretization scheme for the predictions that may be used in the domain decomposition criteria
        spatial_discretization_interface : SpatialDiscretization
            the spatial discretization scheme that is used at the boundary interfaces
        time_integration : TimeIntegration
            the time integration method for the right-hand side source term
        two_step_domdecomp_approx : bool
            true if the model error estimators for the domain decomposition criteria are based on two sets of
            values (at the current time step and predictions at the next time step), false if not
        hierarchical : bool
            true if the moment model is strictly hierarchical
        type_model_error_estimator : str
            the type of model error estimator that is used for computing domain decomposition criteria
        time_step_splitting: bool
            true if a splitting scheme is used for the finite volume updates, false if not
        order_diff = int
            the number of moments that is added or subtracted when model refinement or model coarsening, respectively,
            is indicated
        tols_coarsening = list of floats
            list containing the tolerance values for model coarsening
            Note: if type_model_error_estimator == 'heuristics_plus_discretization', then
            first entry of the list is used for dom_decomp_val_res1, the second is used for dom_decomp_val_res2,
            and the third is used for the heuristics
        tols_refinement = list of floats
            list containing the tolerance values for model refinement
            Note: if type_model_error_estimator == 'heuristics_plus_discretization', then
            first entry of the list is used for the last entry of the source term, and the remaining entries
            are used for the remaining heuristic values (the tolerance for all remaining heuristic values is one value) 
        
        Returns
        -------

        """

        super().__init__(pde_type,
                         mesh,
                         boundary_condition,
                         initial_condition,
                         time_integration,
                         CFL_number)

        self.spatial_discretization_interior = spatial_discretization_interior
        self.spatial_discretization_predictor = spatial_discretization_predictor
        self.spatial_discretization_interface = spatial_discretization_interface
        self.path_conservative_interface_coupling = path_conservation

        self.boundary_interfaces_discretized = [np.floor_divide(self.mesh.resolution,3),np.floor_divide(2*self.mesh.resolution,3)]

        self.max_order = max_order
        self.min_order = min_order
        self.max_number_of_variables = self.pde_type.compute_number_of_variables(self.max_order)
        self.min_number_of_variables = self.pde_type.compute_number_of_variables(self.min_order)

        self.orders = [self.max_order,self.max_order,self.max_order]
        self.numbers_of_variables = [self.max_number_of_variables,self.max_number_of_variables,self.max_number_of_variables]

        self.orders_cellwise = np.full(self.mesh.resolution+2, self.max_order, dtype=int)
        self.numbers_of_variables_cellwise = np.full(self.mesh.resolution+2, self.max_number_of_variables, dtype=int)

        self.dom_decomp_val_res1 = np.zeros(self.mesh.resolution)
        self.dom_decomp_val_res2 = np.zeros(self.mesh.resolution)

        self.subdomainReconstruction = CellwiseSubdomainReconstruction1D(boundary_condition)
        if smoothing:
            self.subdomainReconstruction = SmoothedSubdomainReconstruction1D(boundary_condition,
                                                                             smooth_par,
                                                                             self.max_order,
                                                                             self.max_number_of_variables,
                                                                             self.mesh.resolution,
                                                                             interpolation)

        self.type_model_error_estimator = type_model_error_estimator

        self.domainDecomposition = OneStepModelErrorApproximation(self.pde_type)
        if two_step_domdecomp_approx:
            self.domainDecomposition = TwoStepModelErrorApproximation(self.spatial_discretization_predictor,
                                                                      self.time_integration,
                                                                      self.pde_type,
                                                                      self.mesh.resolution,
                                                                      self.max_order,
                                                                      self.max_number_of_variables)

        # self.interface_coupling = PaddedBufferCell1D(path_conservation,spatial_discretization_interface)

        if time_step_splitting:
            self.AdaptiveFVMStep1D = SplittedNonHierarchicalAdaptiveFVMStep1D(self.spatial_discretization_interior,
                                                                            self.spatial_discretization_interface,
                                                                            self.time_integration,
                                                                            self.max_number_of_variables,
                                                                            self.mesh.resolution,
                                                                            self.pde_type)
            if hierarchical:
                if self.type_model_error_estimator == 'heuristics_plus_discretization' and hierarchical:
                    self.AdaptiveFVMStep1D = SplittedHierarchicalAdaptiveFVMStep1D(self.spatial_discretization_interior,
                                                                            self.spatial_discretization_interface,
                                                                            self.time_integration,
                                                                            self.max_number_of_variables,
                                                                            self.mesh.resolution,
                                                                            self.pde_type)
                elif self.type_model_error_estimator == 'model_difference':
                    self.AdaptiveFVMStep1D = SplittedNonHierarchicalAdaptiveFVMStep1D(self.spatial_discretization_interior,
                                                                            self.spatial_discretization_interface,
                                                                            self.time_integration,
                                                                            self.max_number_of_variables,
                                                                            self.mesh.resolution,
                                                                            self.pde_type)
                else:
                    print('This type of model-error estimator has not been implemented yet!')
        else:
            print("Only time splitting methods implemented so far!")

        self.number_of_breakdown_estimators_coarsening,self.number_of_breakdown_estimators_refinement =\
            self.pde_type.get_number_of_breakdown_estimators(self.max_number_of_variables,self.type_model_error_estimator)

        self.breakdown_estimators_coarsening = np.zeros((self.mesh.resolution,self.number_of_breakdown_estimators_coarsening))
        self.breakdown_estimators_refinement = np.zeros((self.mesh.resolution,self.number_of_breakdown_estimators_refinement))

        self.tols_coarsening = tols_coarsening
        self.tols_refinement = tols_refinement

        if self.type_model_error_estimator == 'heuristics_plus_discretization':
            for i in range(self.max_order):
                self.tols_refinement.append(tols_refinement[-1])

        self.order_diff = order_diff

    def _update_boundary_conditions(self,
                                    values: np.ndarray,
                                    boundary: str) -> np.ndarray:

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

    def run_simulation(self,
                       t_end: float) -> np.ndarray:
        
        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution #TODO: include the possibility of nonuniform grids

        print("Orders at the beginning of the simulation:",self.orders)
        
        values = self._get_initial_conditions(self.mesh.cell_center_positions)
        
        dom_decomp_val_decrease_transport = np.zeros(self.mesh.resolution)
        dom_decomp_val_decrease_source = np.zeros(self.mesh.resolution)

        step_count = 0
        t = 0

        values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values,'left')
        values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values,'right')
        max_speed = self.pde_type.compute_max_wavespeed(max(self.orders),values)

        delta_t = self.CFL*delta_x/max_speed 

        if self.type_model_error_estimator == 'model_difference':
            values = self.AdaptiveFVMStep1D.runFVMStep(values,
                                                       delta_t,
                                                       delta_x,
                                                       self.orders,
                                                       self.numbers_of_variables,
                                                       self.boundary_interfaces_discretized)
        elif self.type_model_error_estimator == 'heuristics_plus_discretization':
            values,dom_decomp_val_decrease_transport,dom_decomp_val_decrease_source = self.AdaptiveFVMStep1D.runFVMStep(values,
                                                                                                                        delta_t,
                                                                                                                        delta_x,
                                                                                                                        self.orders,
                                                                                                                        self.numbers_of_variables,
                                                                                                                        self.boundary_interfaces_discretized)
        else:
            print('this model error estimator has not been implemented yet')

        step_count += 1

        print()
        print('time: '+str(t))
        print('step size: '+str(delta_t))
        print()

        t += delta_t
        while t < t_end:
            print("current orders: ", self.orders)

            # update boundary conditions
            values[0,:self.numbers_of_variables[0]] = self._update_boundary_conditions(values,'left')
            values[self.mesh.resolution+1,:self.numbers_of_variables[-1]] = self._update_boundary_conditions(values,'right')
            # max_speed = self.pde_type.compute_max_wavespeed(self.max_order,values)
            max_speed = self.pde_type.compute_max_wavespeed(max(self.orders),values)

            delta_t = self.CFL*delta_x/max_speed 

            # prev_values = np.copy(values)
            # predicted_values = np.copy(values)

            # padded_vectors_left,padded_vectors_right = self._construct_padded_vectors(prev_values)

            if self.type_model_error_estimator == 'model_difference':
                flags_increase_max,flags_decrease_min,estimators_increase_max,estimators_decrease_max =\
                    self.domainDecomposition.compute_domain_decomposition(values,
                                                                            delta_t,
                                                                            delta_x,
                                                                            self.orders,
                                                                            self.max_order,
                                                                            self.numbers_of_variables,
                                                                            self.mesh.resolution,
                                                                            self.boundary_interfaces_discretized,
                                                                            self.tols_coarsening,
                                                                            self.tols_refinement,
                                                                            self.number_of_breakdown_estimators_coarsening,
                                                                            self.number_of_breakdown_estimators_refinement)

                # self.breakdown_estimators = np.stack((estimators_decrease_max, estimators_increase_max),axis=1)
                self.breakdown_estimators_coarsening = estimators_decrease_max
                self.breakdown_estimators_refinement = estimators_increase_max
                self.breakdown_criteria_flags = self.decompose_domain(flags_decrease_min,
                                                                      flags_increase_max)
                
            elif self.type_model_error_estimator == 'heuristics_plus_discretization':

                flags_increase_max,flags_decrease_min,estimators_increase_max_heur,estimators_decrease_max_heur =\
                    self.domainDecomposition.compute_domain_decomposition(values,
                                                                            delta_t,
                                                                            delta_x,
                                                                            self.orders,
                                                                            self.max_order,
                                                                            self.numbers_of_variables,
                                                                            self.mesh.resolution,
                                                                            self.boundary_interfaces_discretized,
                                                                            [self.tols_coarsening[0]],
                                                                            self.tols_refinement,
                                                                            1,
                                                                            self.number_of_breakdown_estimators_refinement)

                # dom_decomp_val_increase_source_last_entry = 0
                dom_decomp_val_decrease_transport = dom_decomp_val_decrease_transport.reshape(-1, 1)
                dom_decomp_val_decrease_source = dom_decomp_val_decrease_source.reshape(-1, 1)
                self.breakdown_estimators_coarsening = np.stack((estimators_decrease_max_heur,
                                                        dom_decomp_val_decrease_transport,
                                                        dom_decomp_val_decrease_source), axis = 1)
                self.breakdown_estimators_refinement = estimators_increase_max_heur

                for i in range(self.mesh.resolution):
                    if flags_increase_max[i] == 0:
                        # if self.breakdown_estimators_refinement[i,-1] > self.tols_refinement[1]:
                        #     flags_increase_max[i] = 1
                        #     continue
                        # handled = False
                        # for j in range(self.breakdown_estimators_refinement.shape[0]-1):
                        #     if self.breakdown_estimators_refinement[i,j] > self.tols_refinement[0]:
                        #         flags_increase_max[i] = 1 #Generalize this so that it can include general increases
                        #         handled = True
                        #         break                                
                        # if handled:
                        #     continue
                        # if self.breakdown_estimators_coarsening[i,1] < self.tols_coarsening[1] and \
                        #     self.breakdown_estimators_coarsening[i,2] < self.tols_coarsening[2] and\
                        #         self.breakdown_estimators_coarsening[i,0] < self.tols_coarsening[0]:
                        #     flags_decrease_min[i] = -1 #generalize this so that it can include general decreases
                        if flags_decrease_min[i] == 0:
                            if dom_decomp_val_decrease_transport[i] < self.tols_coarsening[0] or dom_decomp_val_decrease_source[i] < self.tols_coarsening[1]:
                                flags_decrease_min[i] = -1
                    
                
                self.breakdown_criteria_flags = self.decompose_domain(flags_decrease_min,
                                                                        flags_increase_max)   

            values,self.orders,self.numbers_of_variables,self.orders_cellwise,self.numbers_of_variables_cellwise,\
                self.boundary_interfaces_discretized = self.subdomainReconstruction.reconstruct_subdomains(values,
                                                                         self.orders_cellwise,
                                                                         self.numbers_of_variables_cellwise,
                                                                         self.max_order,
                                                                         self.mesh.resolution,
                                                                         not self.path_conservative_interface_coupling,
                                                                         self.breakdown_criteria_flags)

            if self.type_model_error_estimator == 'model_difference':
                values = self.AdaptiveFVMStep1D.runFVMStep(values,
                                                           delta_t,
                                                            delta_x,
                                                            self.orders,
                                                            self.numbers_of_variables,
                                                            self.boundary_interfaces_discretized)
            elif self.type_model_error_estimator == 'heuristics_plus_discretization':
                values,dom_decomp_val_decrease_transport,dom_decomp_val_decrease_source = self.AdaptiveFVMStep1D.runFVMStep(values,
                                                                                                                            delta_t,
                                                                                                                            delta_x,
                                                                                                                            self.orders,
                                                                                                                            self.numbers_of_variables,
                                                                                                                            self.boundary_interfaces_discretized)
            else:
                print('this model error estimator has not been implemented yet')


            step_count += 1

            print()
            print('time: '+str(t))
            print('step size: '+str(delta_t))
            print()

            t+=delta_t

        print(self.boundary_interfaces_discretized)
        values = simulation_data = self._post_processing(values)
        return simulation_data     

    def decompose_domain(self,
                        flags_decrease: np.ndarray,
                        flags_increase: np.ndarray) -> np.ndarray:       
        
        """
        given the model-error estimators for model-coarsening (flags_decrease) and the model-error estimators 
        for model-refinement (flags_increase), compute the final domain decomposition by updating the
        cellwise orders

        Parameters
        ----------
        flags_decrease: np.ndarray
            the computed flags for model-coarsening (decreasing the order)
        flags_inrease: np.ndarray
            the computed flags for model-refinement (increasing the order)
            
        Returns
        -------
        domain_decomposition_flags : np.ndarray
            numpy array where entry with index i contains the computed order difference
            (positive if increase, negative if decrease)

        """

        domain_decomposition_flags = np.zeros(self.mesh.resolution,dtype=int)
        for i in range(self.mesh.resolution):
            if flags_increase[i] > 0:
                if self.orders_cellwise[i+1] <= self.max_order - self.order_diff:
                    domain_decomposition_flags[i] = self.order_diff
            else:
                if flags_decrease[i] < 0:
                    if self.orders_cellwise[i+1] >= self.min_order + self.order_diff:
                        domain_decomposition_flags[i] = -self.order_diff

        return domain_decomposition_flags


class AdaptiveFVMStep1D(ABC):

    """
    This abstract class represents one time step of a spatially adaptive finite volume discretization in
    one-dimensional physical space.

    ...

    Attributes
    ----------
    spat_disc_interior : SpatialDiscretization
        the spatial discretization that is used in the interior of the subdomains
    spat_disc_interface : SpatialDiscretization
        the spatial discretization that is used at the boundary interfaces between the different-order subdomains
    max_number_of_variables : int
        the number of variables in the maximum-order moment model
    mesh_resolution : int
        the resolution of the spatial mesh
    pde_type : PDE
        the PDE (moment model) that is simulated
    time_integration : TimeIntegration
        the time integration method

    Instance methods
    ----------------    
    def __init__(self,spatial_discretization_interior,spatial_discretization_interface,time_integration,
             max_number_of_variables,mesh_resolution,pde_type):
        initializes the AdaptiveFVMStep1D object
    
    Abstract methods
    ----------------
    def runFVMStep(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces):
        runs one time step of the finite volume scheme
    """ 

    def __init__(self,
             spatial_discretization_interior: spatialDiscretization.SpatialDiscretization,
             spatial_discretization_interface: spatialDiscretization.SpatialDiscretization,
             time_integration: timeIntegration.TimeIntegration,
             max_number_of_variables: int,
             mesh_resolution: int,
             pde_type: pde.PDE):

        """
        Initializes the AdaptiveFVMStep1D and constructs all the attributes.

        Parameters
        ----------
        spat_disc_interior : SpatialDiscretization
            the spatial discretization that is used in the interior of the subdomains
        spat_disc_interface : SpatialDiscretization
            the spatial discretization that is used at the boundary interfaces between the different-order subdomains
        max_number_of_variables : int
            the number of variables in the maximum-order moment model
        mesh_resolution : int
            the resolution of the spatial mesh
        pde_type : PDE
            the PDE (moment model) that is simulated
        time_integration : TimeIntegration
            the time integration method

        """

        self.spat_disc_interior = spatial_discretization_interior
        self.spat_disc_interface = spatial_discretization_interface
        self.max_number_of_variables = max_number_of_variables
        self.mesh_resolution = mesh_resolution
        self.pde_type = pde_type
        self.time_integration = time_integration

    @abstractmethod
    def runFVMStep(self,
                   values: np.ndarray,
                   delta_t: float,
                   delta_x: float,
                   orders: list,
                   numbers_of_variables: int,
                   boundary_interfaces: list) -> np.ndarray:
        """
        Runs one time step of the adaptive finite volume scheme.

        Parameters
        ----------
        values : np.ndarray
            the initial values of the state variables
        delta_t : float
            the time step size
        delta_x : float
            the mesh cell size
        orders : list of integers
            the orders in each subdomain
        numbers_of_variables : int
            the number of variables in the current moment model
        boundary_interfaces : list of integers
            the indices that denote the boundary interfaces between different-order subdomains


        Returns
        -------
        new_values : np.ndarray
            the values at the end of the time step

        """
        pass

class SplittedAdaptiveFVMStep1D(AdaptiveFVMStep1D):

    """
    This abstract class represents one time step of a splitted spatially adaptive finite volume discretization 
    in one-dimensional physical space.

    ...

    Attributes
    ----------
    spat_disc_interior : SpatialDiscretization
        the spatial discretization that is used in the interior of the subdomains
    spat_disc_interface : SpatialDiscretization
        the spatial discretization that is used at the boundary interfaces between the different-order subdomains
    max_number_of_variables : int
        the number of variables in the maximum-order moment model
    mesh_resolution : int
        the resolution of the spatial mesh
    pde_type : PDE
        the PDE (moment model) that is simulated
    time_integration : TimeIntegration
        the time integration method

    Methods inherited from parent abstract class AdaptiveFVMStep1D
    -------------------------------------------------------------    
    def __init__(self,spatial_discretization_interior,spatial_discretization_interface,time_integration,
             max_number_of_variables,mesh_resolution,pde_type):
        initializes the AdaptiveFVMStep1D object
    
    Implemented methods from parent abstract class AdaptiveFVMStep1D
    ---------------------------------------------------------------------------
    def runFVMStep(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces):
        runs one time step of the finite volume scheme

    Abstract methods
    ----------------
    simulate_transport_step(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        simulates the transport step in the splitting scheme
    simulate_source_step(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        simulates the source step in the splitting scheme    
    """ 

    def runFVMStep(self,
                   values: np.ndarray,
                   delta_t: float,
                   delta_x: float,
                   orders: list,
                   numbers_of_variables: int,
                   boundary_interfaces: list) -> np.ndarray:
        new_values = self.simulate_transport_step(values,
                                        delta_t,
                                        delta_x,
                                        orders,
                                        numbers_of_variables,
                                        boundary_interfaces)
        new_values = self.simulate_source_step(new_values,
                                       delta_t,
                                       delta_x,
                                       orders,
                                       numbers_of_variables,
                                       boundary_interfaces)
        
        
        return new_values

    @abstractmethod
    def simulate_transport_step(self,
                        values: np.ndarray,
                        delta_t: float,
                        delta_x: float,
                        orders: list,
                        numbers_of_variables: int,
                        boundary_interfaces: list) -> np.ndarray:
        
        """
        Runs one transport step in the splitting scheme.

        Parameters
        ----------
        values : np.ndarray
            the initial values of the state variables
        delta_t : float
            the time step size
        delta_x : float
            the mesh cell size
        orders : list of integers
            the orders in each subdomain
        numbers_of_variables : int
            the number of variables in the current moment model
        boundary_interfaces : list of integers
            the indices that denote the boundary interfaces between different-order subdomains


        Returns
        -------
        values : np.ndarray
            the values at the end of the transport time step

        """
        
        pass

    @abstractmethod
    def simulate_source_step(self,
                        values: np.ndarray,
                        delta_t: float,
                        delta_x: float,
                        orders: list,
                        numbers_of_variables: int,
                        boundary_interfaces: list) -> np.ndarray:
        """
        Runs one source step in the splitting scheme.

        Parameters
        ----------
        values : np.ndarray
            the initial values of the state variables
        delta_t : float
            the time step size
        delta_x : float
            the mesh cell size
        orders : list of integers
            the orders in each subdomain
        numbers_of_variables : int
            the number of variables in the current moment model
        boundary_interfaces : list of integers
            the indices that denote the boundary interfaces between different-order subdomains


        Returns
        -------
        values : np.ndarray
            the values at the end of the source time step

        """
        
        pass

class SplittedHierarchicalAdaptiveFVMStep1D(SplittedAdaptiveFVMStep1D):

    """
    This class represents one time step of a splitted spatially adaptive finite volume discretization 
    for strictly hierarchical moment models in one-dimensional physical space.

    ...

    Attributes
    ----------
    spat_disc_interior : SpatialDiscretization
        the spatial discretization that is used in the interior of the subdomains
    spat_disc_interface : SpatialDiscretization
        the spatial discretization that is used at the boundary interfaces between the different-order subdomains
    max_number_of_variables : int
        the number of variables in the maximum-order moment model
    mesh_resolution : int
        the resolution of the spatial mesh
    pde_type : PDE
        the PDE (moment model) that is simulated
    time_integration : TimeIntegration
        the time integration method

    Methods inherited from abstract class AdaptiveFVMStep1D
    -------------------------------------------------------------    
    def __init__(self,spatial_discretization_interior,spatial_discretization_interface,time_integration,
             max_number_of_variables,mesh_resolution,pde_type):
        initializes the AdaptiveFVMStep1D object
    
    Overriden methods from parent abstract class SplittedAdaptiveFVMStep1D
    -------------------------------------------------------------------------
    def runFVMStep(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces):
        runs one time step of the finite volume scheme

    Implemented methods from parent abstract class SplittedAdaptiveFVMStep1D
    ------------------------------------------------------------------------
    simulate_transport_step(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        simulates the transport step in the splitting scheme
    simulate_source_step(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        simulates the source step in the splitting scheme    
    """ 

    def runFVMStep(self,
                   values: np.ndarray,
                   delta_t: float,
                   delta_x: float,
                   orders: list,
                   numbers_of_variables: int,
                   boundary_interfaces: list) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
        """
        Runs one time step of the adaptive finite volume scheme.

        Parameters
        ----------
        values : np.ndarray
            the initial values of the state variables
        delta_t : float
            the time step size
        delta_x : float
            the mesh cell size
        orders : list of integers
            the orders in each subdomain
        numbers_of_variables : int
            the number of variables in the current moment model
        boundary_interfaces : list of integers
            the indices that denote the boundary interfaces between different-order subdomains


        Returns
        -------
        new_values : np.ndarray
            the values at the end of the time step
        dom_decomp_val_transport : np.ndarray
            the values for the transport residual in the highest-order moments
        dom_decomp_val_source : np.ndarray
            the values for the source term residual in the highest-order moments

        """

        dom_decomp_val_transport,new_values = self.simulate_transport_step(values,
                                        delta_t,
                                        delta_x,
                                        orders,
                                        numbers_of_variables,
                                        boundary_interfaces)
        dom_decomp_val_source,new_values = self.simulate_source_step(new_values,
                                       delta_t,
                                       delta_x,
                                       orders,
                                       numbers_of_variables,
                                       boundary_interfaces)
        
        return new_values,dom_decomp_val_transport,dom_decomp_val_source

    def simulate_transport_step(self,
                   values: np.ndarray,
                   delta_t: float,
                   delta_x: float,
                   orders: list,
                   numbers_of_variables: int,
                   boundary_interfaces: list) -> tuple[np.ndarray,np.ndarray]:
        
        """
        Runs one transport step in the splitting scheme.

        Parameters
        ----------
        values : np.ndarray
            the initial values of the state variables
        delta_t : float
            the time step size
        delta_x : float
            the mesh cell size
        orders : list of integers
            the orders in each subdomain
        numbers_of_variables : int
            the number of variables in the current moment model
        boundary_interfaces : list of integers
            the indices that denote the boundary interfaces between different-order subdomains


        Returns
        -------
        dom_decomp_val_transport : np.ndarray
            the values for the transport residual in the highest-order moments
        values : np.ndarray
            the values at the end of the transport time step

        """

        fluctuations_min = np.zeros((self.mesh_resolution+1,self.max_number_of_variables))
        fluctuations_plus = np.zeros((self.mesh_resolution+1,self.max_number_of_variables))
        res1_min = np.zeros((self.mesh_resolution+1,self.max_number_of_variables))
        res1_plus = np.zeros((self.mesh_resolution+1,self.max_number_of_variables))
        dom_decomp_val_transport = np.zeros(self.mesh_resolution)

        prev_values = np.copy(values)

        right_boundary_subdomain = 0

        for m in range(len(boundary_interfaces)):
            order = orders[m]
            n_variables = numbers_of_variables[m]
            def system_matrix(cell_values):
                return self.pde_type.compute_system_matrix(order,
                                                           cell_values)               

            left_boundary_subdomain = right_boundary_subdomain + 1
            right_boundary_subdomain = boundary_interfaces[m]                 

            # if padding:
            #     if order < orders[m+1]:
            #         prev_values[right_boundary_subdomain,:] = padded_vectors_left[m,:]
            #     else:
            #         prev_values[right_boundary_subdomain+1,:] = padded_vectors_right[m,:]    

            if m > 0:
                j = left_boundary_subdomain
                generalized_roe_minus, generalized_roe_plus = self.spat_disc_interface.compute_generalized_roe_and_viscosity(
                    prev_values[j-1,:n_variables],
                    prev_values[j,:n_variables],
                    system_matrix,
                    delta_t,
                    delta_x)
                fluctuations_plus[j-1,:n_variables] =\
                    generalized_roe_plus@(prev_values[j,:n_variables]-prev_values[j-1,:n_variables])
                fluctuations_min[j-1,:n_variables] =\
                    generalized_roe_minus@(prev_values[j,:n_variables]-prev_values[j-1,:n_variables]) 

                left_boundary_subdomain += 1

            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                generalized_roe_minus, generalized_roe_plus = self.spat_disc_interior.compute_generalized_roe_and_viscosity(
                    prev_values[i-1,:n_variables],
                    prev_values[i,:n_variables],
                    system_matrix,
                    delta_t,
                    delta_x)
                fluctuations_plus[i-1,:n_variables] =\
                    generalized_roe_plus@(prev_values[i,:n_variables]-prev_values[i-1,:n_variables])
                fluctuations_min[i-1,:n_variables] =\
                    generalized_roe_minus@(prev_values[i,:n_variables]-prev_values[i-1,:n_variables]) 
                res1_min[i-1,:n_variables-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables-1]-prev_values[i-1,n_variables-1])
                res1_plus[i-1,:n_variables-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables-1]-prev_values[i-1,n_variables-1])
  
        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        def system_matrix(cell_values):
            return self.pde_type.compute_system_matrix(order,cell_values)
        
        for i in range(right_boundary_subdomain+1,self.mesh_resolution+2):
            generalized_roe_minus, generalized_roe_plus = self.spat_disc_interior.compute_generalized_roe_and_viscosity(
                prev_values[i-1,:n_variables],
                prev_values[i,:n_variables],
                system_matrix,
                delta_t,
                delta_x)
            fluctuations_plus[i-1,:n_variables] =\
                generalized_roe_plus@(prev_values[i,:n_variables]-prev_values[i-1,:n_variables])
            fluctuations_min[i-1,:n_variables] =\
                generalized_roe_minus@(prev_values[i,:n_variables]-prev_values[i-1,:n_variables]) 
            res1_min[i-1,:n_variables-1] = generalized_roe_minus[:-1,-1]*(prev_values[i,n_variables-1]-prev_values[i-1,n_variables-1])
            res1_plus[i-1,:n_variables-1] = generalized_roe_plus[:-1,-1]*(prev_values[i,n_variables-1]-prev_values[i-1,n_variables-1])

        right_boundary_subdomain = 0

        for m in range(len(boundary_interfaces)):
            n_variables = numbers_of_variables[m]

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = boundary_interfaces[m]
            
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                values[i,:n_variables] = prev_values[i,:n_variables]\
                    - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables]+fluctuations_min[i,:n_variables]) 
                
                dom_decomp_val_transport[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables-1]+res1_min[i,:n_variables-1])

        n_variables = numbers_of_variables[-1]

        for i in range(right_boundary_subdomain+1,self.mesh_resolution+1):
            values[i,:n_variables] = prev_values[i,:n_variables]\
                - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables]+fluctuations_min[i,:n_variables]) 
            
            dom_decomp_val_transport[i-1] = np.linalg.norm(res1_plus[i-1,:n_variables-1]+res1_min[i,:n_variables-1])

        dom_decomp_val_transport = dom_decomp_val_transport/delta_x

        return dom_decomp_val_transport,values

    def simulate_source_step(self,
                   values: np.ndarray,
                   delta_t: float,
                   delta_x: float,
                   orders: list,
                   numbers_of_variables: int,
                   boundary_interfaces: list) -> tuple[np.ndarray,np.ndarray]:

        """
        Runs one source step in the splitting scheme.

        Parameters
        ----------
        values : np.ndarray
            the initial values of the state variables
        delta_t : float
            the time step size
        delta_x : float
            the mesh cell size
        orders : list of integers
            the orders in each subdomain
        numbers_of_variables : int
            the number of variables in the current moment model
        boundary_interfaces : list of integers
            the indices that denote the boundary interfaces between different-order subdomains


        Returns
        -------
        dom_decomp_val_source : np.ndarray
            the values for the source term residual in the highest-order moments
        values : np.ndarray
            the values at the end of the source time step

        """

        dom_decomp_val_source = np.zeros(self.mesh_resolution)

        prev_values = np.copy(values)

        right_boundary_subdomain = 0

        for m in range(len(boundary_interfaces)):
            order = orders[m]
            n_variables = numbers_of_variables[m]

            def source_term(cell_values,delta_t):
                return self.pde_type.compute_source_term(order,cell_values,delta_t)

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = boundary_interfaces[m]
            
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                values[i,:n_variables] = self.time_integration.integrate(values[i,:n_variables],source_term,delta_t)
                dom_decomp_val_source[i-1] = np.abs(values[i,n_variables-1]-prev_values[i,n_variables-1])/delta_t
            
            # if padding:
            #     if order < orders[m+1]:
            #         values[right_boundary_subdomain,n_variables:] = \
            #             values[right_boundary_subdomain+1,n_variables:] # update boundary interface boundary condition 
            #     else:
            #         values[right_boundary_subdomain+1,numbers_of_variables[m+1]:] = \
            #             values[right_boundary_subdomain,numbers_of_variables[m+1]:] # update boundary interface boundary condition 

        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        def source_term(cell_values,delta_t):
            return self.pde_type.compute_source_term(order,cell_values,delta_t)         
                    
        for i in range(right_boundary_subdomain+1,self.mesh_resolution+1): 
            values[i,:n_variables] = self.time_integration.integrate(values[i,:n_variables],source_term,delta_t)
            dom_decomp_val_source[i-1] = np.abs(values[i,n_variables-1]-prev_values[i,n_variables-1])/delta_t

        # if padding:
        #     if order < orders[m-1]:
        #         values[right_boundary_subdomain+1,n_variables:] = \
        #             values[right_boundary_subdomain,n_variables:] # update boundary interface boundary condition 
        #     else:
        #         values[right_boundary_subdomain,numbers_of_variables[m-1]:] = \
        #             values[right_boundary_subdomain+1,numbers_of_variables[m-1]:] # update boundary interface boundary condition 

        return dom_decomp_val_source,values    

class SplittedNonHierarchicalAdaptiveFVMStep1D(SplittedAdaptiveFVMStep1D):
    """
    This class represents one time step of a splitted spatially adaptive finite volume discretization 
    for moment models that are not stricly hierarchical, in one-dimensional physical space.

    ...

    Attributes
    ----------
    spat_disc_interior : SpatialDiscretization
        the spatial discretization that is used in the interior of the subdomains
    spat_disc_interface : SpatialDiscretization
        the spatial discretization that is used at the boundary interfaces between the different-order subdomains
    max_number_of_variables : int
        the number of variables in the maximum-order moment model
    mesh_resolution : int
        the resolution of the spatial mesh
    pde_type : PDE
        the PDE (moment model) that is simulated
    time_integration : TimeIntegration
        the time integration method

    Methods inherited from abstract class AdaptiveFVMStep1D
    -------------------------------------------------------------    
    def __init__(self,spatial_discretization_interior,spatial_discretization_interface,time_integration,
             max_number_of_variables,mesh_resolution,pde_type):
        initializes the AdaptiveFVMStep1D object
    
    Methods inherited from abstract parent class SplittedAdaptiveFVMStep1D
    -------------------------------------------------------------------------
    def runFVMStep(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces):
        runs one time step of the finite volume scheme

    Implemented methods from abstract parent class SplittedAdaptiveFVMStep1D
    ------------------------------------------------------------------------
    simulate_transport_step(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        simulates the transport step in the splitting scheme
    simulate_source_step(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        simulates the source step in the splitting scheme    
    """ 
    def simulate_transport_step(self,
                   values: np.ndarray,
                   delta_t: float,
                   delta_x: float,
                   orders: list,
                   numbers_of_variables: int,
                   boundary_interfaces: list) -> tuple[np.ndarray,np.ndarray]:
        
        fluctuations_min = np.zeros((self.mesh_resolution+1,self.max_number_of_variables))
        fluctuations_plus = np.zeros((self.mesh_resolution+1,self.max_number_of_variables))

        prev_values = np.copy(values)

        right_boundary_subdomain = 0

        for m in range(len(boundary_interfaces)):
            order = orders[m]
            n_variables = numbers_of_variables[m]
            def system_matrix(cell_values):
                return self.pde_type.compute_system_matrix(order,cell_values)               

            left_boundary_subdomain = right_boundary_subdomain + 1
            right_boundary_subdomain = boundary_interfaces[m]                 

            # if padding:
            #     if order < orders[m+1]:
            #         prev_values[right_boundary_subdomain,:] = padded_vectors_left[m,:]
            #     else:
            #         prev_values[right_boundary_subdomain+1,:] = padded_vectors_right[m,:]    

            if m > 0:
                j = left_boundary_subdomain
                fluctuations_min[j-1,:n_variables], fluctuations_plus[j-1,:n_variables],  =\
                    self.spat_disc_interface.compute_fluctuation(
                                                        prev_values[j-1,:n_variables],
                                                        prev_values[j,:n_variables],
                                                        system_matrix,
                                                        delta_t,
                                                        delta_x) 
                left_boundary_subdomain += 1
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                fluctuations_min[i-1,:n_variables], fluctuations_plus[i-1,:n_variables] =\
                    self.spat_disc_interior.compute_fluctuation(
                                                        prev_values[i-1,:n_variables],
                                                        prev_values[i,:n_variables],
                                                        system_matrix,
                                                        delta_t,
                                                        delta_x)   
        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        def system_matrix(cell_values):
            return self.pde_type.compute_system_matrix(order,cell_values)
        
        for i in range(right_boundary_subdomain+1,self.mesh_resolution+2):
                fluctuations_min[i-1,:n_variables], fluctuations_plus[i-1,:n_variables] =\
                    self.spat_disc_interior.compute_fluctuation(
                                                        prev_values[i-1,:n_variables],
                                                        prev_values[i,:n_variables],
                                                        system_matrix,
                                                        delta_t,
                                                        delta_x) 

        right_boundary_subdomain = 0

        for m in range(len(boundary_interfaces)):
            n_variables = numbers_of_variables[m]

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = boundary_interfaces[m]
            
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                values[i,:n_variables] = prev_values[i,:n_variables]\
                    - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables]+fluctuations_min[i,:n_variables]) 

        n_variables = numbers_of_variables[-1]

        for i in range(right_boundary_subdomain+1,self.mesh_resolution+1):
            values[i,:n_variables] = prev_values[i,:n_variables]\
                - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables]+fluctuations_min[i,:n_variables]) 

        return values      

    def simulate_source_step(self,
                   values: np.ndarray,
                   delta_t: float,
                   delta_x: float,
                   orders: list,
                   numbers_of_variables: int,
                   boundary_interfaces: list) -> tuple[np.ndarray,np.ndarray]:

        right_boundary_subdomain = 0

        for m in range(len(boundary_interfaces)):
            order = orders[m]
            n_variables = numbers_of_variables[m]

            def source_term(cell_values,delta_t):
                return self.pde_type.compute_source_term(order,cell_values,delta_t)

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = boundary_interfaces[m]
            
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                values[i,:n_variables] = self.time_integration.integrate(values[i,:n_variables],source_term,delta_t)
            
            # if padding:
            #     if order < orders[m+1]:
            #         values[right_boundary_subdomain,n_variables:] = \
            #             values[right_boundary_subdomain+1,n_variables:] # update boundary interface boundary condition 
            #     else:
            #         values[right_boundary_subdomain+1,numbers_of_variables[m+1]:] = \
            #             values[right_boundary_subdomain,numbers_of_variables[m+1]:] # update boundary interface boundary condition 

        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        def source_term(cell_values,delta_t):
            return self.pde_type.compute_source_term(order,cell_values,delta_t)         
                    
        for i in range(right_boundary_subdomain+1,self.mesh_resolution+1): 
            values[i,:n_variables] = self.time_integration.integrate(values[i,:n_variables],source_term,delta_t)

        # if padding:
        #     if order < orders[m-1]:
        #         values[right_boundary_subdomain+1,n_variables:] = \
        #             values[right_boundary_subdomain,n_variables:] # update boundary interface boundary condition 
        #     else:
        #         values[right_boundary_subdomain,numbers_of_variables[m-1]:] = \
        #             values[right_boundary_subdomain+1,numbers_of_variables[m-1]:] # update boundary interface boundary condition 

        return values 



class DomainDecomposition1D(ABC):

    """
    This abstract class represents a domain decomposition strategy that updates the orders of the 
    moment models in each mesh cell, for 1D moment models.

    ...

    Attributes
    ----------
    pde_type : PDE
        the PDE (moment model) that is simulated

    Instance methods
    ----------------  
    def __init__(self,pde_type):
        initializes the DomainDecomposition1D object
    def compute_coarsening_criterion(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,n,
            boundary_interfaces,tols_coarsening,increase_criterion_flags,number_of_coarsening_estimators_pde)
        computes the values of the coarsening estimator, and the flags for model coarsening
    def compute_refinement_criterion(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,n,
            boundary_interfaces,tols_coarsening,number_of_coarsening_estimators_pde)
        computes the values of the refinement estimator, and the flags for model refinement

    abstract methods
    ----------------
    def compute_domain_decomposition(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,
                            mesh_resolution,boundary_interfaces_discretized,tols_coarsening,tols_refinement,
                            number_of_coarsening_estimators_pde,number_of_refinement_estimators_pde)
        computes the updated model orders in each mesh cell
  
    """ 

    def __init__(self,
                pde_type: pde.PDE):
        """
        Initializes the DomainDecomposition1D object and constructs all the attributes.

        Parameters
        ----------
        pde_type : PDE
            the PDE (moment model) that is simulated

        """

        self.pde_type = pde_type

    def compute_coarsening_criterion(self,
                                    values: np.ndarray,
                                    delta_t: float,
                                    delta_x: float,
                                    orders: list,
                                    max_order: int,
                                    numbers_of_variables: list,
                                    n: int,
                                    boundary_interfaces: list,
                                    tols_coarsening: list,
                                    increase_criterion_flags: np.ndarray,
                                    number_of_coarsening_estimators_pde: int) -> tuple[np.ndarray,np.ndarray]:

        """
        Computes the model error estimator for model coarsening, compares the error estimators with
        threshold values, and flags cells for model coarsening.

        Parameters
        ----------
        values: np.ndarray
            the current values of the state variables in each mesh cell
        delta_t: float
            the current time step
        delta_x: float
            the mesh cell size
        orders: list
            the orders in the subdomains
        max_order: int
            the maximum order in the adaptive simulation
        numbers_of_variables: list
            the numbers of variables in the subdomains
        n: int
            the number of mesh cells
        boundary_interfaces: list
            the boundary interfaces between different-order subdomains
        tols_coarsening: list
            the threshold values for model coarsening
        increase_criterion_flags: np.ndarray
            the flags for model refinement
        number_of_coarsening_estimators_pde: int
            the number of estimators for model coarsening that are computed

        Returns
        -------
        model_error_estimators_coarsening: np.ndarray
            numpy 2D array containing the values for all the estimators for model coarsening, in each mesh cell
        coarsening_criterion_flags: np.ndarray
            numpy 1D array, containing the flags for model coarsening (-1 if flagged for model coarsening)

        """

        coarsening_criterion_flags = np.zeros(n,dtype=int)
        model_error_estimators_coarsening = np.zeros((n,number_of_coarsening_estimators_pde))

        input_values = np.copy(values)

        r = 0

        order = orders[0]
        n_variables = numbers_of_variables[0]
        for m in range(len(boundary_interfaces)):
            n_variables_prev = n_variables
            order = orders[m]
            n_variables = numbers_of_variables[m]
            n_variables_next = numbers_of_variables[m+1]
            l = r+1
            r = boundary_interfaces[m]  

            n_left = min(n_variables_prev,n_variables)
            n_right = min(n_variables,n_variables_next)

            left_boundary_value = input_values[l,:]
            right_boundary_value = input_values[r,:]
            left_boundary_value[:n_left] = input_values[l-1,:n_left]
            right_boundary_value[:n_right] = input_values[r+1,:n_right]

            if order > 3:
                if increase_criterion_flags[l-1] == 0:
                    model_error_estimators_coarsening[l-1,:] = self.pde_type.compute_coarsening_estimator(left_boundary_value,
                                                                                                values[l,:],
                                                                                                values[l+1,:],
                                                                                                order,
                                                                                                max_order,
                                                                                                delta_t,
                                                                                                delta_x)                                                       
                for i in range(l+1,r):
                    if increase_criterion_flags[i-1] == 0:
                        model_error_estimators_coarsening[i-1,:] = self.pde_type.compute_coarsening_estimator(values[i-1,:],
                                                                                                    values[i,:],
                                                                                                    values[i+1,:],
                                                                                                    order,
                                                                                                    max_order,
                                                                                                    delta_t,
                                                                                                    delta_x)          
                if increase_criterion_flags[r-1] == 0:
                    model_error_estimators_coarsening[r-1,:] = self.pde_type.compute_coarsening_estimator(values[r-1,:],
                                                                                                values[r,:],
                                                                                                right_boundary_value,
                                                                                                order,
                                                                                                max_order,
                                                                                                delta_t,
                                                                                                delta_x) 
        n_variables_prev = n_variables
        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        n_left = min(n_variables_prev,n_variables)

        l = r+1  

        left_boundary_value = input_values[l,:]            
        left_boundary_value[:n_left] = input_values[l-1,:n_left]  

        if order > 3:
            if increase_criterion_flags[l-1] == 0:
                model_error_estimators_coarsening[l-1,:] = self.pde_type.compute_coarsening_estimator(
                                                                                            left_boundary_value,
                                                                                            values[l,:],
                                                                                            values[l+1,:],
                                                                                            order,
                                                                                            max_order,
                                                                                            delta_t,
                                                                                            delta_x)                                                      
            for i in range(l+1,n+1):
                if increase_criterion_flags[i-1] == 0:
                    model_error_estimators_coarsening[i-1,:] = self.pde_type.compute_coarsening_estimator(values[i-1,:],
                                                                                                values[i,:],
                                                                                                values[i+1,:],
                                                                                                order,
                                                                                                max_order,
                                                                                                delta_t,
                                                                                                delta_x) 

        for i in range(n):
            if increase_criterion_flags[i] == 0:
                coarsening = True
                for j in range(number_of_coarsening_estimators_pde):
                    if model_error_estimators_coarsening[i,j] > tols_coarsening[0]:
                        coarsening = False 
                        break
                if coarsening:
                    coarsening_criterion_flags[i] = -1

        return model_error_estimators_coarsening, coarsening_criterion_flags

    def compute_refinement_criterion(self,
                                    values: np.ndarray,
                                    delta_t: float,
                                    delta_x: float,
                                    orders: list,
                                    max_order: int,
                                    numbers_of_variables: list,
                                    n: int,
                                    boundary_interfaces: list,
                                    tols_refinement: list,
                                    number_of_refinement_estimators_pde: int) -> tuple[np.ndarray,np.ndarray]:

        """
        Computes the model error estimator for model refinement, compares the error estimators with
        threshold values, and flags cells for model refinement.

        Parameters
        ----------
        values: np.ndarray
            the current values of the state variables in each mesh cell
        delta_t: float
            the current time step
        delta_x: float
            the mesh cell size
        orders: list
            the orders in the subdomains
        max_order: int
            the maximum order in the adaptive simulation
        numbers_of_variables: list
            the numbers of variables in the subdomains
        n: int
            the number of mesh cells
        boundary_interfaces: list
            the boundary interfaces between different-order subdomains
        tols_refinement: list
            the threshold values for model refinement
        number_of_refinement_estimators_pde: int
            the number of estimators for model refinement that are computed

        Returns
        -------
        model_error_estimators_refinement: np.ndarray
            numpy 2D array containing the values for all the estimators for model refinement, in each mesh cell
        refinement_criterion_flags: np.ndarray
            numpy 1D array, containing the flags for model refinement (1 if flagged for model refinement)

        """

        refinement_criterion_flags = np.zeros(n,dtype=int)
        model_error_estimators_refinement = np.zeros((n,number_of_refinement_estimators_pde))

        input_values = np.copy(values) 

        r = 0

        order = orders[0]
        n_variables = numbers_of_variables[0]
        for m in range(len(boundary_interfaces)):
            n_variables_prev = n_variables
            order = orders[m]
            n_variables = numbers_of_variables[m]
            n_variables_next = numbers_of_variables[m+1]
            l = r+1
            r = boundary_interfaces[m]  

            n_left = min(n_variables_prev,n_variables)
            n_right = min(n_variables,n_variables_next)

            left_boundary_value = input_values[l,:]
            right_boundary_value = input_values[r,:]
            left_boundary_value[:n_left] = input_values[l-1,:n_left]
            right_boundary_value[:n_right] = input_values[r+1,:n_right]

            model_error_estimators_refinement[l-1,:] = self.pde_type.compute_refinement_estimator(
                                                                                    left_boundary_value,
                                                                                    values[l,:],
                                                                                    values[l+1,:],
                                                                                    order,
                                                                                    max_order,
                                                                                    delta_t,
                                                                                    delta_x) 

            for i in range(l+1,r):
                model_error_estimators_refinement[i-1,:] = self.pde_type.compute_refinement_estimator(values[i-1,:],
                                                                                            values[i,:],
                                                                                            values[i+1,:],
                                                                                            order,
                                                                                            max_order,
                                                                                            delta_t,
                                                                                            delta_x)     
            model_error_estimators_refinement[r-1,:] = self.pde_type.compute_refinement_estimator(values[r-1,:],
                                                                                        values[r,:],
                                                                                        right_boundary_value,
                                                                                        order,
                                                                                        max_order,
                                                                                        delta_t,
                                                                                        delta_x) 

        n_variables_prev = n_variables
        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        n_left = min(n_variables_prev,n_variables)

        l = r+1  

        left_boundary_value = input_values[l,:]            
        left_boundary_value[:n_left] = input_values[l-1,:n_left]  

        model_error_estimators_refinement[l-1,:] = self.pde_type.compute_refinement_estimator(left_boundary_value,
                                                                                values[l,:],
                                                                                values[l+1,:],
                                                                                order,
                                                                                max_order,
                                                                                delta_t,
                                                                                delta_x)                 
                                                
        for i in range(l+1,n+1):
            model_error_estimators_refinement[i-1,:] = self.pde_type.compute_refinement_estimator(values[i-1,:],
                                                                        values[i,:],
                                                                        values[i+1,:],
                                                                        order,
                                                                        max_order,
                                                                        delta_t,
                                                                        delta_x) 
            
        for i in range(n):
            for j in range(number_of_refinement_estimators_pde):
                if model_error_estimators_refinement[i,j] > tols_refinement[j]:
                    refinement_criterion_flags[i] = 1
                    break
           
        return model_error_estimators_refinement, refinement_criterion_flags

    @abstractmethod
    def compute_domain_decomposition(self,
                                     values: np.ndarray,
                                     delta_t: float,
                                     delta_x: float,
                                     orders: list,
                                     max_order: int,
                                     numbers_of_variables: list,
                                     mesh_resolution: int,
                                     boundary_interfaces_discretized: list,
                                     tols_coarsening: list,
                                     tols_refinement: list,
                                     number_of_coarsening_estimators_pde: int,
                                     number_of_refinement_estimators_pde: int) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:

        """
        Computes the model error estimators for model coarsening and model refinemtn, 
        compares the error estimators with threshold values, and flags cells for model coarsening and model refinement.

        Parameters
        ----------
        values: np.ndarray
            the current values of the state variables in each mesh cell
        delta_t: float
            the current time step
        delta_x: float
            the mesh cell size
        orders: list
            the orders in the subdomains
        max_order: int
            the maximum order in the adaptive simulation
        numbers_of_variables: list
            the numbers of variables in the subdomains
        n: int
            the number of mesh cells
        boundary_interfaces: list
            the boundary interfaces between different-order subdomains
        tols_coarsening: list
            the threshold values for model coarsening
        tols_refinement: list
            the threshold values for model refinement
        increase_criterion_flags: np.ndarray
            the flags for model refinement
        number_of_coarsening_estimators_pde: int
            the number of estimators for model coarsening that are computed

        Returns
        -------
        flags_refinement: np.ndarray
            numpy 1D array, containing the cellwise flags for model refinement (1 if flagged for model coarsening)
        flags_coarsening: np.ndarray
            numpy 1D array, containing the cellwise flags for model coarsening (-1 if flagged for model coarsening)
        estimators_refinement: np.ndarray
            numpy 2D array containing the values for all the estimators for model refinement, in each mesh cell
        estimators_coarsening: np.ndarray
            numpy 2D array containing the values for all the estimators for model coarsening, in each mesh cell
            
        """
        
        pass

class OneStepModelErrorApproximation(DomainDecomposition1D):

    """
    This class represents a domain decomposition strategy that updates the orders of the 
    moment models in each mesh cell, for 1D moment models, and using only one evaluation in time (the current discrete time).

    ...

    Attributes
    ----------
    pde_type : PDE
        the PDE (moment model) that is simulated

    Inherited methods from abstract parent class DomainDecomposition1D
    ------------------------------------------------------------------
    def __init__(self,pde_type):
        initializes the DomainDecomposition1D object
    def compute_coarsening_criterion(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,n,
            boundary_interfaces,tols_coarsening,increase_criterion_flags,number_of_coarsening_estimators_pde)
        computes the values of the coarsening estimator, and the flags for model coarsening
    def compute_refinement_criterion(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,n,
            boundary_interfaces,tols_coarsening,number_of_coarsening_estimators_pde)
        computes the values of the refinement estimator, and the flags for model refinement

    Implemented methods from abstract parent class DomainDecomposition1D
    ----------------
    def compute_domain_decomposition(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,
                            mesh_resolution,boundary_interfaces_discretized,tols_coarsening,tols_refinement,
                            number_of_coarsening_estimators_pde,number_of_refinement_estimators_pde)
        computes the updated model orders in each mesh cell
  
    """ 

    def compute_domain_decomposition(self,
                                     values: np.ndarray,
                                     delta_t: float,
                                     delta_x: float,
                                     orders: list,
                                     max_order: int,
                                     numbers_of_variables: list,
                                     mesh_resolution: int,
                                     boundary_interfaces_discretized: list,
                                     tols_coarsening: list,
                                     tols_refinement: list,
                                     number_of_coarsening_estimators_pde: int,
                                     number_of_refinement_estimators_pde: int) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:

        estimators_refinement,flags_refinement = self.compute_refinement_criterion(
                                                            values,
                                                            delta_t,
                                                            delta_x,
                                                            orders,
                                                            max_order,
                                                            numbers_of_variables,
                                                            mesh_resolution,
                                                            boundary_interfaces_discretized,
                                                            tols_refinement,
                                                            number_of_refinement_estimators_pde)

        estimators_coarsening,flags_coarsening = self.compute_coarsening_criterion(
                                                            values,
                                                            delta_t,
                                                            delta_x,
                                                            orders,
                                                            max_order,
                                                            numbers_of_variables,
                                                            mesh_resolution,
                                                            boundary_interfaces_discretized,
                                                            tols_coarsening,
                                                            flags_refinement,
                                                            number_of_coarsening_estimators_pde)

        return flags_refinement,flags_coarsening,estimators_refinement,estimators_coarsening    

class TwoStepModelErrorApproximation(DomainDecomposition1D):

    """
    This class represents a domain decomposition strategy that updates the orders of the moment models 
    in each mesh cell, for 1D moment models, and using 2 evaluations in time (the current discrete time
    and the subsequent discrete times, which includes predictions).

    ...

    Attributes
    ----------
    spatial_discretization_predictor: SpatialDiscretization
        the spatial discretization that is used to compute the predictions
    time_integrator: TimeIntegration
        the time integrator that is used to compute the predictions
    pde_type : PDE
        the PDE (moment model) that is simulated
    mesh_resolution: int
        the number of mesh cells
    max_order: int
        the maximum order in the adaptive simulation
    max_number_of_variables: int
        the maximum number of variables in the adaptive simulation

    Inherited methods from abstract parent class DomainDecomposition1D
    ------------------------------------------------------------------
    def compute_coarsening_criterion(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,n,
            boundary_interfaces,tols_coarsening,increase_criterion_flags,number_of_coarsening_estimators_pde)
        computes the values of the coarsening estimator, and the flags for model coarsening
    def compute_refinement_criterion(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,n,
            boundary_interfaces,tols_coarsening,number_of_coarsening_estimators_pde)
        computes the values of the refinement estimator, and the flags for model refinement

    Implemented methods from abstract parent class DomainDecomposition1D
    -------------------------------------------------------------------
    def compute_domain_decomposition(self,values,delta_t,delta_x,orders,max_order,numbers_of_variables,
                            mesh_resolution,boundary_interfaces_discretized,tols_coarsening,tols_refinement,
                            number_of_coarsening_estimators_pde,number_of_refinement_estimators_pde)
        computes the updated model orders in each mesh cell

    Overriden methods from abstract parent class DomainDecomposition1D
    def __init__(self,spatial_discretization_predictor,time_integrator,pde_type,mesh_resolution,max_order,max_number_of_variables):
        initializes the TwoStepModelErrorApproximation object

    Instance methods
    ----------------
    def predict_values(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        predicts the state variables values in each mesh cell at the next discrete time
    def _transport_step_augmented(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        performs one time step of the transport part for computing the predictions
    def _source_step_augmented(self,values,delta_t,delta_x,orders,numbers_of_variables,boundary_interfaces_discretized)
        performs one time step of the source part for computing the predictions
    """ 

    def __init__(self,
             spatial_discretization_predictor: spatialDiscretization.SpatialDiscretization,
             time_integrator: timeIntegration.TimeIntegration,
             pde_type: pde.PDE,
             mesh_resolution: int,
             max_order: int,
             max_number_of_variables: int):

        """
        Initializes the TwoStepModelErrorApproximation object and constructs all the attributes.

        Parameters
        ----------
        spatial_discretization_predictor: SpatialDiscretization
            the spatial discretization that is used to compute the predictions
        time_integrator: TimeIntegration
            the time integrator that is used to compute the predictions
        pde_type : PDE
            the PDE (moment model) that is simulated
        mesh_resolution: int
            the number of mesh cells
        max_order: int
            the maximum order in the adaptive simulation
        max_number_of_variables: int
            the maximum number of variables in the adaptive simulation

        """

        self.spatial_discretization_predictor = spatial_discretization_predictor
        self.time_integrator = time_integrator
        self.pde_type = pde_type
        self.mesh_resolution = mesh_resolution
        self.max_order = max_order
        self.max_number_of_variables = max_number_of_variables

    def compute_domain_decomposition(self,
                                     values: np.ndarray,
                                     delta_t: float,
                                     delta_x: float,
                                     orders: list,
                                     max_order: int,
                                     numbers_of_variables: list,
                                     mesh_resolution: int,
                                     boundary_interfaces_discretized: list,
                                     tols_coarsening: list,
                                     tols_refinement: list,
                                     number_of_coarsening_estimators_pde: int,
                                     number_of_refinement_estimators_pde: int) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
        
        prev_values = np.copy(values)

        predicted_values = self.predict_values(values,
                                               delta_t,
                                               delta_x,
                                               orders,
                                               numbers_of_variables,
                                               boundary_interfaces_discretized)
        
        estimators_refinement_prev,flags_refinement_prev = self.compute_refinement_criterion(
                                                            prev_values,
                                                            delta_t,
                                                            delta_x,
                                                            orders,
                                                            max_order,
                                                            numbers_of_variables,
                                                            mesh_resolution,
                                                            boundary_interfaces_discretized,
                                                            tols_refinement,
                                                            number_of_refinement_estimators_pde)

        estimators_refinement_next,flags_refinement_next = self.compute_refinement_criterion(
                                                            predicted_values,
                                                            delta_t,
                                                            delta_x,
                                                            orders,
                                                            max_order,
                                                            numbers_of_variables,
                                                            mesh_resolution,
                                                            boundary_interfaces_discretized,
                                                            tols_refinement,
                                                            number_of_refinement_estimators_pde)

        flags_refinement = np.maximum(flags_refinement_next,flags_refinement_prev)
        estimators_refinement = np.maximum(estimators_refinement_next,estimators_refinement_prev)

        estimators_coarsening_prev,flags_coarsening_prev = self.compute_coarsening_criterion(
                                                            prev_values,
                                                            delta_t,
                                                            delta_x,
                                                            orders,
                                                            max_order,
                                                            numbers_of_variables,
                                                            mesh_resolution,
                                                            boundary_interfaces_discretized,
                                                            tols_coarsening,
                                                            flags_refinement,
                                                            number_of_coarsening_estimators_pde)

        estimators_coarsening_next,flags_coarsening_next = self.compute_coarsening_criterion(
                                                            predicted_values,
                                                            delta_t,
                                                            delta_x,
                                                            orders,
                                                            max_order,
                                                            numbers_of_variables,
                                                            mesh_resolution,
                                                            boundary_interfaces_discretized,
                                                            tols_coarsening,
                                                            flags_refinement,
                                                            number_of_coarsening_estimators_pde)

        flags_coarsening = np.minimum(flags_coarsening_next,flags_coarsening_prev)
        estimators_coarsening = np.maximum(estimators_coarsening_next,estimators_coarsening_prev)

        return flags_refinement,flags_coarsening,estimators_refinement,estimators_coarsening

    def predict_values(self,
                       values: np.ndarray,
                       delta_t: float,
                       delta_x: float,
                       orders: list,
                       numbers_of_variables: list,
                       boundary_interfaces_discretized: list) -> np.ndarray:
        
        """
        Predicts the state variables values in each mesh cell at the next discrete time.

        Parameters
        ----------
        values: np.ndarray
            the state variable values at the current discrete times, in each mesh cell
        delta_t: float
            the time step size
        delta_x: float
            the mesh cell size
        orders: list
            the orders in each subdomain
        numbers_of_variables: list
            the numbers of variables in each subdomain
        boundary_interfaces_discretized: list
            the boundary interfaces between the different-order moment models

        Returns
        -------
        predicted_values: np.ndarray
            the predicted state variable values at the next discrete time, in each mesh cell

        """

        predicted_values_transport = self._transport_step_augmented(values,
                                                                    delta_t,
                                                                    delta_x,
                                                                    orders,
                                                                    numbers_of_variables,
                                                                    boundary_interfaces_discretized)
        predicted_values = self._source_step_augmented(predicted_values_transport,
                                                        delta_t,
                                                        delta_x,
                                                        orders,
                                                        numbers_of_variables,
                                                        boundary_interfaces_discretized)
        return predicted_values

    def _transport_step_augmented(self,
                       values: np.ndarray,
                       delta_t: float,
                       delta_x: float,
                       orders: list,
                       numbers_of_variables: list,
                       boundary_interfaces_discretized: list) -> np.ndarray:   

        """
        Performs one time step of the transport part for computing the predictions.

        Parameters
        ----------
        values: np.ndarray
            the state variable values at the current discrete times, in each mesh cell
        delta_t: float
            the time step size
        delta_x: float
            the mesh cell size
        orders: list
            the orders in each subdomain
        numbers_of_variables: list
            the numbers of variables in each subdomain
        boundary_interfaces_discretized: list
            the boundary interfaces between the different-order moment models

        Returns
        -------
        predicted_values_transport: np.ndarray
            the predicted state variable values after the transport step

        """

        nr_augmented_variables = 1 #Generalize this such that it includes the general case

        fluctuations_min = np.zeros((self.mesh_resolution+1,self.max_number_of_variables))
        fluctuations_plus = np.zeros((self.mesh_resolution+1,self.max_number_of_variables))

        prev_values = np.copy(values)
        predicted_values_transport = np.copy(values)

        right_boundary_subdomain = 0
        for m in range(len(boundary_interfaces_discretized)):
            order = orders[m]
            n_variables = numbers_of_variables[m]

            left_boundary_subdomain = right_boundary_subdomain + 1
            right_boundary_subdomain = boundary_interfaces_discretized[m]                 

            if order < self.max_order:
                nr_augmented_variables = 1
                def system_matrix(cell_values):
                    return self.pde_type.compute_system_matrix_augmented(order,cell_values)
            else:
                def system_matrix(cell_values):
                    return self.pde_type.compute_system_matrix(order,cell_values)
                nr_augmented_variables = 0                                             

            left_boundary_value = np.zeros(self.max_number_of_variables)
            right_boundary_value = np.zeros(self.max_number_of_variables)
            left_boundary_value[:n_variables] = prev_values[left_boundary_subdomain-1,:n_variables]
            right_boundary_value[:n_variables] = prev_values[right_boundary_subdomain+1,:n_variables]

            y,fluctuations_plus[left_boundary_subdomain-1,:n_variables+nr_augmented_variables] =\
                self.spatial_discretization_predictor.compute_fluctuation(
                                                    left_boundary_value[:n_variables+nr_augmented_variables],
                                                    values[left_boundary_subdomain,:n_variables+nr_augmented_variables],
                                                    system_matrix,
                                                    delta_t,
                                                    delta_x)  
            
            for i in range(left_boundary_subdomain+1,right_boundary_subdomain+1):
                fluctuations_min[i-1,:n_variables+nr_augmented_variables],\
                    fluctuations_plus[i-1,:n_variables+nr_augmented_variables] =\
                    self.spatial_discretization_predictor.compute_fluctuation(
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
                self.spatial_discretization_predictor.compute_fluctuation(
                                                    values[i,:n_variables+nr_augmented_variables],
                                                    right_boundary_value[:n_variables+nr_augmented_variables],
                                                    system_matrix,
                                                    delta_t,
                                                    delta_x)              
            predicted_values_transport[i,:n_variables+nr_augmented_variables] =\
                values[i,:n_variables+nr_augmented_variables]\
                - delta_t/delta_x*(fluctuations_plus[i-1,:n_variables+nr_augmented_variables]+\
                                   fluctuations_min[i,:n_variables+nr_augmented_variables]) 

        order = orders[-1]
        n_variables = numbers_of_variables[-1]

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
            self.spatial_discretization_predictor.compute_fluctuation(
                                                left_boundary_value[:n_variables+nr_augmented_variables],
                                                values[left_boundary_subdomain,:n_variables+nr_augmented_variables],
                                                system_matrix,
                                                delta_t,
                                                delta_x)  
        # fluctuations_plus[left_boundary_subdomain-1,:] = np.zeros(self.max_number_of_variables)
        for i in range(left_boundary_subdomain+1,self.mesh_resolution+2):
            fluctuations_min[i-1,:n_variables+nr_augmented_variables],\
                fluctuations_plus[i-1,:n_variables+nr_augmented_variables] =\
                self.spatial_discretization_predictor.compute_fluctuation(
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
                       values: np.ndarray,
                       delta_t: float,
                       delta_x: float,
                       orders: list,
                       numbers_of_variables: list,
                       boundary_interfaces_discretized: list) -> np.ndarray:

        """
        Performs one time step of the source part for computing the predictions.

        Parameters
        ----------
        values: np.ndarray
            the state variable values at the current discrete times, in each mesh cell
        delta_t: float
            the time step size
        delta_x: float
            the mesh cell size
        orders: list
            the orders in each subdomain
        numbers_of_variables: list
            the numbers of variables in each subdomain
        boundary_interfaces_discretized: list
            the boundary interfaces between the different-order moment models

        Returns
        -------
        predicted_values: np.ndarray
            the predicted state variable values after the source step

        """

        predicted_values = np.copy(values)

        nr_augmented_variables = 1 #Generalize this such that it includes the general case

        right_boundary_subdomain = 0

        for m in range(len(boundary_interfaces_discretized)):
            order = orders[m]
            n_variables = numbers_of_variables[m]

            if order < self.max_order:
                nr_augmented_variables = 1
                def source_term(cell_values,delta_t):
                    return self.pde_type.compute_source_term_augmented(order,cell_values,delta_t)
            else:
                def source_term(cell_values,delta_t):
                    return self.pde_type.compute_source_term(order,cell_values,delta_t)        
                nr_augmented_variables = 0       

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = boundary_interfaces_discretized[m]

            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                predicted_values[i,:n_variables+nr_augmented_variables] \
                    = self.time_integrator.integrate(values[i,:n_variables+nr_augmented_variables],source_term,delta_t)

        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        if order < self.max_order:
            nr_augmented_variables = 1
            def source_term(cell_values,delta_t):
                return self.pde_type.compute_source_term_augmented(order,cell_values,delta_t)
        else:
            def source_term(cell_values,delta_t):
                return self.pde_type.compute_source_term(order,cell_values,delta_t)        
            nr_augmented_variables = 0         
                    
        for i in range(right_boundary_subdomain+1,self.mesh_resolution+1): 
            predicted_values[i,:n_variables+nr_augmented_variables] \
                = self.time_integrator.integrate(values[i,:n_variables+nr_augmented_variables],source_term,delta_t)
            
        return predicted_values 



class AdaptiveSimulationInterfaceCoupling1D(ABC):
    
    """
    TODO: implement the spatial coupling in implementations of this abstract class, 
          instead of implementing the spatial coupling directly in the simulation class.
    """
    
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def compute_fluctuations_interface(self):
        pass

class PaddedBufferCell1D(AdaptiveSimulationInterfaceCoupling1D,ABC):

    def __init__(self):
        pass

    def construct_padded_vector(self):
        pass

    def compute_fluctuations_interface(self):
        pass



class SubdomainReconstruction1D(ABC):
    
    """
    This abstract class represents a subdomain reconstruction method for a 1D domain decomposition.

    ...

    Attributes
    ----------
    boundary_condition : PDE
        the boundary condition
    reconstruct_subdomains : function
        function that takes the cellwise orders and reconstructs the subdomains and boundary interfaces

    Instance methods
    ----------------  
    def __init__(self,boundary_condition)
        initializes the SubdomainReconstruction1D object
    def update_orders_cellwise(self,orders_cellwise,numbers_of_variables_cellwise,mesh_resolution,domain_decomposition_flags)
        update the orders of the moment model in each mesh cell
    def set_removed_and_added_moments(self,orders,numbers_of_variables,boundary_interfaces,values,max_order,nvar_min_order,padding,n)
        sets values for the moments that have been removed (in case of model coarsening) or added (in case of model refinement)

    abstract methods
    ----------------
    def _reconstruct_subdomains(self)
        initializes the function that reconstructs the subdomains in the domain decomposition
    def _reconstruct_subdomains_nonPeriodicBoundary(self,values,orders_cellwise,numbers_of_variables_cellwise,
                                    max_order,mesh_resolution,padding,domain_decomposition_flags)
        reconstructs the subdomains when the boundary condition is not periodic
    def _reconstruct_subdomains_periodicBoundary(self,values,orders_cellwise,numbers_of_variables_cellwise,
                                    max_order,mesh_resolution,padding,domain_decomposition_flags)
        reconstructs the subdomains when the boundary condition is periodic    
    """ 
    
    def __init__(self,
             boundary_condition: str):
        """
        Initializes the SubdomainReconstruction1D object and constructs all the attributes.

        Parameters
        ----------
        boundary_condition : str
            the boundary condition

        """
        self.boundary_condition = boundary_condition
        self.reconstruct_subdomains = self._reconstruct_subdomains()

    @abstractmethod
    def _reconstruct_subdomains(self):
        """
        Initializes the function that reconstructs the subdomains in the domain decomposition.

        Parameters
        ----------
        None

        """
        pass

    @abstractmethod
    def _reconstruct_subdomains_nonPeriodicBoundary(self,
                               values: np.ndarray,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               max_order: int,
                               mesh_resolution: int,
                               padding: bool,
                               domain_decomposition_flags: np.ndarray) -> tuple[np.ndarray,list,list,np.ndarray,np.ndarray,list]:    
        """
        Reconstructs the subdomains when the boundary condition is not periodic.

        Parameters
        ----------
        values: np.ndarray
            the current state variable values
        orders_cellwise : np.ndarray
            the current cellwise orders
        numbers_of_variables_cellwise: np.ndarray
            the current cellwise numbers of variables
        max_order: int
            the maximum order in the adaptive simulation
        mesh_resolution: int
            the number of mesh cells
        padding: bool
            whether the values of the padded buffer cell are set to zero (false) or not (true)
        domain_decomposition_flags: np.ndarray
            the pre-computed flags for model coarsening and refinement in each mesh cell

        Returns
        -------
        values_copy: np.ndarray
            the values after the added and removed moments are handled
        orders_out: list of integers
            the orders in the subdomains
        numbers_of_variables_out: list of integers
            the numbers of variables in the subdomains
        orders_cellwise: np.ndarray
            the updated orders in each mesh cell
        numbers_of_variables_cellwise: np.ndarray
            the updated numbers of variables in each mesh cell
        boundary_interfaces: list of integers
            the updated boundary interfaces between the different-order subdomains
        """
        pass
    
    @abstractmethod
    def _reconstruct_subdomains_periodicBoundary(self,
                            values: np.ndarray,
                            orders_cellwise: np.ndarray,
                            numbers_of_variables_cellwise: np.ndarray,
                            max_order: int,
                            mesh_resolution: int,
                            padding: bool,
                            domain_decomposition_flags) -> tuple[np.ndarray,list,list,np.ndarray,np.ndarray,list]: 
        """
        Reconstructs the subdomains when the boundary condition is periodic.

        Parameters
        ----------
        values: np.ndarray
            the current state variable values
        orders_cellwise : np.ndarray
            the current cellwise orders
        numbers_of_variables_cellwise: np.ndarray
            the current cellwise numbers of variables
        max_order: int
            the maximum order in the adaptive simulation
        mesh_resolution: int
            the number of mesh cells
        padding: bool
            whether the values of the padded buffer cell are set to zero (false) or not (true)
        domain_decomposition_flags: np.ndarray
            the pre-computed flags for model coarsening and refinement in each mesh cell

        Returns
        -------
        values_copy: np.ndarray
            the values after the added and removed moments are handled
        orders_out: list of integers
            the orders in the subdomains
        numbers_of_variables_out: list of integers
            the numbers of variables in the subdomains
        orders_cellwise: np.ndarray
            the updated orders in each mesh cell
        numbers_of_variables_cellwise: np.ndarray
            the updated numbers of variables in each mesh cell
        boundary_interfaces: list of integers
            the updated boundary interfaces between the different-order subdomains
        """
        
        pass

    def update_orders_cellwise(self,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               mesh_resolution: int,
                               domain_decomposition_flags: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
        """
        Updates the orders of the moment model in each mesh cell.

        Parameters
        ----------
        orders_cellwise : np.ndarray
            the current cellwise orders
        numbers_of_variables_cellwise: np.ndarray
            the current cellwise numbers of variables
        mesh_resolution: int
            the number of mesh cells
        domain_decomposition_flags: np.ndarray
            the pre-computed flags for model coarsening and refinement in each mesh cell

        Returns
        -------
        orders_cellwise: np.ndarray
            the updated orders in each mesh cell
        numbers_of_variables_cellwise: np.ndarray
            the updated numbers of variables in each mesh cell
        """
        for i in range(mesh_resolution):
            orders_cellwise[i+1] += int(domain_decomposition_flags[i])
            numbers_of_variables_cellwise[i+1] += int(domain_decomposition_flags[i])   

        return orders_cellwise, numbers_of_variables_cellwise

    def set_removed_and_added_moments(self,
                                    orders: list,
                                    numbers_of_variables: list,
                                    boundary_interfaces: list,
                                    values: np.ndarray,
                                    max_order: int,
                                    nvar_min_order: int,
                                    padding: bool) -> np.ndarray:
        """
        Sets values for the moments that have been removed (in case of model coarsening) or
        added (in case of model refinement).

        Parameters
        ----------
        orders: list of integers
            the orders in the subdomains
        numbers_of_variables: list
            the numbers of variables in the subdomains
        boundary_interfaces: list
            the boundary interfaces between different-order subdomains
        values: np.ndarray
            the current state variable values
        max_order: int
            the maximum order in the adaptive simulation
        nvar_min_order: int
            the difference between the number of variables and the order of the moment model
        padding: bool
            whether the values of the padded buffer cell are set to zero (false) or not (true)

        Returns
        -------
        values: np.ndarray
            the state variable values after handling the removed and added moments, in each mesh cell
        """
        
        # Set undefined moments and padded moments to zero
        right_boundary = -1
        for m in range(len(boundary_interfaces)):
            left_boundary = right_boundary+1
            right_boundary = boundary_interfaces[m]
            values[left_boundary:right_boundary+1,numbers_of_variables[m]:max_order+nvar_min_order]=0
            if padding:
                if orders[m] < orders[m+1]:
                    values[right_boundary,numbers_of_variables[m]:] = \
                        values[right_boundary+1,numbers_of_variables[m]:] # update boundary interface boundary condition 
                else:
                    values[right_boundary+1,numbers_of_variables[m+1]:] = \
                        values[right_boundary,numbers_of_variables[m+1]:] # update boundary interface boundary condition 
        values[right_boundary+1,numbers_of_variables[-1]:max_order+nvar_min_order] = 0
        if padding:
            if orders[-2] < orders[-1]:
                values[right_boundary,numbers_of_variables[-2]:] = \
                    values[right_boundary+1,numbers_of_variables[-2]:] # update boundary interface boundary condition 
            else:
                values[right_boundary+1,numbers_of_variables[-1]:] = \
                    values[right_boundary,numbers_of_variables[-1]:] # update boundary interface boundary condition 

        return values

class CellwiseSubdomainReconstruction1D(SubdomainReconstruction1D):

    """
    This class represents a cellwise subdomain reconstruction method for a 1D domain decomposition.

    ...

    Attributes
    ----------
    boundary_condition : PDE
        the boundary condition
    reconstruct_subdomains : function
        function that takes the cellwise orders and reconstructs the subdomains and boundary interfaces

    Methods inherited from abstract parent class SubdomainReconstruction1D
    -----------------------------------------------------------------------
    def __init__(self,boundary_condition)
        initializes the CellwiseSubdomainReconstruction1D
    def update_orders_cellwise(self,orders_cellwise,numbers_of_variables_cellwise,mesh_resolution,domain_decomposition_flags)
        update the orders of the moment model in each mesh cell
    def set_removed_and_added_moments(self,orders,numbers_of_variables,boundary_interfaces,values,max_order,nvar_min_order,padding,n)
        sets values for the moments that have been removed (in case of model coarsening) or added (in case of model refinement)

    Implemented methods from abstract parent class SubdomainReconstruction1D
    ------------------------------------------------------------------------
    def _reconstruct_subdomains(self)
        initializes the function that reconstructs the subdomains in the domain decomposition
    def _reconstruct_subdomains_nonPeriodicBoundary(self,values,orders_cellwise,numbers_of_variables_cellwise,
                                    max_order,mesh_resolution,padding,domain_decomposition_flags)
        reconstructs the subdomains when the boundary condition is not periodic
    def _reconstruct_subdomains_periodicBoundary(self,values,orders_cellwise,numbers_of_variables_cellwise,
                                    max_order,mesh_resolution,padding,domain_decomposition_flags)
        reconstructs the subdomains when the boundary condition is periodic

    Instance methods
    ----------------
    def compute_boundary_interfaces(self,orders_cellwise,numbers_of_variables_cellwise,mesh_resolution)
        computes the boundary interfaces between different-order subdomains
    """ 

    def _reconstruct_subdomains(self) -> np.ndarray:
        
        _reconstruct_subdomains_fun = self._reconstruct_subdomains_periodicBoundary if self.boundary_condition == 'PERIODIC'\
            else self._reconstruct_subdomains_nonPeriodicBoundary

        return _reconstruct_subdomains_fun

    def _reconstruct_subdomains_nonPeriodicBoundary(self,
                               values: np.ndarray,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               max_order: int,
                               mesh_resolution: int,
                               padding: bool,
                               domain_decomposition_flags: np.ndarray) -> tuple[np.ndarray,list,list,np.ndarray,np.ndarray,list]:

        orders_cellwise, numbers_of_variables_cellwise = self.update_orders_cellwise(orders_cellwise,
                                                                                     numbers_of_variables_cellwise,
                                                                                     mesh_resolution,
                                                                                     domain_decomposition_flags)
        orders_cellwise[0] = orders_cellwise[1]
        numbers_of_variables_cellwise[0] = numbers_of_variables_cellwise[1]
        orders_cellwise[-1] = orders_cellwise[-2]
        numbers_of_variables_cellwise[-1] = numbers_of_variables_cellwise[-2]

        boundary_interfaces, orders_out, numbers_of_variables_out = self.compute_boundary_interfaces(orders_cellwise,
                                                                                                numbers_of_variables_cellwise,
                                                                                                mesh_resolution)

        values_copy = np.copy(values)

        values_copy = self.set_removed_and_added_moments(orders_out,
                                                  numbers_of_variables_out,
                                                  boundary_interfaces,
                                                  values,
                                                  max_order,
                                                  numbers_of_variables_cellwise[0]-orders_cellwise[0],
                                                  padding)

        if len(boundary_interfaces) == 0:
            orders_out.append(int(orders_cellwise[0]))
            numbers_of_variables_out.append(int(numbers_of_variables_cellwise[0]))
            boundary_interfaces.append(np.floor_divide(mesh_resolution,2))

        return values_copy,orders_out,numbers_of_variables_out,orders_cellwise,numbers_of_variables_cellwise,boundary_interfaces

    def _reconstruct_subdomains_periodicBoundary(self,
                               values: np.ndarray,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               max_order: int,
                               mesh_resolution: int,
                               padding: bool,
                               domain_decomposition_flags: np.ndarray) -> tuple[np.ndarray,list,list,np.ndarray,np.ndarray,list]:

        orders_cellwise, numbers_of_variables_cellwise = self.update_orders_cellwise(orders_cellwise,
                                                                                     numbers_of_variables_cellwise,
                                                                                     mesh_resolution,
                                                                                     domain_decomposition_flags)
        
        max_order_boundary = max(orders_cellwise[-2],orders_cellwise[1])
        max_number_of_variables_boundary = max(numbers_of_variables_cellwise[-2],numbers_of_variables_cellwise[1])
        orders_cellwise[-2] = max_order_boundary
        orders_cellwise[-1] = max_order_boundary
        orders_cellwise[0] = max_order_boundary
        orders_cellwise[1] = max_order_boundary
        numbers_of_variables_cellwise[-2] = max_number_of_variables_boundary
        numbers_of_variables_cellwise[-1] = max_number_of_variables_boundary
        numbers_of_variables_cellwise[0] = max_number_of_variables_boundary
        numbers_of_variables_cellwise[1] = max_number_of_variables_boundary

        boundary_interfaces, orders_out, numbers_of_variables_out = self.compute_boundary_interfaces(orders_cellwise,
                                                                                                numbers_of_variables_cellwise,
                                                                                                mesh_resolution)

        values_copy = np.copy(values)

        values_copy = self.set_removed_and_added_moments(orders_out,
                                                  numbers_of_variables_out,
                                                  boundary_interfaces,
                                                  values,
                                                  max_order,
                                                  numbers_of_variables_cellwise[0]-orders_cellwise[0],
                                                  padding)

        if len(boundary_interfaces) == 0:
            orders_out.append(int(orders_cellwise[0]))
            orders_out.append(int(orders_cellwise[0]))
            numbers_of_variables_out.append(int(numbers_of_variables_cellwise[0]))
            numbers_of_variables_out.append(int(numbers_of_variables_cellwise[0]))
            boundary_interfaces.append(np.floor_divide(mesh_resolution,3))
            boundary_interfaces.append(2*np.floor_divide(mesh_resolution,3))

        return values_copy,orders_out,numbers_of_variables_out,orders_cellwise,numbers_of_variables_cellwise, boundary_interfaces

    def compute_boundary_interfaces(self,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               mesh_resolution: int) -> tuple[list,list,list]:
        """
        Computes the positions of the boundary interfaces between different-order moment models.

        Parameters
        ----------
        orders_cellwise: np.ndarray
            the updated orders in each mesh cell
        numbers_of_variables_cellwise: np.ndarray
            the updated numbers of variables in each mesh cell
        mesh_resolution: int
            the number of mesh cells

        Returns
        -------
        boundary_interfaces: list of integers
            the boundary interfaces between the different-order subdomains
        orders_out: list of integers
            the orders in the subdomains
        numbers_of_variables_out: list of integers
            the numbers of variables in the subdomains
        """
        boundary_interfaces = []
        orders_out = []
        numbers_of_variables_out = []

        orders_out.append(int(orders_cellwise[0]))
        numbers_of_variables_out.append(int(numbers_of_variables_cellwise[0]))
        for i in range(1,mesh_resolution+1):
            if orders_cellwise[i] != orders_cellwise[i+1]:
                orders_out.append(int(orders_cellwise[i+1]))
                numbers_of_variables_out.append((int(numbers_of_variables_cellwise[i+1])))
                boundary_interfaces.append(int(i))

        return boundary_interfaces, orders_out, numbers_of_variables_out

class SmoothedSubdomainReconstruction1D(SubdomainReconstruction1D):

    """
    This class represents a smoothed subdomain reconstruction method for a 1D domain decomposition.

    ...

    Attributes
    ----------
    boundary_condition : PDE
        the boundary condition
    reconstruct_subdomains : function
        function that takes the cellwise orders and reconstructs the subdomains and boundary interfaces
    smooth_par: int
        the maximum number of different-order subdomains
    orders_subdomains : np.ndarray
        the orders in the fixed subdomains (i.e., the groups of cells that take the same order)
    increase_flags_subdomains : np.ndarray
        the current model refinement flags in each fixed subdomain
    breakdown_criteria_flags_subdomains : np.ndarray
        the flags for model refinement and model coarsening in each subdomain
    numbers_of_variables_subdomains : np.ndarray
        the numbers of variables in the fixed subdomains (i.e., the groups of cells that take the same order)
    boundary_interfaces_subdomains : np.ndarray 
        the fixed boundary interfaces between the fixed subdomains
    nvar_min_order : int
        the difference between the number of variables and the order of the moment model
    n_cells_subdomain : int
        the numbers of cells in the fixed subdomains (at least, in the interior ones)
    subdomain_start : int
        the start of the second subdomain
    interpolation : bool
        whether the subdomains with added moments (after model refinement) are interpolated or not
    interpolate_subdomains : function
        function that interpolates the added moments

    Methods inherited from abstract parent class SubdomainReconstruction1D
    -----------------------------------------------------------------------
    def update_orders_cellwise(self,orders_cellwise,numbers_of_variables_cellwise,mesh_resolution,domain_decomposition_flags)
        update the orders of the moment model in each mesh cell
    def set_removed_and_added_moments(self,orders,numbers_of_variables,boundary_interfaces,values,max_order,nvar_min_order,padding,n)
        sets values for the moments that have been removed (in case of model coarsening) or added (in case of model refinement)

    Implemented methods from abstract parent class SubdomainReconstruction1D
    ------------------------------------------------------------------------
    def _reconstruct_subdomains(self)
        initializes the function that reconstructs the subdomains in the domain decomposition
    def _reconstruct_subdomains_nonPeriodicBoundary(self,values,orders_cellwise,numbers_of_variables_cellwise,
                                    max_order,mesh_resolution,padding,domain_decomposition_flags)
        reconstructs the subdomains when the boundary condition is not periodic
    def _reconstruct_subdomains_periodicBoundary(self,values,orders_cellwise,numbers_of_variables_cellwise,
                                    max_order,mesh_resolution,padding,domain_decomposition_flags)
        reconstructs the subdomains when the boundary condition is periodic

    Instance methods
    ----------------
    def compute_boundary_interfaces(self,orders_cellwise,numbers_of_variables_cellwise,mesh_resolution)
        computes the boundary interfaces between different-order subdomains
    def _interpolate_subdomains(self)
        interpolates the added moments after model refinement 
    def _interpolate_subdomains_interior(self,values)
        interpolates the added moments after model refinement in the interior subdomains
    def _interpolate_subdomains_periodicBoundary(self,values,mesh_resolution)
        interpolates the added moments for a periodic boundary
    def _interpolate_subdomains_nonPeriodicBoundary(self,values,mesh_resolution)
        interpolates the added moments for a non-periodic boundary
    def _linear_interpolate(self,value_left,value_right,delta_i)
        linearly interpolates between a left index and a right index
    def compute_boundary_interfaces_nonPeriodic(self,orders_cellwise,numbers_of_variables_cellwise,mesh_resolution,domain_decomposition_flags)
        computes the boundary interfaces between different-order subdomains for a non-periodic boundary
    def compute_boundary_interfaces_periodic(self,orders_cellwise,numbers_of_variables_cellwise,mesh_resolution,domain_decomposition_flags)
        computes the boundary interfaces between different-order subdomains for a periodic boundary
    """ 

    def __init__(self,
             boundary_condition: str,
             smooth_par: int,
             start_order : int,
             start_nr_of_variables: int,
             mesh_resolution: int,
             interpolation: bool):
            
            self.smooth_par = smooth_par
            self.orders_subdomains = np.full(shape=self.smooth_par+1,fill_value=start_order,dtype=int)
            self.increase_flags_subdomains = np.full(shape=self.smooth_par+1,fill_value=0,dtype=int) 
            self.breakdown_criteria_flags_subdomains = np.full(shape=self.smooth_par+1,fill_value=0,dtype=int)
            # self.breakdown_criteria_flags = np.full(shape=mesh_resolution,dtype=int,fill_value=0)
            self.numbers_of_variables_subdomains = np.full(shape=self.smooth_par+1,fill_value=start_nr_of_variables,dtype=int)
            self.boundary_interfaces_subdomains = np.zeros(self.smooth_par,dtype=int)  

            self.nvar_min_order = start_nr_of_variables - start_order

            self.n_cells_subdomain = int(np.floor(mesh_resolution/self.smooth_par))
            self.subdomain_start = int(np.ceil(self.n_cells_subdomain/2))+1
            self.boundary_interfaces_subdomains[0] = self.subdomain_start - 1
            for i in range(1,self.smooth_par):
                self.boundary_interfaces_subdomains[i] = self.subdomain_start + i*self.n_cells_subdomain - 1

            super().__init__(boundary_condition)
            if self.boundary_condition != 'PERIODIC':
                self.increase_flags_subdomains = np.full(shape=self.smooth_par,fill_value=0,dtype=int)  
                self.orders_subdomains = np.full(shape=self.smooth_par,fill_value=start_order,dtype=int)
                self.increase_flags_subdomains = np.full(shape=self.smooth_par,fill_value=0,dtype=int)
                self.numbers_of_variables_subdomains = np.full(shape=self.smooth_par,fill_value=start_nr_of_variables,dtype=int)
                self.boundary_interfaces_subdomains = np.zeros(self.smooth_par-1,dtype=int)
                for i in range(self.smooth_par-1):
                    self.boundary_interfaces_subdomains[i] = (i+1)*self.n_cells_subdomain - 1 

            self.interpolation = interpolation

            self.interpolate_subdomains = self._interpolate_subdomains()

    def _reconstruct_subdomains(self) -> np.ndarray:
        
        _reconstruct_subdomains_fun = self._reconstruct_subdomains_periodicBoundary if self.boundary_condition == 'PERIODIC'\
            else self._reconstruct_subdomains_nonPeriodicBoundary

        return _reconstruct_subdomains_fun

    def _reconstruct_subdomains_nonPeriodicBoundary(self,
                               values: np.ndarray,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               max_order: int,
                               mesh_resolution: int,
                               padding: bool,
                               domain_decomposition_flags: np.ndarray) -> tuple[np.ndarray,list,list,np.ndarray,np.ndarray,list]:

        print("Current orders subdomains: ",self.orders_subdomains)

        orders_cellwise, numbers_of_variables_cellwise = self.update_orders_cellwise(orders_cellwise,
                                                                                     numbers_of_variables_cellwise,
                                                                                     mesh_resolution,
                                                                                     domain_decomposition_flags)
           
        orders_cellwise[0] = orders_cellwise[1]
        numbers_of_variables_cellwise[0] = numbers_of_variables_cellwise[1]
        orders_cellwise[-1] = orders_cellwise[-2]
        numbers_of_variables_cellwise[-1] = numbers_of_variables_cellwise[-2]

        boundary_interfaces, orders_out, numbers_of_variables_out =\
            self.compute_boundary_interfaces_nonPeriodic(orders_cellwise,
                                                        numbers_of_variables_cellwise,
                                                        mesh_resolution,
                                                        domain_decomposition_flags)

        if len(boundary_interfaces) == 0:
            orders_out.append(self.orders_subdomains[1])
            orders_out.append(self.orders_subdomains[1])
            numbers_of_variables_out.append(int(self.numbers_of_variables_subdomains[1]))
            numbers_of_variables_out.append(int(self.numbers_of_variables_subdomains[1]))
            boundary_interfaces.append(np.floor_divide(mesh_resolution,3))
            boundary_interfaces.append(np.floor_divide(2*mesh_resolution,3))

        values_copy = np.copy(values)

        values_copy = self.set_removed_and_added_moments(orders_out,
                                                  numbers_of_variables_out,
                                                  boundary_interfaces,
                                                  values,
                                                  max_order,
                                                  numbers_of_variables_cellwise[0]-orders_cellwise[0],
                                                  padding)

        if self.interpolation:
            values_copy = self.interpolate_subdomains(values_copy,
                                                      mesh_resolution)
     
        return values_copy,orders_out,numbers_of_variables_out,orders_cellwise,numbers_of_variables_cellwise,boundary_interfaces

    def _reconstruct_subdomains_periodicBoundary(self,
                               values: np.ndarray,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               max_order: int,
                               mesh_resolution: int,
                               padding: bool,
                               domain_decomposition_flags: np.ndarray) -> tuple[np.ndarray,list,list,np.ndarray,np.ndarray,list]:

        print("Current orders subdomains: ",self.orders_subdomains)

        orders_cellwise, numbers_of_variables_cellwise = self.update_orders_cellwise(orders_cellwise,
                                                                                     numbers_of_variables_cellwise,
                                                                                     mesh_resolution,
                                                                                     domain_decomposition_flags)  
        
        max_order_boundary = max(orders_cellwise[-2],orders_cellwise[1])
        max_number_of_variables_boundary = max(numbers_of_variables_cellwise[-2],numbers_of_variables_cellwise[1])
        orders_cellwise[-2] = max_order_boundary
        orders_cellwise[-1] = max_order_boundary
        orders_cellwise[0] = max_order_boundary
        orders_cellwise[1] = max_order_boundary
        numbers_of_variables_cellwise[-2] = max_number_of_variables_boundary
        numbers_of_variables_cellwise[-1] = max_number_of_variables_boundary
        numbers_of_variables_cellwise[0] = max_number_of_variables_boundary
        numbers_of_variables_cellwise[1] = max_number_of_variables_boundary

        boundary_interfaces, orders_out, numbers_of_variables_out =\
            self.compute_boundary_interfaces_periodic(orders_cellwise,
                                                    numbers_of_variables_cellwise,
                                                    mesh_resolution,
                                                    domain_decomposition_flags)

        if len(boundary_interfaces) == 0:
            orders_out.append(self.orders_subdomains[1])
            orders_out.append(self.orders_subdomains[1])
            numbers_of_variables_out.append(int(self.numbers_of_variables_subdomains[1]))
            numbers_of_variables_out.append(int(self.numbers_of_variables_subdomains[1]))
            boundary_interfaces.append(np.floor_divide(mesh_resolution,3))
            boundary_interfaces.append(np.floor_divide(2*mesh_resolution,3))

        values_copy = np.copy(values)

        values_copy = self.set_removed_and_added_moments(orders_out,
                                                  numbers_of_variables_out,
                                                  boundary_interfaces,
                                                  values,
                                                  max_order,
                                                  numbers_of_variables_cellwise[0]-orders_cellwise[0],
                                                  padding)

        if self.interpolation:
            values_copy = self.interpolate_subdomains(values_copy,
                                                      mesh_resolution)

        return values_copy,orders_out,numbers_of_variables_out,orders_cellwise,numbers_of_variables_cellwise,boundary_interfaces

    def _interpolate_subdomains(self) -> np.ndarray:

        """
        Initializes the function that interpolates the added moments in the domain decomposition.

        Parameters
        ----------
        None

        """

        _interpolate_subdomains_fun = self._interpolate_subdomains_periodicBoundary if self.boundary_condition == 'PERIODIC' \
            else self._interpolate_subdomains_nonPeriodicBoundary

        return _interpolate_subdomains_fun   

    def _interpolate_subdomains_interior(self,
                                         values: np.ndarray) -> np.ndarray:
        """
        Interplates the added moments in the interior subdomains.

        Parameters
        ----------
        values: np.ndarray
            the current state variable values

        Returns
        -------
        interpolated_values_interior : np.ndarray
            the state variable values after interpolating the added moments in the interior subdomains

        """
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

    def _interpolate_subdomains_periodicBoundary(self,
                                                 values: np.ndarray,
                                                 mesh_resolution: int) -> np.ndarray:
        """
        Interplates the added moments for a periodic boundary.

        Parameters
        ----------
        values: np.ndarray
            the current state variable values
        mesh_resolution: int
            the number of mesh cells

        Returns
        -------
        interpolated_values : np.ndarray
            the state variable values after interpolating the added moments

        """
        interpolated_values = self._interpolate_subdomains_interior(values)

        increase_flag = self.increase_flags_subdomains[0]
        boundary_left = self.boundary_interfaces_subdomains[-1]
        boundary_right = self.boundary_interfaces_subdomains[0]+1
        nr_variables = self.numbers_of_variables_subdomains[0]
        if increase_flag > 0:  
            interp_val_bound = self._linear_interpolate(
                values[boundary_left,nr_variables-increase_flag:nr_variables],
                values[boundary_right,nr_variables-increase_flag:nr_variables],
                mesh_resolution-boundary_left+boundary_right+1)
            interpolated_values[boundary_left+1:-1,nr_variables-increase_flag:nr_variables] = interp_val_bound[:mesh_resolution-boundary_left,:]
            interpolated_values[1:boundary_right,nr_variables-increase_flag:nr_variables] = interp_val_bound[mesh_resolution-boundary_left+1:,:]

        return interpolated_values            

    def _interpolate_subdomains_nonPeriodicBoundary(self,
                                                 values: np.ndarray,
                                                 mesh_resolution: int) -> np.ndarray:
        """
        Interplates the added moments for a non-periodic boundary.

        Parameters
        ----------
        values: np.ndarray
            the current state variable values
        mesh_resolution: int
            the number of mesh cells

        Returns
        -------
        interpolated_values : np.ndarray
            the state variable values after interpolating the added moments

        """

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
                                             mesh_resolution-boundary_left+3)

        return interpolated_values                     

    def _linear_interpolate(self,
                            value_left: np.ndarray,
                            value_right: np.ndarray,
                            delta_i: int) -> np.ndarray:

        """
        Interplates the added moments for a periodic boundary.

        Parameters
        ----------
        value_left: np.ndarray
            the state variables values at the left boundary of the subdomain in which the interpolation is performed
        value_right: np.ndarray
            the state variables values at the right boundary of the subdomain in which the interpolation is performed
        delta_i: int
            the difference between the index of the right boundary cell and the index of the left boundary cell

        Returns
        -------
        interpolated_values : np.ndarray
            the state variable values after interpolating the added moments between the left boundary cell 
            and the right boundary cell

        """ 

        slope = (value_right - value_left)/(delta_i)
        interpolated_values = value_left + np.outer(np.arange(1,delta_i), slope)
        # interpolated_values = np.zeros((delta_i-1,1))

        return interpolated_values

    def compute_boundary_interfaces_nonPeriodic(self,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               mesh_resolution: int,
                               domain_decomposition_flags: np.ndarray):

        """
        Computes the positions of the boundary interfaces between different-order moment models
        for a non-periodic boundary.

        Parameters
        ----------
        orders_cellwise: np.ndarray
            the updated orders in each mesh cell
        numbers_of_variables_cellwise: np.ndarray
            the updated numbers of variables in each mesh cell
        mesh_resolution: int
            the number of mesh cells
        domain_decomposition_flags: np.ndarray
            the flags for model coarsening and model refinement in each cell

        Returns
        -------
        boundary_interfaces: list of integers
            the boundary interfaces between the different-order subdomains
        orders_out: list of integers
            the orders in the subdomains
        numbers_of_variables_out: list of integers
            the numbers of variables in the subdomains
        """
        
        boundary_interfaces = []
        orders_out = []
        numbers_of_variables_out = []

        for i in range(self.smooth_par-1):
            local_order = np.max(orders_cellwise[i*self.n_cells_subdomain:(i+1)*self.n_cells_subdomain])
            local_number_of_variables = local_order + self.nvar_min_order
            max_increase_flags = int(np.max(
                domain_decomposition_flags[max(i*self.n_cells_subdomain-1,0):(i+1)*self.n_cells_subdomain-1]))
            self.increase_flags_subdomains[i] = max_increase_flags
            self.orders_subdomains[i] = local_order
            self.numbers_of_variables_subdomains[i] = local_number_of_variables
            orders_cellwise[i*self.n_cells_subdomain:(i+1)*self.n_cells_subdomain] = local_order
            numbers_of_variables_cellwise[i*self.n_cells_subdomain:(i+1)*self.n_cells_subdomain] = local_number_of_variables

        local_order = np.max(orders_cellwise[(self.smooth_par-1)*self.n_cells_subdomain:])
        self.orders_subdomains[-1] = local_order
        orders_cellwise[(self.smooth_par-1)*self.n_cells_subdomain:] = local_order
        self.numbers_of_variables_subdomains[-1] = local_order + self.nvar_min_order
        numbers_of_variables_cellwise[(self.smooth_par-1)*self.n_cells_subdomain:] = self.numbers_of_variables_subdomains[-1]
        max_increase_flags = int(np.max(domain_decomposition_flags[(self.smooth_par-1)*self.n_cells_subdomain-1:]))
        self.increase_flags_subdomains[-1] = max_increase_flags

        orders_out.append(int(self.orders_subdomains[0]))
        numbers_of_variables_out.append(int(self.numbers_of_variables_subdomains[0]))
        for i in range(self.smooth_par-1):
            if self.orders_subdomains[i] != self.orders_subdomains[i+1]:
                orders_out.append(int(self.orders_subdomains[i+1]))
                numbers_of_variables_out.append(int(self.numbers_of_variables_subdomains[i+1]))
                boundary_interfaces.append(self.n_cells_subdomain*(i+1)-1)

        return boundary_interfaces, orders_out, numbers_of_variables_out

    def compute_boundary_interfaces_periodic(self,
                               orders_cellwise: np.ndarray,
                               numbers_of_variables_cellwise: np.ndarray,
                               mesh_resolution: int,
                               domain_decomposition_flags: np.ndarray):

        """
        Computes the positions of the boundary interfaces between different-order moment models
        for a periodic boundary.

        Parameters
        ----------
        orders_cellwise: np.ndarray
            the updated orders in each mesh cell
        numbers_of_variables_cellwise: np.ndarray
            the updated numbers of variables in each mesh cell
        mesh_resolution: int
            the number of mesh cells
        domain_decomposition_flags: np.ndarray
            the flags for model coarsening and model refinement in each cell

        Returns
        -------
        boundary_interfaces: list of integers
            the boundary interfaces between the different-order subdomains
        orders_out: list of integers
            the orders in the subdomains
        numbers_of_variables_out: list of integers
            the numbers of variables in the subdomains
        """

        boundary_interfaces = []
        orders_out = []
        numbers_of_variables_out = []

        local_order_start = np.max(orders_cellwise[:self.subdomain_start])
        max_increase_flags_start = int(np.max(domain_decomposition_flags[:self.subdomain_start]))
        for i in range(self.smooth_par-1):
            local_order = np.max(orders_cellwise[self.subdomain_start+i*self.n_cells_subdomain:\
                                                      self.subdomain_start+(i+1)*self.n_cells_subdomain])
            local_number_of_variables = local_order + self.nvar_min_order
            max_increase_flags = int(np.max(domain_decomposition_flags[\
                self.subdomain_start+i*self.n_cells_subdomain:self.subdomain_start+(i+1)*self.n_cells_subdomain]))
            self.increase_flags_subdomains[i+1] = max_increase_flags
            self.orders_subdomains[i+1] = local_order
            self.numbers_of_variables_subdomains[i+1] = local_number_of_variables
            orders_cellwise[self.subdomain_start+i*self.n_cells_subdomain:\
                                    self.subdomain_start+(i+1)*self.n_cells_subdomain] = local_order
            numbers_of_variables_cellwise[self.subdomain_start+i*self.n_cells_subdomain:\
                                                    self.subdomain_start+(i+1)*self.n_cells_subdomain] = local_number_of_variables
        local_order_end = np.max(orders_cellwise[self.subdomain_start+(self.smooth_par-1)*self.n_cells_subdomain:])
        max_boundary_order = max(local_order_start,local_order_end)
        max_increase_flags_end = int(np.max(domain_decomposition_flags[self.subdomain_start+(self.smooth_par-1)*self.n_cells_subdomain:]))
        self.increase_flags_subdomains[0] = max(max_increase_flags_start,max_increase_flags_end)
        self.increase_flags_subdomains[-1] = max(max_increase_flags_start,max_increase_flags_end)
        self.orders_subdomains[0] = max_boundary_order
        self.orders_subdomains[-1] = max_boundary_order
        self.numbers_of_variables_subdomains[0] = max_boundary_order + self.nvar_min_order
        self.numbers_of_variables_subdomains[-1] = max_boundary_order + self.nvar_min_order

        orders_out.append(int(self.orders_subdomains[0]))
        numbers_of_variables_out.append(int(self.numbers_of_variables_subdomains[0]))
        for i in range(self.smooth_par):
            if self.orders_subdomains[i] != self.orders_subdomains[i+1]:
                orders_out.append(int(self.orders_subdomains[i+1]))
                numbers_of_variables_out.append(int(self.numbers_of_variables_subdomains[i+1]))
                boundary_interfaces.append(self.subdomain_start+self.n_cells_subdomain*i-1)

        return boundary_interfaces, orders_out, numbers_of_variables_out



#TODO: make this a child of the ClassicalSimulation1D
class Micro_macro(Simulation):
    """
    This class represents a micro-macro simulation.

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
                 pde_type: pde.PDE,
                 mesh: mesh.RectangularMesh,
                 boundary_condition: str,
                 initial_condition: str,
                 spatial_discretization: spatialDiscretization.SpatialDiscretization,
                 time_integration: timeIntegration.TimeIntegration):
 
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
                       g = 1) -> np.ndarray:

        delta_x = (self.mesh.boundaries[1] - self.mesh.boundaries[0])/self.mesh.resolution

        micro_moments = self._get_initial_conditions(self.mesh.cell_center_positions)
        macro_moments = np.zeros((self.mesh.resolution+2, self.macro_order+2))

        CFL = 0.5

        def micro_system_matrix(cell_values):
            return self.pde_type.compute_system_matrix(self.micro_order,cell_values)

        def micro_source_term(cell_values):
            return self.pde_type.compute_source_term(self.micro_order,cell_values)
        
        def macro_system_matrix(cell_values):
            return self.pde_type.compute_system_matrix(self.macro_order, cell_values)

        def macro_source_term(cell_values):
            return self.pde_type.compute_source_term(self.macro_order, cell_values)

        t = 0
        step = 0

        while t < t_end:

            # MICRO STEP
            micro_moments[0,:] = self._update_boundary_conditions(micro_moments,'left')
            micro_moments[-1,:] = self._update_boundary_conditions(micro_moments,'right')
            
            # Calculate step size using CFL condition
            wave_speed_sqrt = micro_moments[:,0]*int(g)
            for i in range(self.micro_order):
                wave_speed_sqrt += np.divide(micro_moments[:,i+2]*micro_moments[:,i+2],micro_moments[:,0]*micro_moments[:,0])
            max_speed =  np.max(np.abs(np.divide(micro_moments[:,1],micro_moments[:,0]))+ np.sqrt(wave_speed_sqrt))

            micro_delta_t = CFL*delta_x/max_speed

            values = np.copy(micro_moments)

            # Calculate the space derivative term
            for i in range(1,self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    values[i-1,:],
                    values[i,:],
                    micro_system_matrix,
                    'positive',
                    micro_delta_t,
                    delta_x) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    values[i,:],
                    values[i+1,:],
                    micro_system_matrix,
                    'negative',
                    micro_delta_t,
                    delta_x) 
                
                micro_moments[i,:] = values[i,:] - micro_delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
                micro_moments[i,:] = self.time_integration.integrate(micro_moments[i,:],micro_source_term,micro_delta_t)

            t += micro_delta_t

            # RESTRICTION
            macro_moments = micro_moments[:, :self.macro_order+2]

            # MACRO STEP
            macro_moments[0,:] = self._update_boundary_conditions(macro_moments,'left')
            macro_moments[-1,:] = self._update_boundary_conditions(macro_moments,'right')

            # Calculate step size using CFL condition
            wave_speed_sqrt = macro_moments[:,0]*int(g)
            for i in range(self.macro_order):
                wave_speed_sqrt += np.divide(macro_moments[:,i+2]*macro_moments[:,i+2],macro_moments[:,0]*macro_moments[:,0])
            max_speed =  np.max(np.abs(np.divide(macro_moments[:,1],macro_moments[:,0]))+np.sqrt(wave_speed_sqrt))

            macro_delta_t = CFL*delta_x/max_speed

            # Calculate the space derivative term
            values = np.copy(macro_moments)

            for i in range(1,self.mesh.resolution+1):
                fluctuation_plus = self.spatial_discretization.compute_fluctuation(
                    values[i-1,:],
                    values[i,:],
                    macro_system_matrix,
                    'positive',
                    macro_delta_t,
                    delta_x) 
                fluctuation_minus = self.spatial_discretization.compute_fluctuation(
                    values[i,:],
                    values[i+1,:],
                    macro_system_matrix,
                    'negative',
                    macro_delta_t,
                    delta_x) 

                macro_moments[i,:] = values[i,:] - macro_delta_t/delta_x*(fluctuation_plus+fluctuation_minus)
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
                               cell_centers_x: np.ndarray) -> np.ndarray:

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
                         values) -> np.ndarray:

        data_array = np.zeros((self.mesh.resolution,self.number_of_variables+1)) # rewrite this such that it can be generalized to other PDE models

        for i in range(self.mesh.resolution):
            data_array[i,0] = self.mesh.cell_center_positions[i]
        data_array[:,1] = values[1:-1,0]
        data_array[:,2] = np.divide(values[1:-1,1],data_array[:,1])
        for j in range(self.micro_order): #TODO: this is unnecessary routine here
            data_array[:,j+3] = np.divide(values[1:-1,j+2],data_array[:,1])

        return data_array
