import numpy as np

def compute_breakdown_criteria_full(self,
                                values: np.ndarray,
                                n: int,
                                delta_x: float,
                                delta_t: float,
                                max_order: int,
                                orders_cellwise: list,
                                numbers_of_variables_cellwise: list,
                                dom_decomp_val_res1: np.ndarray,
                                dom_decomp_val_res2: np.ndarray,
                                tolerance_up_source = 0.04,
                                tolerance_up_height_gradient = 0.5,
                                tolerance_up_momentum_gradient = 0.5,
                                tolerance_up_moment_gradient = 0.5,
                                tolerance_down_last_moment = 0.0001,
                                tolerance_down_res1 = 0.0001,
                                tolerance_down_res2 = 0.01) -> tuple[np.ndarray,np.ndarray]:        
    """
    Computes the breakown criteria adaptive simulation

    Parameters
    ----------
    values : list of numpy 1D arrays
        the values of the variables in each mesh cell
    n : integer
        number of grid cells
    delta_x : float
        grid cell size
    delta_t : float
        time step size
    max_order : integer
        maximum order of the model
    orders_cellwise : list of integers
        the order in each cell
    number_of_variables_cellwise : list of integers
        the number of variables in each cell
    dom_decomp_val_res1 : np.ndarray
        numpy array containing the current values of res1
    dom_decomp_val_res2 : np.ndarray
        numpy array containing the current values of res2
    tolerance_up_source : float
        threshold for the increase-criterion corresponding to the last entry of the source term vector
    tolerance_up_height_gradient : float
        threshold for the increase-criterion corresponding to the height gradient
    tolerance_up_momentum_gradient : float
        threshold for the increase-criterion corresponding to the momentum gradient       
    tolerance_up_moment_gradient : float
        threshold for the increase-criterion corresponding to the moment gradient
    tolerance_down_last_moment : float
        threshold for the decrease-criterion corresponding to the last moment
    tolerance_down_res1 : float
        threshold for the decrease-criterion corresponding to res1
    tolerance_down_res2 : float
        threshold for the decrease-criterion corresponding to res2

    Returns
    -------
    breakdown_estimators : np.ndarray
        values for each breakdown estimator in each grid cell
    breakdown_criterion_flags : np.ndarray
        flags for increasing or reducing the order in each grid cell
        this array is filled with the values of the changes in order in each grid cell
    """
    breakdown_criterion_flags = np.zeros(n)
    breakdown_estimators = np.zeros((n,max_order+4))

    for i in range(n):
        breakdown_estimators[i,0] = 1*np.abs(self.compute_source_term_lastentry(orders_cellwise[i+1],values[i+1,:numbers_of_variables_cellwise[i+1]],True))
        breakdown_estimators[i,1] = np.abs(values[i+1,numbers_of_variables_cellwise[i+1]-1]/values[i+1,0])
        breakdown_estimators[i,2] = 1*np.abs((values[i+2,0] - values[i,0]))/(2*delta_x)
        breakdown_estimators[i,3] = 1*np.abs((values[i+2,1] - values[i,1]))/(2*delta_x)
        for j in range(orders_cellwise[i+1]):
            breakdown_estimators[i,4+j] = 1*np.abs((values[i+2,2+j])-values[i,2+j])/(2*delta_x)
    breakdown_estimators[0,0] = 1*np.abs(self.compute_source_term_lastentry(orders_cellwise[1],values[1,:numbers_of_variables_cellwise[1]],True))
    breakdown_estimators[0,2] = 1*np.abs((values[2,0] - values[1,0]))/delta_x
    breakdown_estimators[0,3] = 1*np.abs((values[2,1] - values[1,1]))/delta_x
    for j in range(orders_cellwise[i+1]):
        breakdown_estimators[0,4+j] = 1*np.abs((values[2,2+j])-values[1,2+j])/delta_x    

    for i in range(n):
        if orders_cellwise[i+1] < max_order and \
            (breakdown_estimators[i,0] > tolerance_up_source or\
                breakdown_estimators[i,2] > tolerance_up_height_gradient or\
                    breakdown_estimators[i,3] > tolerance_up_momentum_gradient or\
                        np.any(breakdown_estimators[i,4:orders_cellwise[i+1]+4] > tolerance_up_moment_gradient)):
            breakdown_criterion_flags[i] = 1

        elif (breakdown_criterion_flags[i] !=1 and\
                orders_cellwise[i+1] > 0 and\
                    dom_decomp_val_res1[i] < tolerance_down_res1\
                        and dom_decomp_val_res2[i] < tolerance_down_res2\
                            and breakdown_estimators[i,1] < tolerance_down_last_moment):
            breakdown_criterion_flags[i] = -1
    
    return breakdown_estimators, breakdown_criterion_flags



def compute_refinement_criterion(self,
                                values: np.ndarray,
                                orders: list,
                                max_order: int,
                                numbers_of_variables: list,
                                n: int,
                                boundary_interfaces: list,
                                delta_x: float,
                                tolerance_increase = 0.0015) -> np.ndarray:
    
    increase_criterion_flags = np.zeros(n,dtype=int)
    breakdown_estimators_increase = np.zeros(n)

    backward_differences = np.zeros(n)
    forward_differences = np.zeros(n)

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

        backward_differences[l-1] = np.abs(self.compute_breakdown_quantity_increase(values[l,:],
                                                                                    left_boundary_value,
                                                                                    values[l,:],
                                                                                    order,
                                                                                    max_order)) 
        forward_differences[l-1] = np.abs(self.compute_breakdown_quantity_increase(values[l,:],
                                                                                    values[l,:],
                                                                                    values[l+1,:],
                                                                                    order,
                                                                                    max_order))                  
                                                
        for i in range(l+1,r):
            backward_differences[i-1] = np.abs(self.compute_breakdown_quantity_increase(values[i,:],
                                                                                        values[i-1,:],
                                                                                        values[i,:],
                                                                                        order,
                                                                                        max_order)) 
            forward_differences[i-1] = np.abs(self.compute_breakdown_quantity_increase(values[i,:],
                                                                                        values[i,:],
                                                                                        values[i+1,:],
                                                                                        order,
                                                                                        max_order))        
        backward_differences[r-1] = np.abs(self.compute_breakdown_quantity_increase(values[r,:],
                                                                                    values[r-1,:],
                                                                                    values[r,:],
                                                                                    order,
                                                                                    max_order)) 
        forward_differences[r-1] = np.abs(self.compute_breakdown_quantity_increase(values[r,:],
                                                                                    values[r,:],
                                                                                    right_boundary_value,
                                                                                    order,
                                                                                    max_order)) 

    n_variables_prev = n_variables
    order = orders[-1]
    n_variables = numbers_of_variables[-1]

    n_left = min(n_variables_prev,n_variables)

    l = r+1  

    left_boundary_value = input_values[l,:]            
    left_boundary_value[:n_left] = input_values[l-1,:n_left]  

    backward_differences[l-1] = np.abs(self.compute_breakdown_quantity_increase(values[l,:],
                                                                                left_boundary_value,
                                                                                values[l,:],
                                                                                order,
                                                                                max_order)) 
    forward_differences[l-1] = np.abs(self.compute_breakdown_quantity_increase(values[l,:],
                                                                                values[l,:],
                                                                                values[l+1,:],
                                                                                order,
                                                                                max_order))                  
                                            
    for i in range(l+1,n+1):
        backward_differences[i-1] = np.abs(self.compute_breakdown_quantity_increase(values[i,:],
                                                                                    values[i-1,:],
                                                                                    values[i,:],
                                                                                    order,
                                                                                    max_order)) 
        forward_differences[i-1] = np.abs(self.compute_breakdown_quantity_increase(values[i,:],
                                                                                    values[i,:],
                                                                                    values[i+1,:],
                                                                                    order,
                                                                                    max_order)) 

    breakdown_estimators_increase = np.maximum(forward_differences,backward_differences)/delta_x

    for i in range(n):
        if breakdown_estimators_increase[i] > tolerance_increase: 
            increase_criterion_flags[i] = 2            
    return breakdown_estimators_increase, increase_criterion_flags

def compute_coarsening_criterion(self,
                                values: np.ndarray,
                                orders: list,
                                numbers_of_variables: list,
                                n: int,
                                boundary_interfaces: list,
                                delta_x: float,
                                increase_criterion_flags: np.ndarray,
                                tolerance_decrease = 0.001) -> np.ndarray:

    decrease_criterion_flags = np.zeros(n,dtype=int)
    breakdown_estimators_decrease = np.zeros(n)

    backward_differences = np.zeros(n)
    forward_differences = np.zeros(n)

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
            backward_differences[l-1] = np.abs(self.compute_breakdown_quantity_decrease(values[l,:],
                                                                                        left_boundary_value,
                                                                                        values[l,:],
                                                                                        order)) 
            forward_differences[l-1] = np.abs(self.compute_breakdown_quantity_decrease(values[l,:],
                                                                                        values[l,:],
                                                                                        values[l+1,:],
                                                                                        order))                                                        
            for i in range(l+1,r):
                backward_differences[i-1] = np.abs(self.compute_breakdown_quantity_decrease(values[i,:],
                                                                                            values[i-1,:],
                                                                                            values[i,:],
                                                                                            order)) 
                forward_differences[i-1] = np.abs(self.compute_breakdown_quantity_decrease(values[i,:],
                                                                                            values[i,:],
                                                                                            values[i+1,:],
                                                                                            order))           
            backward_differences[r-1] = np.abs(self.compute_breakdown_quantity_decrease(values[r,:],
                                                                                        values[r-1,:],
                                                                                        values[r,:],
                                                                                        order)) 
            forward_differences[r-1] = np.abs(self.compute_breakdown_quantity_decrease(values[r,:],
                                                                                        values[r,:],
                                                                                        right_boundary_value,
                                                                                        order)) 
    n_variables_prev = n_variables
    order = orders[-1]
    n_variables = numbers_of_variables[-1]

    n_left = min(n_variables_prev,n_variables)

    l = r+1  

    left_boundary_value = input_values[l,:]            
    left_boundary_value[:n_left] = input_values[l-1,:n_left]  

    if order > 3:
        backward_differences[l-1] = np.abs(self.compute_breakdown_quantity_decrease(values[l,:],
                                                                                    left_boundary_value,
                                                                                    values[l,:],
                                                                                    order)) 
        forward_differences[l-1] = np.abs(self.compute_breakdown_quantity_decrease(values[l,:],
                                                                                    values[l,:],
                                                                                    values[l+1,:],
                                                                                    order))                                                        
        for i in range(l+1,n+1):
            backward_differences[i-1] = np.abs(self.compute_breakdown_quantity_decrease(values[i,:],
                                                                                        values[i-1,:],
                                                                                        values[i,:],
                                                                                        order)) 
            forward_differences[i-1] = np.abs(self.compute_breakdown_quantity_decrease(values[i,:],
                                                                                        values[i,:],
                                                                                        values[i+1,:],
                                                                                        order)) 

    breakdown_estimators_decrease = np.maximum(backward_differences,forward_differences)/delta_x

    for i in range(n):
        if increase_criterion_flags[i] == 0 and breakdown_estimators_decrease[i] < tolerance_decrease:
            decrease_criterion_flags[i] = -2

    return breakdown_estimators_decrease, decrease_criterion_flags

def compute_coarsening_estimator_model_difference_old(self,
                                        value_central: np.ndarray, 
                                        value_left: np.ndarray,
                                        value_right: np.ndarray,
                                        order: int,
                                        delta_t: float,
                                        delta_x: float) -> float:
    
    coarsening_quantity = 0

    if order == 4:
        coarsening_quantity = 6/value_central[0]*(value_right[3]-value_left[3])
    elif order == 5:
        coarsening_quantity = 4*value_central[3]*(value_right[1]-value_left[1])+4*(value_right[4]-value_left[4])
        # coarsening_quantity = 4*value_central[3]*(value_right[1]-value_left[1])
    else:
        coarsening_quantity = (order-1)*(value_right[order-1]-value_left[order-1])+\
                            (order-1)/2*(2*value_central[order-2]*(value_right[1]-value_left[1])+\
                            value_central[order-3]*(value_right[2]-value_left[2]))
        # coarsening_quantity = (order-1)/2*(2*value_central[order-2]*(value_right[1]-value_left[1])+\
        #                     value_central[order-3]*(value_right[2]-value_left[2]))

    coarsening_quantity = coarsening_quantity/delta_x

    return coarsening_quantity


def compute_refinement_estimator_model_difference_old(self,
                                value_central: np.ndarray, 
                                value_left: np.ndarray,
                                value_right: np.ndarray,
                                order,
                                max_order,
                                delta_t: float,
                                delta_x: float):
    
    refinement_quantity = 0
    if order == 2:
        refinement_quantity = 6/value_central[0]*(value_right[3]-value_left[3])
    elif order == 3:
        refinement_quantity = 4*value_central[3]*(value_right[1]-value_left[1])+4*(value_right[4]-value_left[4])
    elif order == max_order:
        refinement_quantity = (order+1)/2*(2*value_central[order]*(value_right[1]-value_left[1])+\
                            value_central[order-1]*(value_right[2]-value_left[2]))  
    else:
        refinement_quantity = (order+1)*(value_right[order+1]-value_left[order+1])+\
                            (order+1)/2*(2*value_central[order]*(value_right[1]-value_left[1])+\
                            value_central[order-1]*(value_right[2]-value_left[2]))
    
    refinement_quantity = refinement_quantity/delta_x
    
    return refinement_quantity


def compute_refinement_criterion_old(self,
                                values: np.ndarray,
                                orders: list,
                                max_order: int,
                                numbers_of_variables: list,
                                n: int,
                                boundary_interfaces: list,
                                delta_t: float,
                                delta_x: float,
                                tolerance_increase = 0.0015) -> np.ndarray:
    
    increase_criterion_flags = np.zeros(n,dtype=int)
    breakdown_estimators_increase = np.zeros(n)

    backward_differences = np.zeros(n)
    forward_differences = np.zeros(n)

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

        backward_differences[l-1] = np.abs(self.compute_refinement_estimator(values[l,:],
                                                                                    left_boundary_value,
                                                                                    values[l,:],
                                                                                    order,
                                                                                    max_order,
                                                                                    delta_t,
                                                                                    delta_x)) 
        forward_differences[l-1] = np.abs(self.compute_refinement_estimator(values[l,:],
                                                                                    values[l,:],
                                                                                    values[l+1,:],
                                                                                    order,
                                                                                    max_order,
                                                                                    delta_t,
                                                                                    delta_x))                  
                                                
        for i in range(l+1,r):
            backward_differences[i-1] = np.abs(self.compute_refinement_estimator(values[i,:],
                                                                                        values[i-1,:],
                                                                                        values[i,:],
                                                                                        order,
                                                                                        max_order,
                                                                                        delta_t,
                                                                                        delta_x)) 
            forward_differences[i-1] = np.abs(self.compute_refinement_estimator(values[i,:],
                                                                                        values[i,:],
                                                                                        values[i+1,:],
                                                                                        order,
                                                                                        max_order,
                                                                                        delta_t,
                                                                                        delta_x))        
        backward_differences[r-1] = np.abs(self.compute_refinement_estimator(values[r,:],
                                                                                    values[r-1,:],
                                                                                    values[r,:],
                                                                                    order,
                                                                                    max_order,
                                                                                    delta_t,
                                                                                    delta_x)) 
        forward_differences[r-1] = np.abs(self.compute_refinement_estimator(values[r,:],
                                                                                    values[r,:],
                                                                                    right_boundary_value,
                                                                                    order,
                                                                                    max_order,
                                                                                    delta_t,
                                                                                    delta_x)) 

    n_variables_prev = n_variables
    order = orders[-1]
    n_variables = numbers_of_variables[-1]

    n_left = min(n_variables_prev,n_variables)

    l = r+1  

    left_boundary_value = input_values[l,:]            
    left_boundary_value[:n_left] = input_values[l-1,:n_left]  

    backward_differences[l-1] = np.abs(self.compute_refinement_estimator(values[l,:],
                                                                                left_boundary_value,
                                                                                values[l,:],
                                                                                order,
                                                                                max_order,
                                                                                delta_t,
                                                                                delta_x)) 
    forward_differences[l-1] = np.abs(self.compute_refinement_estimator(values[l,:],
                                                                                values[l,:],
                                                                                values[l+1,:],
                                                                                order,
                                                                                max_order,
                                                                                delta_t,
                                                                                delta_x))                  
                                            
    for i in range(l+1,n+1):
        backward_differences[i-1] = np.abs(self.compute_refinement_estimator(values[i,:],
                                                                                    values[i-1,:],
                                                                                    values[i,:],
                                                                                    order,
                                                                                    max_order,
                                                                                    delta_t,
                                                                                    delta_x)) 
        forward_differences[i-1] = np.abs(self.compute_refinement_estimator(values[i,:],
                                                                                    values[i,:],
                                                                                    values[i+1,:],
                                                                                    order,
                                                                                    max_order,
                                                                                    delta_t,
                                                                                    delta_x)) 

    breakdown_estimators_increase = np.maximum(forward_differences,backward_differences)

    for i in range(n):
        if breakdown_estimators_increase[i] > tolerance_increase: 
            increase_criterion_flags[i] = 2            
    return breakdown_estimators_increase, increase_criterion_flags

def compute_refinement_criterion(self,
                                values: np.ndarray,
                                orders: list,
                                max_order: int,
                                numbers_of_variables: list,
                                n: int,
                                boundary_interfaces: list,
                                delta_t: float,
                                delta_x: float,
                                tols_refinement) -> np.ndarray:

    increase_criterion_flags = np.zeros(n,dtype=int)
    model_error_estimators_increase = np.zeros(n)

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

        model_error_estimators_increase[l-1], increase_criterion_flags[l-1] = self.compute_refinement_estimator(
                                                                                left_boundary_value,
                                                                                values[l,:],
                                                                                values[l+1,:],
                                                                                order,
                                                                                max_order,
                                                                                delta_t,
                                                                                delta_x,
                                                                                tols_refinement) 

        for i in range(l+1,r):
            model_error_estimators_increase[i-1], increase_criterion_flags[i-1] = self.compute_refinement_estimator(values[i-1,:],
                                                                                        values[i,:],
                                                                                        values[i+1,:],
                                                                                        order,
                                                                                        max_order,
                                                                                        delta_t,
                                                                                        delta_x,
                                                                                        tols_refinement)     
        model_error_estimators_increase[r-1], increase_criterion_flags[r-1] = self.compute_refinement_estimator(values[r-1,:],
                                                                                    values[r,:],
                                                                                    right_boundary_value,
                                                                                    order,
                                                                                    max_order,
                                                                                    delta_t,
                                                                                    delta_x,
                                                                                    tols_refinement) 

    n_variables_prev = n_variables
    order = orders[-1]
    n_variables = numbers_of_variables[-1]

    n_left = min(n_variables_prev,n_variables)

    l = r+1  

    left_boundary_value = input_values[l,:]            
    left_boundary_value[:n_left] = input_values[l-1,:n_left]  

    model_error_estimators_increase[l-1], increase_criterion_flags[l-1] = self.compute_refinement_estimator(left_boundary_value,
                                                                            values[l,:],
                                                                            values[l+1,:],
                                                                            order,
                                                                            max_order,
                                                                            delta_t,
                                                                            delta_x,
                                                                            tols_refinement)                 
                                            
    for i in range(l+1,n+1):
        model_error_estimators_increase[i-1], increase_criterion_flags[i-1] = self.compute_refinement_estimator(values[i-1,:],
                                                                    values[i,:],
                                                                    values[i+1,:],
                                                                    order,
                                                                    max_order,
                                                                    delta_t,
                                                                    delta_x,
                                                                    tols_refinement) 
        
    model_error_estimators_increase = model_error_estimators_increase/delta_x
        
    return model_error_estimators_increase, increase_criterion_flags

def compute_coarsening_criterion_old(self,
                                values: np.ndarray,
                                delta_t: float,
                                delta_x: float,
                                orders: list,
                                max_order: int,
                                numbers_of_variables: list,
                                n: int,
                                boundary_interfaces: list,
                                tols_decrease,
                                increase_criterion_flags: np.ndarray) -> np.ndarray:

    decrease_criterion_flags = np.zeros(n,dtype=int)
    breakdown_estimators_decrease = np.zeros(n)

    backward_differences = np.zeros(n)
    forward_differences = np.zeros(n)

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
            backward_differences[l-1] = np.abs(self.compute_coarsening_estimator(values[l,:],
                                                                                        left_boundary_value,
                                                                                        values[l,:],
                                                                                        order,
                                                                                        delta_t,
                                                                                        delta_x)) 
            forward_differences[l-1] = np.abs(self.compute_coarsening_estimator(values[l,:],
                                                                                        values[l,:],
                                                                                        values[l+1,:],
                                                                                        order,
                                                                                        delta_t,
                                                                                        delta_x))                                                        
            for i in range(l+1,r):
                backward_differences[i-1] = np.abs(self.compute_coarsening_estimator(values[i,:],
                                                                                            values[i-1,:],
                                                                                            values[i,:],
                                                                                            order,
                                                                                            delta_t,
                                                                                            delta_x)) 
                forward_differences[i-1] = np.abs(self.compute_coarsening_estimator(values[i,:],
                                                                                            values[i,:],
                                                                                            values[i+1,:],
                                                                                            order))           
            backward_differences[r-1] = np.abs(self.compute_coarsening_estimator(values[r,:],
                                                                                        values[r-1,:],
                                                                                        values[r,:],
                                                                                        order,
                                                                                        delta_t,
                                                                                        delta_x)) 
            forward_differences[r-1] = np.abs(self.compute_coarsening_estimator(values[r,:],
                                                                                        values[r,:],
                                                                                        right_boundary_value,
                                                                                        order,
                                                                                        delta_t,
                                                                                        delta_x)) 
    n_variables_prev = n_variables
    order = orders[-1]
    n_variables = numbers_of_variables[-1]

    n_left = min(n_variables_prev,n_variables)

    l = r+1  

    left_boundary_value = input_values[l,:]            
    left_boundary_value[:n_left] = input_values[l-1,:n_left]  

    if order > 3:
        backward_differences[l-1] = np.abs(self.compute_coarsening_estimator(values[l,:],
                                                                                    left_boundary_value,
                                                                                    values[l,:],
                                                                                    order,
                                                                                    delta_t,
                                                                                    delta_x)) 
        forward_differences[l-1] = np.abs(self.compute_coarsening_estimator(values[l,:],
                                                                                    values[l,:],
                                                                                    values[l+1,:],
                                                                                    order,
                                                                                    delta_t,
                                                                                    delta_x))                                                        
        for i in range(l+1,n+1):
            backward_differences[i-1] = np.abs(self.compute_coarsening_estimator(values[i,:],
                                                                                        values[i-1,:],
                                                                                        values[i,:],
                                                                                        order,
                                                                                        delta_t,
                                                                                        delta_x)) 
            forward_differences[i-1] = np.abs(self.compute_coarsening_estimator(values[i,:],
                                                                                        values[i,:],
                                                                                        values[i+1,:],
                                                                                        order,
                                                                                        delta_t,
                                                                                        delta_x)) 

    breakdown_estimators_decrease = np.maximum(backward_differences,forward_differences)

    for i in range(n):
        if increase_criterion_flags[i] == 0 and breakdown_estimators_decrease[i] < tols_decrease:
            decrease_criterion_flags[i] = -2

    return breakdown_estimators_decrease, decrease_criterion_flags

def compute_coarsening_criterion(self,
                                values: np.ndarray,
                                delta_t: float,
                                delta_x: float,
                                orders: list,
                                max_order: int,
                                numbers_of_variables: list,
                                n: int,
                                boundary_interfaces: list,
                                tols_coarsening,
                                increase_criterion_flags: np.ndarray) -> np.ndarray:

    decrease_criterion_flags = np.zeros(n,dtype=int)
    model_error_estimators_decrease = np.zeros(n)

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
                model_error_estimators_decrease[l-1], decrease_criterion_flags[l-1] = self.compute_coarsening_estimator(left_boundary_value,
                                                                                            values[l,:],
                                                                                            values[l+1,:],
                                                                                            order,
                                                                                            delta_t,
                                                                                            delta_x,
                                                                                            tols_coarsening)                                                       
            for i in range(l+1,r):
                if increase_criterion_flags[i-1] == 0:
                    model_error_estimators_decrease[i-1], decrease_criterion_flags[i-1] = self.compute_coarsening_estimator(values[i-1,:],
                                                                                                values[i,:],
                                                                                                values[i+1,:],
                                                                                                order,
                                                                                                delta_t,
                                                                                                delta_x,
                                                                                                tols_coarsening)          
            if increase_criterion_flags[r-1] == 0:
                model_error_estimators_decrease[r-1], decrease_criterion_flags[r-1] = self.compute_coarsening_estimator(values[r-1,:],
                                                                                            values[r,:],
                                                                                            right_boundary_value,
                                                                                            order,
                                                                                            delta_t,
                                                                                            delta_x,
                                                                                            tols_coarsening) 
    n_variables_prev = n_variables
    order = orders[-1]
    n_variables = numbers_of_variables[-1]

    n_left = min(n_variables_prev,n_variables)

    l = r+1  

    left_boundary_value = input_values[l,:]            
    left_boundary_value[:n_left] = input_values[l-1,:n_left]  

    if order > 3:
        if increase_criterion_flags[l-1] == 0:
            model_error_estimators_decrease[l-1], decrease_criterion_flags[l-1] = self.compute_coarsening_estimator(
                                                                                        left_boundary_value,
                                                                                        values[l,:],
                                                                                        values[l+1,:],
                                                                                        order,
                                                                                        delta_t,
                                                                                        delta_x,
                                                                                        tols_coarsening)                                                      
        for i in range(l+1,n+1):
            if increase_criterion_flags[i-1] == 0:
                model_error_estimators_decrease[i-1], decrease_criterion_flags[i-1] = self.compute_coarsening_estimator(values[i-1,:],
                                                                                            values[i,:],
                                                                                            values[i+1,:],
                                                                                            order,
                                                                                            delta_t,
                                                                                            delta_x,
                                                                                            tols_coarsening) 

    return model_error_estimators_decrease, decrease_criterion_flags