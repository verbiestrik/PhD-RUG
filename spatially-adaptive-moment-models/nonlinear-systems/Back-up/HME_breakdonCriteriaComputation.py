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
                                   tolerance_up_flow_gradient = 0.00001,
                                   tolerance_down_last_moment = 0.00001) -> tuple[np.ndarray,np.ndarray]:       
        """
        computes the breakdown criteria for adaptive simulation

        Parameters
        ----------
        values : list of numpy 1D arrays
            the values of the variables in each mesh cell
        delta_x : float
            grid cell size
        delta_t : float
            current time step
        max_order : integer
            max order of the moment model
        orders_cellwise : list of integers
            the order in each cell
        number_of_variables_cellwise : list of integers
            the number of variables in each cell
        dom_decomp_val_res1 : np.ndarray
            value of res1 in each grid cell
        dom_decomp_val_res2 : np.ndarray
            value of res2 in each grid cell
        tolerance_up_flow_gradient : float
            threshold for increase-criterion related to the flow gradients
        tolerance_down_last_moment : float
            threshold for decrease-criterion related to the magnitude of the last moment
        
        Returns
        -------
        breakdown_estimators : np.ndarray
            values of the breakdown estimators in each grid cell
        breakdown_criterion_flags : np.ndarray
            flags for increasing or reducing the order in each grid cell
            this array is filled with the values of the changes in order in each grid cell
        """
        breakdown_criterion_flags = np.zeros(n,dtype=int)
        breakdown_estimators = np.zeros((n,2))

        # # print(np.abs((values[2,0] - values[0,0]))/(2*delta_x*max(0.1,np.abs(values[1,0]))))
        # for i in range(n):
        #     loc_M = orders_cellwise[i+1]-1
        #     breakdown_estimators[i,0] = max(np.abs(values[i+1,loc_M+1]),np.abs(values[i+1,loc_M]))
        #     # breakdown_estimators[i,0] = np.sqrt(values[i+1,loc_M+1]**2+values[i+1,loc_M]**2)
        #     # for j in range(numbers_of_variables_cellwise[i+1]):
        #     #     # breakdown_estimators[i,1] += (max(np.abs((values[i+1,j] - values[i,j])),np.abs((values[i+2,j] - values[i+1,j])))/(delta_x))**2
        #     #     breakdown_estimators[i,1] += (np.abs((values[i+2,j] - values[i,j]))/(2*delta_x))**2 
        #     # j = 0
        #     # breakdown_estimators[i,1] = (np.abs((values[i+2,j] - values[i,j]))/(2*delta_x*max(1,np.abs(values[i+1,j]))))**2
        #     if loc_M == 2:
        #         breakdown_estimators[i,1] = np.abs(delta_t/(4*delta_x)*values[i+1,0]*values[i+1,2]*(values[i+2,2]-values[i,2]))
        #     elif loc_M == 3:
        #         # breakdown_estimators[i,1] = np.abs(delta_t/(2*delta_x)*\
        #         #                                    (\
        #         #                                        values[i+1,3]*values[i+1,2]*(values[i+2,0]-values[i,0])/values[i+1,0]\
        #         #                                         +values[i+1,3]*(values[i+2,2]-values[i,2]))\
        #         #                                         -values[i+1,2]*(values[i+2,3]-values[i,3]
        #         #                                                         )
        #         #                                         )
        #         breakdown_estimators[i,1] = max(
        #             np.abs(delta_t/(2*delta_x)*\
        #                             (\
        #                                 values[i+1,3]*values[i+1,2]*(values[i+2,0]-values[i,0])/values[i+1,0]\
        #                                 -3/2*values[i+1,3]*(values[i+2,2]-values[i,2]))\
        #                                 -values[i+1,2]*(values[i+2,3]-values[i,3]
        #                                                 )
        #                                 ),
        #             np.abs(delta_t/(2*delta_x)*(3*values[i+1,3]/values[i+1,0]*(values[i+2,3]-values[i,3])))
        #         )
        #     elif loc_M == 4:
        #         # breakdown_estimators[i,1] = np.abs(delta_t/(2*delta_x)*\
        #         #                                    (\
        #         #                                        values[i+1,4]*values[i+1,2]*(values[i+2,0]-values[i,0])/values[i+1,0]\
        #         #                                         +values[i+1,4]*(values[i+2,2]-values[i,2]))\
        #         #                                         +3*values[i+1,3]/values[i+1,0]*(values[i+2,3]-values[i,3])\
        #         #                                         -values[i+1,2]*(values[i+2,4]-values[i,4]
        #         #                                                         )
        #         #                                         )
        #         breakdown_estimators[i,1] = max(
        #                                         np.abs(delta_t/(2*delta_x)*\
        #                                            (\
        #                                                values[i+1,4]*values[i+1,2]*(values[i+2,0]-values[i,0])/values[i+1,0]\
        #                                                 -2*values[i+1,4]*(values[i+2,2]-values[i,2]))\
        #                                                 -values[i+1,2]*(values[i+2,4]-values[i,4]\
        #                                                 +3*values[i+1,3]/values[i+1,0]*(values[i+2,3]-values[i,3])
        #                                                                 )
        #                                                 ),
        #                                         np.abs(delta_t/(2*delta_x)*(1/2*values[i+1,2]*values[i+1,3])*(values[i+2,2]-values[i,2])\
        #                                                -3*values[i+1,4]/values[i+1,0]*(values[i+2,3]-values[i,3]))
        #         ) 
        #         breakdown_estimators[i,0] = max(
        #             delta_t/(2*delta_x)*np.abs(4*values[i+1,3]*(values[i+2,1]-values[i,1])\
        #                                         +values[i+1,0]*values[i+1,2]/2*(values[i+2,2]-values[i,2])\
        #                                         +values[i+1,1]*(values[i+2,3]-values[i,3])\
        #                                         +4*(values[i+2,4]-values[i,4])),
        #             delta_t/(2*delta_x)*np.abs(-values[i+1,2]*values[i+1,3]/values[i+1,0]*(values[i+2,0]-values[i,0])\
        #                                         -values[i+1,3]*(values[i+2,2]-values[i,2])\
        #                                         +values[i+1,2]*(values[i+2,3]-values[i,3])\
        #                                         +values[i+1,1]*(values[i+2,4]-values[i,4])),
        #             delta_t/(2*delta_x)*np.abs(6/values[i+1,0])
        #         )                      
        #     else:
        #         # breakdown_estimators[i,1] = np.abs(delta_t/(2*delta_x)*\
        #         #                                    (\
        #         #                                        values[i+1,loc_M]*values[i+1,2]*(values[i+2,0]-values[i,0])/values[i+1,0]\
        #         #                                         -(values[i+1,2]*values[i+1,loc_M-2]/2-values[i+1,loc_M])*(values[i+2,2]-values[i,2]))\
        #         #                                         +3*values[i+1,3]/values[i+1,0]*(values[i+2,3]-values[i,3])\
        #         #                                         -values[i+1,2]*(values[i+2,loc_M]-values[i,loc_M]
        #         #                                                         )
        #         breakdown_estimators[i,1] = max(
        #                                         np.abs(delta_t/(2*delta_x)*\
        #                                            (\
        #                                                values[i+1,loc_M]*values[i+1,2]*(values[i+2,0]-values[i,0])/values[i+1,0]\
        #                                                 -(values[i+1,2]*values[i+1,loc_M-2]/2-values[i+1,loc_M])*(values[i+2,2]-values[i,2]))\
        #                                                 +3*values[i+1,3]/values[i+1,0]*(values[i+2,3]-values[i,3])\
        #                                                 -values[i+1,2]*(values[i+2,loc_M]-values[i,loc_M]
        #                                                                 )
        #                                                 ),
        #                                         np.abs(delta_t/(2*delta_x)*(1/2*values[i+1,2]*values[i+1,loc_M-1])*(values[i+2,2]-values[i,2])\
        #                                                -3*values[i+1,loc_M]/values[i+1,0]*(values[i+2,3]-values[i,3]))
        #         )
        #     if loc_M == 5:
        #         breakdown_estimators[i,0] = max(
        #             delta_t/(2*delta_x)*np.abs(-values[i+1,2]*values[i+1,3]/values[i+1,0]*(values[i+2,0]-values[i,0])\
        #                                         +5*values[i+1,4]*(values[i+2,1]-values[i,1])\
        #                                         +3*values[i+1,3]/2*(values[i+2,2]-values[i,2])\
        #                                         +values[i+1,2]*(values[i+2,3]-values[i,3])\
        #                                         +values[i+1,1]*(values[i+2,4]-values[i,4])\
        #                                         +5*(values[i+2,5]-values[i,5])),
        #             delta_t/(2*delta_x)*np.abs(-values[i+1,2]*values[i+1,4]/values[i+1,0]*(values[i+2,0]-values[i,0])\
        #                                         -values[i+1,4]*(values[i+2,2]-values[i,2])\
        #                                         -3*values[i+1,3]/values[i+1,0]*(values[i+2,3]-values[i,3])\
        #                                         +values[i+1,2]*(values[i+2,4]-values[i,4])\
        #                                         +values[i+1,1]*(values[i+2,5]-values[i,5])),
        #             delta_t/(2*delta_x)*np.abs(4*(values[i+2,4]-values[i,4]))
        #         )   
        #     elif loc_M == 6:
        #         breakdown_estimators[i,0] = max(
        #             delta_t/(2*delta_x)*np.abs(-values[i+1,2]*values[i+1,4]/values[i+1,0]*(values[i+2,0]-values[i,0])\
        #                                         +6*values[i+1,5]*(values[i+2,1]-values[i,1])\
        #                                         +2*values[i+1,4]*(values[i+2,2]-values[i,2])\
        #                                         -3*values[i+1,3]/values[i+1,0]*(values[i+2,3]-values[i,3])\
        #                                         +values[i+1,2]*(values[i+2,4]-values[i,4])\
        #                                         +values[i+1,1]*(values[i+2,5]-values[i,5])\
        #                                         +6*(values[i+2,6]-values[i,6])),
        #             delta_t/(2*delta_x)*np.abs(-values[i+1,2]*values[i+1,5]/values[i+1,0]*(values[i+2,0]-values[i,0])\
        #                                         +(values[i+1,2]*values[i+1,3]/2-values[i+1,5])*(values[i+2,2]-values[i,2])\
        #                                         -3*values[i+1,4]/values[i+1,0]*(values[i+2,3]-values[i,3])\
        #                                         +values[i+1,2]*(values[i+2,5]-values[i,5])\
        #                                         +values[i+1,1]*(values[i+2,6]-values[i,6])),
        #             delta_t/(2*delta_x)*np.abs(5*(values[i+2,5]-values[i,5]))
        #         )  
        #     elif loc_M > 6:
        #         breakdown_estimators[i,0] = max(
        #             delta_t/(2*delta_x)*np.abs(-values[i+1,2]*values[i+1,loc_M-2]/values[i+1,0]*(values[i+2,0]-values[i,0])\
        #                                         +loc_M*values[i+1,loc_M-1]*(values[i+2,1]-values[i,1])\
        #                                         +(values[i+1,2]*values[i+1,loc_M-4]/2+(loc_M-2)/2*values[i+1,loc_M-2])*(values[i+2,2]-values[i,2])\
        #                                         -3*values[i+1,loc_M-3]/values[i+1,0]*(values[i+2,3]-values[i,3])\
        #                                         +values[i+1,2]*(values[i+2,loc_M-2]-values[i,loc_M-2])\
        #                                         +values[i+1,1]*(values[i+2,loc_M-1]-values[i,loc_M-1])\
        #                                         +loc_M*(values[i+2,loc_M]-values[i,loc_M])),
        #             delta_t/(2*delta_x)*np.abs(-values[i+1,2]*values[i+1,loc_M-1]/values[i+1,0]*(values[i+2,0]-values[i,0])\
        #                                         +(values[i+1,2]*values[i+1,loc_M-3]/2-values[i+1,loc_M-1])*(values[i+2,2]-values[i,2])\
        #                                         -3*values[i+1,loc_M-2]/values[i+1,0]*(values[i+2,3]-values[i,3])\
        #                                         +values[i+1,2]*(values[i+2,loc_M-1]-values[i,loc_M-1])\
        #                                         +values[i+1,1]*(values[i+2,loc_M]-values[i,loc_M])),
        #             delta_t/(2*delta_x)*np.abs((loc_M-1)*(values[i+2,loc_M-1]-values[i,loc_M-1]))
        #         )                   
        # breakdown_estimators[:,1] = breakdown_estimators[:,1]

        for i in range(n):
            loc_M = orders_cellwise[i+1]
            if loc_M > 3:
                if loc_M == 4:
                    breakdown_estimators[i,0] = np.abs((loc_M-1)*(values[i+2,loc_M-2]-values[i,loc_M-2]))                 
                elif loc_M == 5:
                    breakdown_estimators[i,0] = np.abs(
                        (loc_M-1)*(values[i+2,loc_M-2]-values[i,loc_M-2])\
                        +(loc_M-1)*values[i+1,loc_M-2]*(values[i+2,1]-values[i,1])
                    )
                else:
                    breakdown_estimators[i,0] = np.abs(
                        (loc_M-1)*(values[i+2,loc_M-2]-values[i,loc_M-2])\
                            +(loc_M-1)*values[i+1,loc_M-2]*(values[i+2,1]-values[i,1])\
                                (loc_M-1)/2*values[i+1,loc_M-3]*(values[i+2,2]-values[i,2])
                    )
            if loc_M < max_order - 1:
                pass

        for i in range(n):
            if orders_cellwise[i+1] < max_order-1 and breakdown_estimators[i,1] > tolerance_up_flow_gradient:
                breakdown_criterion_flags[i] = 2
            else:
                if orders_cellwise[i+1] > 3 and breakdown_estimators[i,0] < tolerance_down_last_moment: 
                    breakdown_criterion_flags[i] = -2
        
        return breakdown_estimators, breakdown_criterion_flags

    def compute_breakdown_criteria_decrease_old(self,
                                   values: np.ndarray,
                                   padded_vectors_left,
                                   padded_vectors_right,
                                   n: int,
                                   delta_x: float,
                                   delta_t: float,
                                   max_order: int,
                                   orders_cellwise: list,
                                   orders: list,
                                   numbers_of_variables :list,
                                   boundary_interfaces,
                                   increase_criterion_flags: np.ndarray,
                                   tolerance_decrease = 0.0005) -> tuple[np.ndarray,np.ndarray]:       
        """
        TODO
        """
        decrease_criterion_flags = np.zeros(n,dtype=int)
        breakdown_estimators_decrease = np.full(n,tolerance_decrease+1)

        # for i in range(n):
        #     if increase_criterion_flags[i] == 0 and orders_cellwise[i+1] > 3:
        #         loc_M = orders_cellwise[i+1]
        #         if loc_M == 4:
        #             breakdown_estimators_decrease[i] = np.abs((loc_M-1)*(values[i+2,loc_M-2]-values[i,loc_M-2]))                 
        #         elif loc_M == 5:
        #             breakdown_estimators_decrease[i] = np.abs(
        #                 (loc_M-1)*(values[i+2,loc_M-2]-values[i,loc_M-2])\
        #                 +(loc_M-1)*values[i+1,loc_M-2]*(values[i+2,1]-values[i,1])
        #             )
        #         else:
        #             breakdown_estimators_decrease[i] = np.abs(
        #                 (loc_M-1)*(values[i+2,loc_M-2]-values[i,loc_M-2])\
        #                     +(loc_M-1)*values[i+1,loc_M-2]*(values[i+2,1]-values[i,1])\
        #                         (loc_M-1)/2*values[i+1,loc_M-3]*(values[i+2,2]-values[i,2])
        #             )
        #         if breakdown_estimators_decrease[i,0] < tolerance_decrease: 
        #             decrease_criterion_flags[i] = -2

        backward_differences = np.zeros(n)
        forward_differences = np.zeros(n)

        right_boundary_subdomain = 0

        # order = orders[0]
        # n_variables = n_variables[0]

        # left_boundary_subdomain = right_boundary_subdomain+1
        # right_boundary_subdomain = boundary_interfaces[0]                 

        # if order > 3:
        #     for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
        #         if increase_criterion_flags[i-1] == 0:
        #             if order == 4:
        #                 backward_differences[i-1] = np.abs((order-1)*(values[i,order-2]-values[i-1,order-2])) 
        #                 forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-2]-values[i,order-2]))                    
        #             elif order == 5:
        #                 backward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                     +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])
        #                 )
        #                 forward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i+1,order-2]-values[i,order-2])\
        #                     +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])
        #                 )
        #             else:
        #                 backward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                             (order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2])
        #                 )
        #                 forward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i+1,order-2]-values[i,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
        #                             (order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2])
        #                 )
        #             if breakdown_estimators_decrease[i-1,0] < tolerance_decrease: 
        #                 decrease_criterion_flags[i-1] = -2
        #     if order > orders[m+1]:
        #         i = right_boundary_subdomain
        #         if increase_criterion_flags[i] == 0:
        #             if order == 4:
        #                 backward_differences[i-1] = np.abs((order-1)*(values[i,order-2]-values[i-1,order-2])) 
        #                 forward_differences[i-1] = np.abs((order-1)*(padded_vectors_right[m,order-2]-values[i,order-2]))                    
        #             elif order == 5:
        #                 backward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                     +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])
        #                 )
        #                 forward_differences[i-1] = np.abs(
        #                     (order-1)*(padded_vectors_right[m,order-2]-values[i,order-2])\
        #                     +(order-1)*values[i,order-2]*(padded_vectors_right[m,1]-values[i,1])
        #                 )
        #             else:
        #                 backward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                             (order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2])
        #                 )
        #                 forward_differences[i-1] = np.abs(
        #                     (order-1)*(padded_vectors_right[i+1,order-2]-values[i,order-2])\
        #                         +(order-1)*values[i,order-2]*(padded_vectors_right[m,1]-values[i,1])\
        #                             (order-1)/2*values[i,order-3]*(padded_vectors_right[m,2]-values[i,2])
        #                 )
        #             if breakdown_estimators_decrease[i-1,0] < tolerance_decrease: 
        #                 decrease_criterion_flags[i-1] = -2  
        #             else:
        #                 decrease_criterion_flags[i-1] = 0

        # for m in range(1,len(boundary_interfaces)):
        #     order = orders[m]
        #     n_variables = n_variables[m]

        #     left_boundary_subdomain = right_boundary_subdomain+1
        #     right_boundary_subdomain = boundary_interfaces[m]                 

        #     if order > 3:
        #         if order > orders[m-1]:
        #             values[left_boundary_subdomain,:] = padded_vectors_left[m-1,:]  
        #         for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
        #             if increase_criterion_flags[i-1] == 0:
        #                 if order == 4:
        #                     backward_differences[i-1] = np.abs((order-1)*(values[i,order-2]-values[i-1,order-2])) 
        #                     forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-2]-values[i,order-2]))                    
        #                 elif order == 5:
        #                     backward_differences[i-1] = np.abs(
        #                         (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])
        #                     )
        #                     forward_differences[i-1] = np.abs(
        #                         (order-1)*(values[i+1,order-2]-values[i,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])
        #                     )
        #                 else:
        #                     backward_differences[i-1] = np.abs(
        #                         (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                             +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                                 (order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2])
        #                     )
        #                     forward_differences[i-1] = np.abs(
        #                         (order-1)*(values[i+1,order-2]-values[i,order-2])\
        #                             +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
        #                                 (order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2])
        #                     )
        #                 if breakdown_estimators_decrease[i-1,0] < tolerance_decrease: 
        #                     decrease_criterion_flags[i-1] = -2
        #                 else:
        #                     decrease_criterion_flags[i-1] = 0
        #         if order > orders[m+1]:
        #             i = right_boundary_subdomain
        #             if increase_criterion_flags[i] == 0:
        #                 if order == 4:
        #                     backward_differences[i-1] = np.abs((order-1)*(values[i,order-2]-values[i-1,order-2])) 
        #                     forward_differences[i-1] = np.abs((order-1)*(padded_vectors_right[m,order-2]-values[i,order-2]))                    
        #                 elif order == 5:
        #                     backward_differences[i-1] = np.abs(
        #                         (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])
        #                     )
        #                     forward_differences[i-1] = np.abs(
        #                         (order-1)*(padded_vectors_right[m,order-2]-values[i,order-2])\
        #                         +(order-1)*values[i,order-2]*(padded_vectors_right[m,1]-values[i,1])
        #                     )
        #                 else:
        #                     backward_differences[i-1] = np.abs(
        #                         (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                             +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                                 (order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2])
        #                     )
        #                     forward_differences[i-1] = np.abs(
        #                         (order-1)*(padded_vectors_right[i+1,order-2]-values[i,order-2])\
        #                             +(order-1)*values[i,order-2]*(padded_vectors_right[m,1]-values[i,1])\
        #                                 (order-1)/2*values[i,order-3]*(padded_vectors_right[m,2]-values[i,2])
        #                     )
        #                 if breakdown_estimators_decrease[i-1,0] < tolerance_decrease: 
        #                     decrease_criterion_flags[i-1] = -2                    

        # order = orders[-1]
        # n_variables = n_variables[-1]

        # left_boundary_subdomain = right_boundary_subdomain+1
        # right_boundary_subdomain = n                 

        # if order > 3:
        #     if order > orders[m-1]:
        #         values[left_boundary_subdomain,:] = padded_vectors_left[m-1,:]  
        #     for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
        #         if increase_criterion_flags[i-1] == 0:
        #             if order == 4:
        #                 backward_differences[i-1] = np.abs((order-1)*(values[i,order-2]-values[i-1,order-2])) 
        #                 forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-2]-values[i,order-2]))                    
        #             elif order == 5:
        #                 backward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                     +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])
        #                 )
        #                 forward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i+1,order-2]-values[i,order-2])\
        #                     +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])
        #                 )
        #             else:
        #                 backward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                             (order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2])
        #                 )
        #                 forward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i+1,order-2]-values[i,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
        #                             (order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2])
        #                 )
        #             if breakdown_estimators_decrease[i-1,0] < tolerance_decrease: 
        #                 decrease_criterion_flags[i-1] = -2

        for m in range(len(boundary_interfaces)):
            order = orders[m]
            n_variables = numbers_of_variables[m]

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = boundary_interfaces[m]                 

            if order > 3:
                if order > orders[m+1]:
                    values[right_boundary_subdomain+1,:] = padded_vectors_right[m,:]  
                else:
                    values[right_boundary_subdomain,:] = padded_vectors_left[m,:]  
                for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                    if increase_criterion_flags[i-1] == 0:
                        if order == 4:
                            backward_differences[i-1] = np.abs((order-1)*(values[i,order-2]-values[i-1,order-2])) 
                            forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-2]-values[i,order-2]))                    
                        elif order == 5:
                            backward_differences[i-1] = np.abs(
                                (order-1)*(values[i,order-2]-values[i-1,order-2])\
                                +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])
                            )
                            forward_differences[i-1] = np.abs(
                                (order-1)*(values[i+1,order-2]-values[i,order-2])\
                                +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])
                            )
                        else:
                            backward_differences[i-1] = np.abs(
                                (order-1)*(values[i,order-2]-values[i-1,order-2])\
                                    +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
                                        +(order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2])
                            )
                            forward_differences[i-1] = np.abs(
                                (order-1)*(values[i+1,order-2]-values[i,order-2])\
                                    +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
                                        +(order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2])
                            )             

        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        left_boundary_subdomain = right_boundary_subdomain+1
        right_boundary_subdomain = n                 

        if order > 3:
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                if increase_criterion_flags[i-1] == 0:
                    if order == 4:
                        backward_differences[i-1] = np.abs((order-1)*(values[i,order-2]-values[i-1,order-2])) 
                        forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-2]-values[i,order-2]))                    
                    elif order == 5:
                        backward_differences[i-1] = np.abs(
                            (order-1)*(values[i,order-2]-values[i-1,order-2])\
                            +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])
                        )
                        forward_differences[i-1] = np.abs(
                            (order-1)*(values[i+1,order-2]-values[i,order-2])\
                            +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])
                        )
                    else:
                        backward_differences[i-1] = np.abs(
                            (order-1)*(values[i,order-2]-values[i-1,order-2])\
                                +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
                                    +(order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2])
                        )
                        forward_differences[i-1] = np.abs(
                            (order-1)*(values[i+1,order-2]-values[i,order-2])\
                                +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
                                    +(order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2])
                        )

        # breakdown_estimators_decrease = (backward_differences+forward_differences)/delta_x
        breakdown_estimators_decrease = np.maximum(backward_differences,forward_differences)/delta_x

        right_boundary_subdomain = 0
        for m in range(1,len(boundary_interfaces)):
            left_boundary_subdomain = right_boundary_subdomain + 1
            right_boundary_subdomain = boundary_interfaces[m]
            for i in range(left_boundary_subdomain,right_boundary_subdomain):
                breakdown_estimators_decrease[i-1] = np.maximum(backward_differences[i-1],forward_differences[i-1])/delta_x
            breakdown_estimators_decrease[right_boundary_subdomain-1] = backward_differences[right_boundary_subdomain-1]/delta_x 

        for i in range(n):
            if breakdown_estimators_decrease[i] < tolerance_decrease: 
                decrease_criterion_flags[i] = -2

        return breakdown_estimators_decrease, decrease_criterion_flags


    def compute_breakdown_criteria_increase_old(self,
                                   values: np.ndarray,
                                   padded_vectors_left,
                                   padded_vectors_right,
                                   n: int,
                                   delta_x: float,
                                   delta_t: float,
                                   max_order: int,
                                   orders_cellwise: list,
                                   orders,
                                   numbers_of_variables,
                                   boundary_interfaces: list,
                                   tolerance_increase = 0.005) -> tuple[np.ndarray,np.ndarray]:       
        """
        TODO
        """
        increase_criterion_flags = np.zeros(n,dtype=int)
        breakdown_estimators_increase = np.zeros(n)

        backward_differences = np.zeros(n)
        forward_differences = np.zeros(n)

        right_boundary_subdomain = 0

        for m in range(len(boundary_interfaces)):
            order = orders[m]
            n_variables = numbers_of_variables[m]

            left_boundary_subdomain = right_boundary_subdomain+1
            right_boundary_subdomain = boundary_interfaces[m]                 

            # if order > orders[m+1]:
            #     values[right_boundary_subdomain+1,:] = padded_vectors_right[m,:]  
            # else:
            #     values[right_boundary_subdomain,:] = padded_vectors_left[m,:]  
            for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
                if order == 2:
                    backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
                    forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                    
                elif order == 3:
                    backward_differences[i-1] = np.abs(
                        4*values[i,3]*(values[i,1]-values[i-1,1])+4*(values[i,4]*values[i-1,4])
                    )
                    forward_differences[i-1] = np.abs(
                        4*values[i,3]*(values[i+1,1]-values[i,1])+4*(values[i+1,4]*values[i,4])
                    )
                elif order == 4:
                    backward_differences[i-1] = np.abs(
                        5*values[i,4]*(values[i,1]-values[i-1,1])+5/2*values[i,3]*(values[i,2]-values[i-1,2])\
                            +5*(values[i,5]*values[i-1,5])
                    )
                    forward_differences[i-1] = np.abs(
                        5*values[i,4]*(values[i+1,1]-values[i,1])+5/2*values[i,3]*(values[i+1,2]-values[i,2])\
                            +5*(values[i+1,5]*values[i,5])
                    )
                elif order == max_order:
                    backward_differences[i-1] = np.abs(
                        (order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
                            +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2])
                    )
                    forward_differences[i-1] = np.abs(
                        (order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
                            +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])
                    )     
                else:
                    backward_differences[i-1] = np.abs(
                        (order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
                            +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2])\
                            +(order+1)*(values[i,order+1]*values[i-1,order+1])
                    )
                    forward_differences[i-1] = np.abs(
                        (order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
                            +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])\
                            +(order+1)*(values[i+1,order+1]*values[i,order+1])
                    )           

        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        left_boundary_subdomain = right_boundary_subdomain+1
        right_boundary_subdomain = n                 

        for i in range(left_boundary_subdomain,right_boundary_subdomain+1):
            if order == 2:
                backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
                forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                    
            elif order == 3:
                backward_differences[i-1] = np.abs(
                    4*values[i,3]*(values[i,1]-values[i-1,1])+4*(values[i,4]*values[i-1,4])
                )
                forward_differences[i-1] = np.abs(
                    4*values[i,3]*(values[i+1,1]-values[i,1])+4*(values[i+1,4]*values[i,4])
                )
            elif order == 4:
                backward_differences[i-1] = np.abs(
                    5*values[i,4]*(values[i,1]-values[i-1,1])+5/2*values[i,3]*(values[i,2]-values[i-1,2])\
                        +5*(values[i,5]*values[i-1,5])
                )
                forward_differences[i-1] = np.abs(
                    5*values[i,4]*(values[i+1,1]-values[i,1])+5/2*values[i,3]*(values[i+1,2]-values[i,2])\
                        +5*(values[i+1,5]*values[i,5])
                )
            elif order == max_order:
                backward_differences[i-1] = np.abs(
                    (order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
                        +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2])
                )
                forward_differences[i-1] = np.abs(
                    (order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
                        +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])
                )  
            else:
                backward_differences[i-1] = np.abs(
                    (order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
                        +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2])\
                        +(order+1)*(values[i,order+1]*values[i-1,order+1])
                )
                forward_differences[i-1] = np.abs(
                    (order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
                        +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])\
                        +(order+1)*(values[i+1,order+1]*values[i,order+1])
                )   

        # breakdown_estimators_increase = (forward_differences+backward_differences)/delta_x
        breakdown_estimators_increase = np.maximum(forward_differences,backward_differences)/delta_x
        for i in range(n):
            if breakdown_estimators_increase[i] > tolerance_increase: 
                increase_criterion_flags[i] = 2            

        return breakdown_estimators_increase, increase_criterion_flags

    def compute_breakdown_criteria_decrease(self,
                                   values: np.ndarray,
                                #    padded_vectors_left,
                                #    padded_vectors_right,
                                   n: int,
                                   delta_x: float,
                                   delta_t: float,
                                   max_order: int,
                                   orders_cellwise: list,
                                   orders: list,
                                   numbers_of_variables :list,
                                   boundary_interfaces,
                                   increase_criterion_flags: np.ndarray,
                                   tolerance_decrease = 0.0005) -> tuple[np.ndarray,np.ndarray]:       
        """
        TODO
        """
        decrease_criterion_flags = np.zeros(n,dtype=int)
        breakdown_estimators_decrease = np.full(n,tolerance_decrease+1)

        backward_differences = np.zeros(n)
        forward_differences = np.zeros(n)

        input_values = np.copy(values)

        # r = 0

        # order = orders[0]
        # for m in range(len(boundary_interfaces)):
        #     prev_order = order
        #     order = orders[m]
        #     next_order = orders[m+1]
        #     n_variables = numbers_of_variables[m]

        #     l = r+1
        #     r = boundary_interfaces[m]                 
        #     if order > 3:
        #         if prev_order < order:
        #             if increase_criterion_flags[l-1] == 0:
        #                 if order == 4:
        #                     forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                    
        #                 elif order == 5:
        #                     backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-values[l-1,1]))
        #                     forward_differences[l-1] = np.abs(4*(values[l+1,4]-values[l,4])+4*values[l,3]*(values[l+1,1]-values[l,1]))
        #                 else:
        #                     backward_differences[l-1] = np.abs((order-1)*values[l,order-2]*(values[l,1]-values[l-1,1])\
        #                                 +(order-1)/2*values[l,order-3]*(values[l,2]-values[l-1,2]))
        #                     forward_differences[l-1] = np.abs((order-1)*(values[l+1,order-1]-values[l,order-1])\
        #                             +(order-1)*values[l,order-2]*(values[l+1,1]-values[l,1])\
        #                                 +(order-1)/2*values[l,order-3]*(values[l+1,2]-values[l,2]))     
        #             l = l+1
        #         # else:
        #         #     l = l-1
        #         if next_order < order:
        #             r2 = r
        #             if increase_criterion_flags[r-1] == 0:
        #                 if order == 4:
        #                     backward_differences[r-1] = np.abs(6/values[r,0]*(values[r,3]-values[r-1,3]))                
        #                 elif order == 5:
        #                     backward_differences[r-1] = np.abs(4*(values[r,4]-values[r-1,4])+4*values[r,3]*(values[r,1]-values[r-1,1]))
        #                     forward_differences[r-1] = np.abs(4*values[r,3]*(values[r+1,1]-values[r,1]))
        #                 else:
        #                     backward_differences[r-1] = np.abs((order-1)*(values[r,order-1]-values[r-1,order-1])\
        #                             +(order-1)*values[r,order-2]*(values[r,1]-values[r-1,1])\
        #                                 +(order-1)/2*values[r,order-3]*(values[r,2]-values[r-1,2]))
        #                     forward_differences[r-1] = np.abs((order-1)*values[r,order-2]*(values[r+1,1]-values[r,1])\
        #                                 +(order-1)/2*values[r,order-3]*(values[r+1,2]-values[r,2])) 
        #         else:
        #             r2 = r+1
        #         for i in range(l,r2):
        #             if increase_criterion_flags[i-1] == 0:
        #                 if order == 4:
        #                     backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
        #                     forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                    
        #                 elif order == 5:
        #                     backward_differences[i-1] = np.abs(4*(values[i,4]-values[i-1,4])\
        #                         +4*values[i,3]*(values[i,1]-values[i-1,1]))
        #                     forward_differences[i-1] = np.abs(4*(values[i+1,4]-values[i,4])\
        #                         +4*values[i,3]*(values[i+1,1]-values[i,1]))
        #                 else:
        #                     backward_differences[i-1] = np.abs((order-1)*(values[i,order-1]-values[i-1,order-1])\
        #                             +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                                 +(order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2]))
        #                     forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-1]-values[i,order-1])\
        #                             +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
        #                                 +(order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2]))             

        # order = orders[-1]
        # n_variables = numbers_of_variables[-1]

        # l = r+1
        # r = n                 

        # if order > 3:
        #     if orders[-2] < orders[-1]:
        #         if increase_criterion_flags[l-1] == 0:
        #             if order == 4:
        #                 forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                    
        #             elif order == 5:
        #                 backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-values[l-1,1]))
        #                 forward_differences[l-1] = np.abs(4*(values[l+1,4]-values[l,4])+4*values[l,3]*(values[l+1,1]-values[l,1]))
        #             else:
        #                 backward_differences[l-1] = np.abs((order-1)*values[l,order-2]*(values[l,1]-values[l-1,1])\
        #                             +(order-1)/2*values[l,order-3]*(values[l,2]-values[l-1,2]))
        #                 forward_differences[l-1] = np.abs((order-1)*(values[l+1,order-1]-values[l,order-1])\
        #                         +(order-1)*values[l,order-2]*(values[l+1,1]-values[l,1])\
        #                             +(order-1)/2*values[l,order-3]*(values[l+1,2]-values[l,2]))     
        #         l = l+1
        #     # else:
        #     #     l = l-1
        #     for i in range(l,r+1):
        #         if increase_criterion_flags[i-1] == 0:
        #             if order == 4:
        #                 backward_differences[i-1] = np.abs((order-1)*(values[i,order-2]-values[i-1,order-2])) 
        #                 forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-2]-values[i,order-2]))                    
        #             elif order == 5:
        #                 backward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                     +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])
        #                 )
        #                 forward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i+1,order-2]-values[i,order-2])\
        #                     +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])
        #                 )
        #             else:
        #                 backward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i,order-2]-values[i-1,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                             +(order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2])
        #                 )
        #                 forward_differences[i-1] = np.abs(
        #                     (order-1)*(values[i+1,order-2]-values[i,order-2])\
        #                         +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
        #                             +(order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2])
        #                 )

        # r = 0

        # order = orders[0]
        # n_variables = numbers_of_variables[0]
        # for m in range(len(boundary_interfaces)):

        #     prev_order = order
        #     prev_n_variables = n_variables
        #     order = orders[m]
        #     n_variables = numbers_of_variables[m]
        #     next_order = orders[m+1]
        #     next_n_variables = numbers_of_variables[m+1]

        #     l = r+1
        #     r = boundary_interfaces[m]  

        #     left_boundary_value = input_values[l-1,:]
        #     right_boundary_value = input_values[r+1,:]
        #     # left_boundary_value[prev_n_variables:] = values[l,prev_n_variables:]
        #     # right_boundary_value[prev_n_variables:] = values[r,prev_n_variables:]
        #     left_boundary_value[n_variables:] = input_values[l,n_variables:]
        #     right_boundary_value[n_variables:] = input_values[r,n_variables:]

        #     if order > 3:
        #         if increase_criterion_flags[l-1] == 0:
        #             if order == 4:
        #                 backward_differences[l-1] = np.abs(6/values[l,0]*(values[l,3]-left_boundary_value[3])) 
        #                 forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                    
        #             elif order == 5:
        #                 backward_differences[l-1] = np.abs(4*(values[l,4]-left_boundary_value[4])\
        #                     +4*values[l,3]*(values[l,1]-left_boundary_value[1]))
        #                 forward_differences[l-1] = np.abs(4*(values[l+1,4]-values[l,4])\
        #                     +4*values[l,3]*(values[l+1,1]-values[l,1]))
        #             else:
        #                 backward_differences[l-1] = np.abs((order-1)*(values[l,order-1]-left_boundary_value[order-1])\
        #                         +(order-1)*values[l,order-2]*(values[l,1]-left_boundary_value[1])\
        #                             +(order-1)/2*values[l,order-3]*(values[l,2]-left_boundary_value[2]))
        #                 forward_differences[l-1] = np.abs((order-1)*(values[l+1,order-1]-values[l,order-1])\
        #                         +(order-1)*values[l,order-2]*(values[l+1,1]-values[l,1])\
        #                             +(order-1)/2*values[l,order-3]*(values[l+1,2]-values[l,2]))                                     
        #         for i in range(l+1,r):
        #             if increase_criterion_flags[i-1] == 0:
        #                 if order == 4:
        #                     backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
        #                     forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                    
        #                 elif order == 5:
        #                     backward_differences[i-1] = np.abs(4*(values[i,4]-values[i-1,4])\
        #                         +4*values[i,3]*(values[i,1]-values[i-1,1]))
        #                     forward_differences[i-1] = np.abs(4*(values[i+1,4]-values[i,4])\
        #                         +4*values[i,3]*(values[i+1,1]-values[i,1]))
        #                 else:
        #                     backward_differences[i-1] = np.abs((order-1)*(values[i,order-1]-values[i-1,order-1])\
        #                             +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                                 +(order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2]))
        #                     forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-1]-values[i,order-1])\
        #                             +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
        #                                 +(order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2]))             
        #         if increase_criterion_flags[r-1] == 0:
        #             if order == 4:
        #                 backward_differences[r-1] = np.abs(6/values[r,0]*(values[r,3]-values[r-1,3])) 
        #                 forward_differences[r-1] = np.abs(6/values[r,0]*(right_boundary_value[3]-values[r,3]))                    
        #             elif order == 5:
        #                 backward_differences[r-1] = np.abs(4*(values[r,4]-values[r-1,4])\
        #                     +4*values[r,3]*(values[r,1]-values[r-1,1]))
        #                 forward_differences[r-1] = np.abs(4*(right_boundary_value[4]-values[r,4])\
        #                     +4*values[r,3]*(right_boundary_value[1]-values[r,1]))
        #             else:
        #                 backward_differences[r-1] = np.abs((order-1)*(values[r,order-1]-values[r-1,order-1])\
        #                         +(order-1)*values[r,order-2]*(values[r,1]-values[r-1,1])\
        #                             +(order-1)/2*values[r,order-3]*(values[r,2]-values[r-1,2]))
        #                 forward_differences[r-1] = np.abs((order-1)*(right_boundary_value[order-1]-values[r,order-1])\
        #                         +(order-1)*values[r,order-2]*(right_boundary_value[1]-values[r,1])\
        #                             +(order-1)/2*values[r,order-3]*(right_boundary_value[2]-values[r,2]))  

        # prev_order = order
        # prev_n_variables = n_variables
        # order = orders[-1]
        # n_variables = numbers_of_variables[-1]

        # l = r+1 

        # left_boundary_value = input_values[l-1,:]
        # # left_boundary_value[n_variables:] = values[l,n_variables:]              
        # left_boundary_value[n_variables:] = input_values[l,n_variables:]  

        # if order > 3:
        #     if increase_criterion_flags[l-1] == 0:
        #         if order == 4:
        #             backward_differences[l-1] = np.abs(6/values[l,0]*(values[l,3]-left_boundary_value[3])) 
        #             forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                    
        #         elif order == 5:
        #             backward_differences[l-1] = np.abs(4*(values[l,4]-left_boundary_value[4])\
        #                 +4*values[l,3]*(values[l,1]-left_boundary_value[1]))
        #             forward_differences[l-1] = np.abs(4*(values[l+1,4]-values[l,4])\
        #                 +4*values[l,3]*(values[l+1,1]-values[l,1]))
        #         else:
        #             backward_differences[l-1] = np.abs((order-1)*(values[l,order-1]-left_boundary_value[order-1])\
        #                     +(order-1)*values[l,order-2]*(values[l,1]-left_boundary_value[1])\
        #                         +(order-1)/2*values[l,order-3]*(values[l,2]-left_boundary_value[2]))
        #             forward_differences[l-1] = np.abs((order-1)*(values[l+1,order-1]-values[l,order-1])\
        #                     +(order-1)*values[l,order-2]*(values[l+1,1]-values[l,1])\
        #                         +(order-1)/2*values[l,order-3]*(values[l+1,2]-values[l,2]))                                     
        #     for i in range(l+1,r+1):
        #         if increase_criterion_flags[i-1] == 0:
        #             if order == 4:
        #                 backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
        #                 forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                    
        #             elif order == 5:
        #                 backward_differences[i-1] = np.abs(4*(values[i,4]-values[i-1,4])\
        #                     +4*values[i,3]*(values[i,1]-values[i-1,1]))
        #                 forward_differences[i-1] = np.abs(4*(values[i+1,4]-values[i,4])\
        #                     +4*values[i,3]*(values[i+1,1]-values[i,1]))
        #             else:
        #                 backward_differences[i-1] = np.abs((order-1)*(values[i,order-1]-values[i-1,order-1])\
        #                         +(order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
        #                             +(order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2]))
        #                 forward_differences[i-1] = np.abs((order-1)*(values[i+1,order-1]-values[i,order-1])\
        #                         +(order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
        #                             +(order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2])) 

        r = 0

        order = orders[0]
        n_variables = numbers_of_variables[0]
        for m in range(len(boundary_interfaces)):

            prev_order = order
            prev_n_variables = n_variables
            order = orders[m]
            n_variables = numbers_of_variables[m]
            next_order = orders[m+1]
            next_n_variables = numbers_of_variables[m+1]

            l = r+1
            r = boundary_interfaces[m]  

            left_boundary_value = input_values[l-1,:]
            right_boundary_value = input_values[r+1,:]
            # left_boundary_value[prev_n_variables:] = values[l,prev_n_variables:]
            # right_boundary_value[prev_n_variables:] = values[r,prev_n_variables:]
            left_boundary_value[n_variables:] = input_values[l,n_variables:]
            right_boundary_value[n_variables:] = input_values[r,n_variables:]

            if order > 3:
                if increase_criterion_flags[l-1] == 0:
                    if order == 4:
                        backward_differences[l-1] = np.abs(6/values[l,0]*(values[l,3]-left_boundary_value[3])) 
                        forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                    
                    elif order == 5:
                        backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-left_boundary_value[1]))
                        forward_differences[l-1] = np.abs(4*values[l,3]*(values[l+1,1]-values[l,1]))
                    else:
                        backward_differences[l-1] = np.abs((order-1)*values[l,order-2]*(values[l,1]-left_boundary_value[1])\
                                    +(order-1)/2*values[l,order-3]*(values[l,2]-left_boundary_value[2]))
                        forward_differences[l-1] = np.abs((order-1)*values[l,order-2]*(values[l+1,1]-values[l,1])\
                                    +(order-1)/2*values[l,order-3]*(values[l+1,2]-values[l,2]))                                     
                for i in range(l+1,r):
                    if increase_criterion_flags[i-1] == 0:
                        if order == 4:
                            backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
                            forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                    
                        elif order == 5:
                            backward_differences[i-1] = np.abs(4*values[i,3]*(values[i,1]-values[i-1,1]))
                            forward_differences[i-1] = np.abs(4*values[i,3]*(values[i+1,1]-values[i,1]))
                        else:
                            backward_differences[i-1] = np.abs((order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
                                        +(order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2]))
                            forward_differences[i-1] = np.abs((order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
                                        +(order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2]))             
                if increase_criterion_flags[r-1] == 0:
                    if order == 4:
                        backward_differences[r-1] = np.abs(6/values[r,0]*(values[r,3]-values[r-1,3])) 
                        forward_differences[r-1] = np.abs(6/values[r,0]*(right_boundary_value[3]-values[r,3]))                    
                    elif order == 5:
                        backward_differences[r-1] = np.abs(4*values[r,3]*(values[r,1]-values[r-1,1]))
                        forward_differences[r-1] = np.abs(4*values[r,3]*(right_boundary_value[1]-values[r,1]))
                    else:
                        backward_differences[r-1] = np.abs((order-1)*values[r,order-2]*(values[r,1]-values[r-1,1])\
                                    +(order-1)/2*values[r,order-3]*(values[r,2]-values[r-1,2]))
                        forward_differences[r-1] = np.abs((order-1)*values[r,order-2]*(right_boundary_value[1]-values[r,1])\
                                    +(order-1)/2*values[r,order-3]*(right_boundary_value[2]-values[r,2]))  

        prev_order = order
        prev_n_variables = n_variables
        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        l = r+1 

        left_boundary_value = input_values[l-1,:]
        # left_boundary_value[n_variables:] = values[l,n_variables:]              
        left_boundary_value[n_variables:] = input_values[l,n_variables:]  

        if order > 3:
            if increase_criterion_flags[l-1] == 0:
                if order == 4:
                    backward_differences[l-1] = np.abs(6/values[l,0]*(values[l,3]-left_boundary_value[3])) 
                    forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                    
                elif order == 5:
                    backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-left_boundary_value[1]))
                    forward_differences[l-1] = np.abs(4*values[l,3]*(values[l+1,1]-values[l,1]))
                else:
                    backward_differences[l-1] = np.abs((order-1)*values[l,order-2]*(values[l,1]-left_boundary_value[1])\
                                +(order-1)/2*values[l,order-3]*(values[l,2]-left_boundary_value[2]))
                    forward_differences[l-1] = np.abs((order-1)*values[l,order-2]*(values[l+1,1]-values[l,1])\
                                +(order-1)/2*values[l,order-3]*(values[l+1,2]-values[l,2]))                                     
            for i in range(l+1,n+1):
                if increase_criterion_flags[i-1] == 0:
                    if order == 4:
                        backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
                        forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                    
                    elif order == 5:
                        backward_differences[i-1] = np.abs(4*values[i,3]*(values[i,1]-values[i-1,1]))
                        forward_differences[i-1] = np.abs(4*values[i,3]*(values[i+1,1]-values[i,1]))
                    else:
                        backward_differences[i-1] = np.abs((order-1)*values[i,order-2]*(values[i,1]-values[i-1,1])\
                                    +(order-1)/2*values[i,order-3]*(values[i,2]-values[i-1,2]))
                        forward_differences[i-1] = np.abs((order-1)*values[i,order-2]*(values[i+1,1]-values[i,1])\
                                    +(order-1)/2*values[i,order-3]*(values[i+1,2]-values[i,2])) 

        # breakdown_estimators_decrease = (backward_differences+forward_differences)/delta_x
        breakdown_estimators_decrease = np.maximum(backward_differences,forward_differences)/delta_x

        for i in range(n):
            if breakdown_estimators_decrease[i] < tolerance_decrease: 
                decrease_criterion_flags[i] = -2

        return breakdown_estimators_decrease, decrease_criterion_flags

    def compute_breakdown_criteria_increase(self,
                                   values: np.ndarray,
                                #    padded_vectors_left,
                                #    padded_vectors_right,
                                   n: int,
                                   delta_x: float,
                                   delta_t: float,
                                   max_order: int,
                                   orders_cellwise: list,
                                   orders,
                                   numbers_of_variables,
                                   boundary_interfaces: list,
                                   tolerance_increase = 0.005) -> tuple[np.ndarray,np.ndarray]:       
        """
        TODO
        """
        increase_criterion_flags = np.zeros(n,dtype=int)
        breakdown_estimators_increase = np.zeros(n)

        backward_differences = np.zeros(n)
        forward_differences = np.zeros(n)

        input_values = np.copy(values)

        # r = 0

        # order = orders[0]
        # for m in range(len(boundary_interfaces)):
        #     prev_order = order
        #     order = orders[m]
        #     next_order = orders[m+1]
        #     n_variables = numbers_of_variables[m]

        #     l = r+1
        #     r = boundary_interfaces[m]                  

        #     if prev_order > order:
        #         if order == 2:
        #             forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                    
        #         elif order == 3:
        #             backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-values[l-1,1]))
        #             forward_differences[l-1] = np.abs(4*values[l,3]*(values[l+1,1]-values[l,1])+4*(values[l+1,4]-values[l,4]))
        #         elif order == 4:
        #             backward_differences[l-1] = np.abs(5*values[l,4]*(values[l,1]-values[l-1,1])+5/2*values[l,3]*(values[l,2]-values[l-1,2]))
        #             forward_differences[l-1] = np.abs(5*values[l,4]*(values[l+1,1]-values[l,1])+\
        #                                               5/2*values[l,3]*(values[l+1,2]-values[l,2])+5*(values[l+1,5]-values[l,5]))
        #         elif order == max_order:
        #             backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-values[l-1,1])\
        #                     +(order+1)/2*values[l,order-1]*(values[l,2]-values[l-1,2]))
        #             forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
        #                     +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2]))     
        #         else:
        #             backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-values[l-1,1])\
        #                     +(order+1)/2*values[l,order-1]*(values[l,2]-values[l-1,2]))
        #             forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
        #                     +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2])\
        #                     +(order+1)*(values[l+1,order+1]-values[l,order+1]))    
        #         l = l+1
        #     # else:
        #     #     l = l-1
        #     if next_order > order:
        #         r2 = r
        #         if order == 2:
        #             backward_differences[r-1] = np.abs(6/values[r,0]*(values[r,3]-values[r-1,3]))             
        #         elif order == 3:
        #             backward_differences[r-1] = np.abs(4*values[r,3]*(values[r,1]-values[r-1,1])+4*(values[r,4]-values[r-1,4]))
        #             forward_differences[r-1] = np.abs(4*values[r,3]*(values[r+1,1]-values[r,1]))
        #         elif order == 4:
        #             backward_differences[r-1] = np.abs(5*values[r,4]*(values[r,1]-values[r-1,1])+5/2*values[r,3]*(values[r,2]-values[r-1,2])\
        #                                                +5*(values[r,5]-values[r-1,5]))
        #             forward_differences[r-1] = np.abs(5*values[r,4]*(values[r+1,1]-values[r,1])+5/2*values[r,3]*(values[r+1,2]-values[r,2]))
        #         elif order == max_order:
        #             backward_differences[r-1] = np.abs((order+1)*values[r,order]*(values[r,1]-values[r-1,1])\
        #                     +(order+1)/2*values[r,order-1]*(values[r,2]-values[r-1,2]))
        #             forward_differences[r-1] = np.abs((order+1)*values[r,order]*(values[r+1,1]-values[r,1])\
        #                     +(order+1)/2*values[r,order-1]*(values[r+1,2]-values[r,2]))     
        #         else:
        #             backward_differences[r-1] = np.abs((order+1)*values[r,order]*(values[r,1]-values[r-1,1])\
        #                     +(order+1)/2*values[r,order-1]*(values[r,2]-values[r-1,2])\
        #                         +(order+1)*(values[r,order+1]-values[r-1,order+1]))
        #             forward_differences[r-1] = np.abs((order+1)*values[r,order]*(values[r+1,1]-values[r,1])\
        #                     +(order+1)/2*values[r,order-1]*(values[r+1,2]-values[r,2])) 
        #     else:
        #         r2 = r+1

        #     for i in range(l,r2):
        #         if order == 2:
        #             backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
        #             forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                
        #         elif order == 3:
        #             backward_differences[i-1] = np.abs(4*values[i,3]*(values[i,1]-values[i-1,1])+4*(values[i,4]-values[i-1,4]))
        #             forward_differences[i-1] = np.abs(4*values[i,3]*(values[i+1,1]-values[i,1])+4*(values[i+1,4]-values[i,4]))
        #         elif order == 4:
        #             backward_differences[i-1] = np.abs(5*values[i,4]*(values[i,1]-values[i-1,1])+\
        #                                                5/2*values[i,3]*(values[i,2]-values[i-1,2])+5*(values[i,5]-values[i-1,5]))
        #             forward_differences[i-1] = np.abs(5*values[i,4]*(values[i+1,1]-values[i,1])+\
        #                                               5/2*values[i,3]*(values[i+1,2]-values[i,2])+5*(values[i+1,5]-values[i,5]))
        #         elif order == max_order:
        #             backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
        #                     +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2]))
        #             forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
        #                     +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2]))     
        #         else:
        #             backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
        #                     +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2])\
        #                         +(order+1)*(values[i,order+1]-values[i-1,order+1]))
        #             forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
        #                     +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])\
        #                     +(order+1)*(values[i+1,order+1]-values[i,order+1]))           

        # order = orders[-1]
        # n_variables = numbers_of_variables[-1]

        # l = r+1
        # r = n                 

        # if orders[-2] > orders[-1]:
        #     if order == 2:
        #         forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                    
        #     elif order == 3:
        #         backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-values[l-1,1]))
        #         forward_differences[l-1] = np.abs(4*values[l,3]*(values[l+1,1]-values[l,1])+4*(values[l+1,4]-values[l,4]))
        #     elif order == 4:
        #         backward_differences[l-1] = np.abs(5*values[l,4]*(values[l,1]-values[l-1,1])+5/2*values[l,3]*(values[l,2]-values[l-1,2]))
        #         forward_differences[l-1] = np.abs(5*values[l,4]*(values[l+1,1]-values[l,1])+\
        #                                             5/2*values[l,3]*(values[l+1,2]-values[l,2])+5*(values[l+1,5]-values[l,5]))
        #     elif order == max_order:
        #         backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-values[l-1,1])\
        #                 +(order+1)/2*values[l,order-1]*(values[l,2]-values[l-1,2]))
        #         forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
        #                 +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2]))     
        #     else:
        #         backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-values[l-1,1])\
        #                 +(order+1)/2*values[l,order-1]*(values[l,2]-values[l-1,2]))
        #         forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
        #                 +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2])\
        #                 +(order+1)*(values[l+1,order+1]-values[l,order+1]))    
        #     l = l+1
        # # else:
        # #     l = l-1  

        # for i in range(l,r+1):
        #     if order == 2:
        #         backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
        #         forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                    
        #     elif order == 3:
        #         backward_differences[i-1] = np.abs(4*values[i,3]*(values[i,1]-values[i-1,1])+4*(values[i,4]-values[i-1,4]))
        #         forward_differences[i-1] = np.abs(4*values[i,3]*(values[i+1,1]-values[i,1])+4*(values[i+1,4]-values[i,4]))
        #     elif order == 4:
        #         backward_differences[i-1] = np.abs(5*values[i,4]*(values[i,1]-values[i-1,1])+\
        #                                             5/2*values[i,3]*(values[i,2]-values[i-1,2])+5*(values[i,5]-values[i-1,5]))
        #         forward_differences[i-1] = np.abs(5*values[i,4]*(values[i+1,1]-values[i,1])+\
        #                                             5/2*values[i,3]*(values[i+1,2]-values[i,2])+5*(values[i+1,5]-values[i,5]))
        #     elif order == max_order:
        #         backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
        #                 +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2]))
        #         forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
        #                 +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2]))     
        #     else:
        #         backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
        #                 +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2])\
        #                     +(order+1)*(values[i,order+1]-values[i-1,order+1]))
        #         forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
        #                 +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])\
        #                 +(order+1)*(values[i+1,order+1]-values[i,order+1]))      




        # r = 0

        # order = orders[0]
        # n_variables = numbers_of_variables[0]
        # for m in range(len(boundary_interfaces)):

        #     prev_order = order
        #     prev_n_variables = n_variables
        #     order = orders[m]
        #     n_variables = numbers_of_variables[m]
        #     next_order = orders[m+1]
        #     next_n_variables = numbers_of_variables[m+1]

        #     l = r+1
        #     r = boundary_interfaces[m]  

        #     left_boundary_value = input_values[l-1,:]
        #     right_boundary_value = input_values[r+1,:]
        #     # left_boundary_value[prev_n_variables:] = values[l,prev_n_variables:]
        #     # right_boundary_value[prev_n_variables:] = values[r,prev_n_variables:]
        #     left_boundary_value[n_variables:] = input_values[l,n_variables:]
        #     right_boundary_value[n_variables:] = input_values[r,n_variables:]

        #     if increase_criterion_flags[l-1] == 0:
        #         if order == 2:
        #             backward_differences[l-1] = np.abs(6/values[l,0]*(values[l,3]-left_boundary_value[3])) 
        #             forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                
        #         elif order == 3:
        #             backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-left_boundary_value[1])+4*(values[l,4]-left_boundary_value[4]))
        #             forward_differences[l-1] = np.abs(4*values[l,3]*(values[l+1,1]-values[l,1])+4*(values[l+1,4]-values[l,4]))
        #         elif order == 4:
        #             backward_differences[l-1] = np.abs(5*values[l,4]*(values[l,1]-left_boundary_value[1])+\
        #                                                5/2*values[l,3]*(values[l,2]-left_boundary_value[2])+5*(values[l,5]-left_boundary_value[5]))
        #             forward_differences[l-1] = np.abs(5*values[l,4]*(values[l+1,1]-values[l,1])+\
        #                                               5/2*values[l,3]*(values[l+1,2]-values[l,2])+5*(values[l+1,5]-values[l,5]))
        #         elif order == max_order:
        #             backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-left_boundary_value[1])\
        #                     +(order+1)/2*values[l,order-1]*(values[l,2]-left_boundary_value[2]))
        #             forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
        #                     +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2]))     
        #         else:
        #             backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-left_boundary_value[1])\
        #                     +(order+1)/2*values[l,order-1]*(values[l,2]-left_boundary_value[2])\
        #                         +(order+1)*(values[l,order+1]-left_boundary_value[order+1]))
        #             forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
        #                     +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2])\
        #                     +(order+1)*(values[l+1,order+1]-values[l,order+1]))    
                                                    
        #     for i in range(l+1,r):
        #         if order == 2:
        #             backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
        #             forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                
        #         elif order == 3:
        #             backward_differences[i-1] = np.abs(4*values[i,3]*(values[i,1]-values[i-1,1])+4*(values[i,4]-values[i-1,4]))
        #             forward_differences[i-1] = np.abs(4*values[i,3]*(values[i+1,1]-values[i,1])+4*(values[i+1,4]-values[i,4]))
        #         elif order == 4:
        #             backward_differences[i-1] = np.abs(5*values[i,4]*(values[i,1]-values[i-1,1])+\
        #                                                5/2*values[i,3]*(values[i,2]-values[i-1,2])+5*(values[i,5]-values[i-1,5]))
        #             forward_differences[i-1] = np.abs(5*values[i,4]*(values[i+1,1]-values[i,1])+\
        #                                               5/2*values[i,3]*(values[i+1,2]-values[i,2])+5*(values[i+1,5]-values[i,5]))
        #         elif order == max_order:
        #             backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
        #                     +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2]))
        #             forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
        #                     +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2]))     
        #         else:
        #             backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
        #                     +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2])\
        #                         +(order+1)*(values[i,order+1]-values[i-1,order+1]))
        #             forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
        #                     +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])\
        #                     +(order+1)*(values[i+1,order+1]-values[i,order+1]))             

        #     if order == 2:
        #         backward_differences[r-1] = np.abs(6/values[r,0]*(values[r,3]-values[r-1,3])) 
        #         forward_differences[r-1] = np.abs(6/values[r,0]*(right_boundary_value[3]-values[r,3]))                
        #     elif order == 3:
        #         backward_differences[r-1] = np.abs(4*values[r,3]*(values[r,1]-values[r-1,1])+4*(values[r,4]-values[r-1,4]))
        #         forward_differences[r-1] = np.abs(4*values[r,3]*(right_boundary_value[1]-values[r,1])+4*(right_boundary_value[4]-values[r,4]))
        #     elif order == 4:
        #         backward_differences[r-1] = np.abs(5*values[r,4]*(values[r,1]-values[r-1,1])+\
        #                                             5/2*values[r,3]*(values[r,2]-values[r-1,2])+5*(values[r,5]-values[r-1,5]))
        #         forward_differences[r-1] = np.abs(5*values[r,4]*(right_boundary_value[1]-values[r,1])+\
        #                                             5/2*values[r,3]*(right_boundary_value[2]-values[r,2])+5*(right_boundary_value[5]-values[r,5]))
        #     elif order == max_order:
        #         backward_differences[r-1] = np.abs((order+1)*values[r,order]*(values[r,1]-values[r-1,1])\
        #                 +(order+1)/2*values[r,order-1]*(values[r,2]-values[r-1,2]))
        #         forward_differences[r-1] = np.abs((order+1)*values[r,order]*(right_boundary_value[1]-values[r,1])\
        #                 +(order+1)/2*values[r,order-1]*(right_boundary_value[2]-values[r,2]))     
        #     else:
        #         backward_differences[r-1] = np.abs((order+1)*values[r,order]*(values[r,1]-values[r-1,1])\
        #                 +(order+1)/2*values[r,order-1]*(values[r,2]-values[r-1,2])\
        #                     +(order+1)*(values[r,order+1]-values[r-1,order+1]))
        #         forward_differences[r-1] = np.abs((order+1)*values[r,order]*(right_boundary_value[1]-values[r,1])\
        #                 +(order+1)/2*values[r,order-1]*(right_boundary_value[2]-values[r,2])\
        #                 +(order+1)*(right_boundary_value[order+1]-values[r,order+1]))   

        # prev_order = order
        # prev_n_variables = n_variables
        # order = orders[-1]
        # n_variables = numbers_of_variables[-1]

        # l = r+1 

        # left_boundary_value = input_values[l-1,:]
        # # left_boundary_value[n_variables:] = values[l,n_variables:]              
        # left_boundary_value[n_variables:] = input_values[l,n_variables:] 

        # if increase_criterion_flags[l-1] == 0:
        #     if order == 2:
        #         backward_differences[l-1] = np.abs(6/values[l,0]*(values[l,3]-left_boundary_value[3])) 
        #         forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                
        #     elif order == 3:
        #         backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-left_boundary_value[1])+4*(values[l,4]-left_boundary_value[4]))
        #         forward_differences[l-1] = np.abs(4*values[l,3]*(values[l+1,1]-values[l,1])+4*(values[l+1,4]-values[l,4]))
        #     elif order == 4:
        #         backward_differences[l-1] = np.abs(5*values[l,4]*(values[l,1]-left_boundary_value[1])+\
        #                                             5/2*values[l,3]*(values[l,2]-left_boundary_value[2])+5*(values[l,5]-left_boundary_value[5]))
        #         forward_differences[l-1] = np.abs(5*values[l,4]*(values[l+1,1]-values[l,1])+\
        #                                             5/2*values[l,3]*(values[l+1,2]-values[l,2])+5*(values[l+1,5]-values[l,5]))
        #     elif order == max_order:
        #         backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-left_boundary_value[1])\
        #                 +(order+1)/2*values[l,order-1]*(values[l,2]-left_boundary_value[2]))
        #         forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
        #                 +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2]))     
        #     else:
        #         backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-left_boundary_value[1])\
        #                 +(order+1)/2*values[l,order-1]*(values[l,2]-left_boundary_value[2])\
        #                     +(order+1)*(values[l,order+1]-left_boundary_value[order+1]))
        #         forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
        #                 +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2])\
        #                 +(order+1)*(values[l+1,order+1]-values[l,order+1]))    
                                                
        # for i in range(l+1,r):
        #     if order == 2:
        #         backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
        #         forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                
        #     elif order == 3:
        #         backward_differences[i-1] = np.abs(4*values[i,3]*(values[i,1]-values[i-1,1])+4*(values[i,4]-values[i-1,4]))
        #         forward_differences[i-1] = np.abs(4*values[i,3]*(values[i+1,1]-values[i,1])+4*(values[i+1,4]-values[i,4]))
        #     elif order == 4:
        #         backward_differences[i-1] = np.abs(5*values[i,4]*(values[i,1]-values[i-1,1])+\
        #                                             5/2*values[i,3]*(values[i,2]-values[i-1,2])+5*(values[i,5]-values[i-1,5]))
        #         forward_differences[i-1] = np.abs(5*values[i,4]*(values[i+1,1]-values[i,1])+\
        #                                             5/2*values[i,3]*(values[i+1,2]-values[i,2])+5*(values[i+1,5]-values[i,5]))
        #     elif order == max_order:
        #         backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
        #                 +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2]))
        #         forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
        #                 +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2]))     
        #     else:
        #         backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
        #                 +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2])\
        #                     +(order+1)*(values[i,order+1]-values[i-1,order+1]))
        #         forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
        #                 +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])\
        #                 +(order+1)*(values[i+1,order+1]-values[i,order+1]))  

        r = 0

        order = orders[0]
        n_variables = numbers_of_variables[0]
        for m in range(len(boundary_interfaces)):

            prev_order = order
            prev_n_variables = n_variables
            order = orders[m]
            n_variables = numbers_of_variables[m]
            next_order = orders[m+1]
            next_n_variables = numbers_of_variables[m+1]

            l = r+1
            r = boundary_interfaces[m]  

            left_boundary_value = input_values[l-1,:]
            right_boundary_value = input_values[r+1,:]
            # left_boundary_value[prev_n_variables:] = values[l,prev_n_variables:]
            # right_boundary_value[prev_n_variables:] = values[r,prev_n_variables:]
            left_boundary_value[n_variables:] = input_values[l,n_variables:]
            right_boundary_value[n_variables:] = input_values[r,n_variables:]

            if order == 2:
                backward_differences[l-1] = np.abs(6/values[l,0]*(values[l,3]-left_boundary_value[3])) 
                forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3]))                
            elif order == 3:
                backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-left_boundary_value[1]))
                forward_differences[l-1] = np.abs(4*values[l,3]*(values[l+1,1]-values[l,1]))
            elif order == 4:
                backward_differences[l-1] = np.abs(5*values[l,4]*(values[l,1]-left_boundary_value[1])+\
                                                    5/2*values[l,3]*(values[l,2]-left_boundary_value[2]))
                forward_differences[l-1] = np.abs(5*values[l,4]*(values[l+1,1]-values[l,1])+\
                                                    5/2*values[l,3]*(values[l+1,2]-values[l,2]))
            elif order == max_order:
                backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-left_boundary_value[1])\
                        +(order+1)/2*values[l,order-1]*(values[l,2]-left_boundary_value[2]))
                forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
                        +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2]))     
            else:
                backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-left_boundary_value[1])\
                        +(order+1)/2*values[l,order-1]*(values[l,2]-left_boundary_value[2]))
                forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
                        +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2]))    
                                                    
            for i in range(l+1,r):
                if order == 2:
                    backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
                    forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                
                elif order == 3:
                    backward_differences[i-1] = np.abs(4*values[i,3]*(values[i,1]-values[i-1,1]))
                    forward_differences[i-1] = np.abs(4*values[i,3]*(values[i+1,1]-values[i,1]))
                elif order == 4:
                    backward_differences[i-1] = np.abs(5*values[i,4]*(values[i,1]-values[i-1,1])+\
                                                       5/2*values[i,3]*(values[i,2]-values[i-1,2]))
                    forward_differences[i-1] = np.abs(5*values[i,4]*(values[i+1,1]-values[i,1])+\
                                                      5/2*values[i,3]*(values[i+1,2]-values[i,2]))
                elif order == max_order:
                    backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
                            +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2]))
                    forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
                            +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2]))     
                else:
                    backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
                            +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2]))
                    forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
                            +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2]))             

            if order == 2:
                backward_differences[r-1] = np.abs(6/values[r,0]*(values[r,3]-values[r-1,3])) 
                forward_differences[r-1] = np.abs(6/values[r,0]*(right_boundary_value[3]-values[r,3]))                
            elif order == 3:
                backward_differences[r-1] = np.abs(4*values[r,3]*(values[r,1]-values[r-1,1]))
                forward_differences[r-1] = np.abs(4*values[r,3]*(right_boundary_value[1]-values[r,1]))
            elif order == 4:
                backward_differences[r-1] = np.abs(5*values[r,4]*(values[r,1]-values[r-1,1])+\
                                                    5/2*values[r,3]*(values[r,2]-values[r-1,2]))
                forward_differences[r-1] = np.abs(5*values[r,4]*(right_boundary_value[1]-values[r,1])+\
                                                    5/2*values[r,3]*(right_boundary_value[2]-values[r,2]))
            elif order == max_order:
                backward_differences[r-1] = np.abs((order+1)*values[r,order]*(values[r,1]-values[r-1,1])\
                        +(order+1)/2*values[r,order-1]*(values[r,2]-values[r-1,2]))
                forward_differences[r-1] = np.abs((order+1)*values[r,order]*(right_boundary_value[1]-values[r,1])\
                        +(order+1)/2*values[r,order-1]*(right_boundary_value[2]-values[r,2]))     
            else:
                backward_differences[r-1] = np.abs((order+1)*values[r,order]*(values[r,1]-values[r-1,1])\
                        +(order+1)/2*values[r,order-1]*(values[r,2]-values[r-1,2]))
                forward_differences[r-1] = np.abs((order+1)*values[r,order]*(right_boundary_value[1]-values[r,1])\
                        +(order+1)/2*values[r,order-1]*(right_boundary_value[2]-values[r,2]))   

        prev_order = order
        prev_n_variables = n_variables
        order = orders[-1]
        n_variables = numbers_of_variables[-1]

        l = r+1 

        left_boundary_value = input_values[l-1,:]
        # left_boundary_value[n_variables:] = values[l,n_variables:]              
        left_boundary_value[n_variables:] = input_values[l,n_variables:] 

        if order == 2:
            left_boundary_value[3] = values[l,3]
            backward_differences[l-1] = np.abs(6/values[l,0]*(values[l,3]-left_boundary_value[3])) 
            forward_differences[l-1] = np.abs(6/values[l,0]*(values[l+1,3]-values[l,3])) 
            print(values[boundary_interfaces[-1]+1:,3])             
        elif order == 3:
            backward_differences[l-1] = np.abs(4*values[l,3]*(values[l,1]-left_boundary_value[1]))
            forward_differences[l-1] = np.abs(4*values[l,3]*(values[l+1,1]-values[l,1]))
        elif order == 4:
            backward_differences[l-1] = np.abs(5*values[l,4]*(values[l,1]-left_boundary_value[1])+\
                                                5/2*values[l,3]*(values[l,2]-left_boundary_value[2]))
            forward_differences[l-1] = np.abs(5*values[l,4]*(values[l+1,1]-values[l,1])+\
                                                5/2*values[l,3]*(values[l+1,2]-values[l,2]))
        elif order == max_order:
            backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-left_boundary_value[1])\
                    +(order+1)/2*values[l,order-1]*(values[l,2]-left_boundary_value[2]))
            forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
                    +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2]))     
        else:
            backward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l,1]-left_boundary_value[1])\
                    +(order+1)/2*values[l,order-1]*(values[l,2]-left_boundary_value[2]))
            forward_differences[l-1] = np.abs((order+1)*values[l,order]*(values[l+1,1]-values[l,1])\
                    +(order+1)/2*values[l,order-1]*(values[l+1,2]-values[l,2]))    
                                                
        for i in range(l+1,n+1):
            if order == 2:
                backward_differences[i-1] = np.abs(6/values[i,0]*(values[i,3]-values[i-1,3])) 
                forward_differences[i-1] = np.abs(6/values[i,0]*(values[i+1,3]-values[i,3]))                
            elif order == 3:
                backward_differences[i-1] = np.abs(4*values[i,3]*(values[i,1]-values[i-1,1]))
                forward_differences[i-1] = np.abs(4*values[i,3]*(values[i+1,1]-values[i,1]))
            elif order == 4:
                backward_differences[i-1] = np.abs(5*values[i,4]*(values[i,1]-values[i-1,1])+\
                                                    5/2*values[i,3]*(values[i,2]-values[i-1,2]))
                forward_differences[i-1] = np.abs(5*values[i,4]*(values[i+1,1]-values[i,1])+\
                                                    5/2*values[i,3]*(values[i+1,2]-values[i,2]))
            elif order == max_order:
                backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
                        +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2]))
                forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
                        +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2]))     
            else:
                backward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i,1]-values[i-1,1])\
                        +(order+1)/2*values[i,order-1]*(values[i,2]-values[i-1,2]))
                forward_differences[i-1] = np.abs((order+1)*values[i,order]*(values[i+1,1]-values[i,1])\
                        +(order+1)/2*values[i,order-1]*(values[i+1,2]-values[i,2])) 

        # breakdown_estimators_increase = (forward_differences+backward_differences)/delta_x
        breakdown_estimators_increase = np.maximum(forward_differences,backward_differences)/delta_x
        for i in range(n):
            if breakdown_estimators_increase[i] > tolerance_increase: 
                increase_criterion_flags[i] = 2            

        return breakdown_estimators_increase, increase_criterion_flags