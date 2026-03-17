from abc import ABC, abstractmethod
import numpy as np
import scipy

#TODO: implement MomentModel as a subclass of PDE and include the possibility of simulating PDEs that are not moment models (and don't have an order)
class PDE(ABC):
    """
    This interface represents a partial differential equation.

    ...

    Attributes
    ----------
    initial_condition : str
        initial condition for the partial differential equation

    
    Abstract methods
    -------
    def compute_system_matrix(self,order,values):
        computes the system matrix of the partial differential equation evaluated in the given values, for the given order.
    def compute_source_term(self,order,values,delta_t):
        computes the source term of the partial differential equation evaluated in the given values, for the given order.
    def get_initial_values(self,order,initial_condition,position):
        calculates the initial values for one specific physical position
    def compute_number_of_variables(self,order):
        computes the number of state variables in the PDE given the order of the moment model
    def compute_all_breakdown_criteria(self,values,number_of_variables,n,delta_x,tolerance_up_height_gradient,tolerance_down_height_gradient,tolerance_up_momentum_gradient,tolerance_down_momentum_gradient,tolerance_up_last_moment,tolerance_down_last_moment)
        computes the values of all breakdown criteria in each mesh cell
    def compute_breakdown_criterion(self,values,number_of_variables,breakdown_criterion,n)
        computes the values of the given breakdown criterion in each mesh cell
    """

    @abstractmethod
    def __init__(self, 
                 initial_condition: str):
        """
        Constructs all the necessary attributes for the PDE object.

        Parameters
        ----------
        initial_condition : str
            initial condition of the PDE
        """

        pass

    @abstractmethod
    def compute_system_matrix(self,
                              order: int,
                              values: np.array) -> np.array:
        """
        Computes the system matrix with a given order of the PDE evaluated in the given values.

        Parameters
        ----------
        order : int
            order of the moment model PDE (TODO: create MomentModel as a subclass of PDE)
        values : numpy 1D array
            values of the variables

        
        Returns
        -------
        A: numpy 2D array
            System matrix

        """


        pass
    
    @abstractmethod
    def compute_source_term(self,
                            order: int,
                            values: np.array) -> np.array:
        """
        Computes the source term with a given order of the PDE evaluated in the given values.

        Parameters
        ----------
        order : int
            order of the moment model PDE (TODO: create MomentModel as a subclass of PDE)
        values : numpy 1D array
            values of the variables
        
        
        Returns
        -------
        S: numpy 1D array
            source term vector

        """

        pass

    @abstractmethod
    def get_initial_values(self,
                           order: int,
                           initial_condition: str,
                           position: float) -> np.array:

        """
        calculates the initial values for one specific physical position

        Parameters
        ----------
        order : int
            order of the moment model PDE (TODO: create MomentModel as a subclass of PDE)
        initial condition : str
            name of the initial condition
        position : float (if 1D) or numpy 1D array of floats (2D)
            the physical position in which the initial values are computed
        
        
        Returns
        -------
        initial_values: numpy 1D array
            initial values for the given initial condition evaluated in the phyiscal position

        """

        pass

    @abstractmethod
    def compute_number_of_variables(self,
                                    order: int) -> int:

        """
        given the order of the moment model expansion, compute the number of state variables in the PDE

        Parameters
        ----------
        order : int
            order of the moment model PDE (TODO: create MomentModel as a subclass of PDE)
        
        
        Returns
        -------
        number_of_variables: int
            number of state variables in the PDE

        """

        pass

    
    @abstractmethod
    def compute_all_breakdown_criteria(self,
                                   values: np.array,
                                   number_of_variables: list,
                                   n,
                                   delta_x,
                                   tolerance_up_height_gradient,
                                   tolerance_down_height_gradient,
                                   tolerance_up_momentum_gradient,
                                   tolerance_down_momentum_gradient,
                                   tolerance_up_last_moment,
                                   tolerance_down_last_moment) -> np.array:
        pass

    @abstractmethod
    def compute_breakdown_criterion(self,
                                   values: list,
                                   number_of_variables: int,
                                   breakdown_criterion: str,
                                   n: int) -> np.array:

        """
        Compute specific breakdown criterion value for quantifying the required modelling complexity

        Parameters
        ----------
        values : list of numpy 1D arrays
            the values of the variables in each mesh cell
        number_of_variables : integer
            the number of variables
        breakdown_criterion : str
            the breakdon criterion that is considered
        n : int
            the number of grid cells
        
        Returns
        -------
        breakdown_criterion_values: np.array
            value for the breakdown criterion in each mesh cell

        """

        pass


class SWME1D(PDE):

    """
    This class represents the one-dimensional Shallow Water Moment Equations (SWME1D).

    ...

    Attributes
    ----------
    initial_condition : str
        initial condition for the SWME1D
    viscosity : float
        value for the dynamic viscosity
    slip_length : float
        value for the slip length
    hyperbolic : boolean
        whether the model is hyperbolic, true (HSWME) or false (SWME)

    
    Implemented methods from interface PDE
    ---------------------------------
    def compute_system_matrix(self,order,values):
        computes the system matrix of the SWME1D evaluated in the given values, for the given order. 
    def compute_source_term(self,order,values):
        computes the source term of the SWME1D evaluated in the given values, for the given order.
    def get_initial_values(self,order,initial_condition,position):
        calculates the initial values for one specific physical position
    def compute_number_of_variables(self,order):
        computes the number of state variables in the PDE given the order of the moment model
    def compute_all_breakdown_criteria(self,values,orders,number_of_variables,n,delta_x,tolerance_up_height_gradient,
                                       tolerance_down_height_gradient,tolerance_up_momentum_gradient,tolerance_down_momentum_gradient,
                                       tolerance_up_moment_gradient,tolerance_down_moment_gradient,tolerance_up_last_moment,
                                       tolerance_down_last_moment,tolerance_up_source,tolerance_down_source)
        computes the values of all breakdown criterion in each mesh cell
    def compute_breakdown_criterion(self,values,breakdown_criterion,n)
        computes the values of the given breakdown criterion in each mesh cell

    Instance methods
    ----------------
    def compute_system_matrix_diff(self,order_low,values):
        calcultates the system matrix diff
    def compute_system_matrix_last_row(self,order,values):
        calculates the system matrix last row
    def _compute_source_matrix_inverse(self,order,values,delta_t):
        calculates the source matrix inverse if the rhs is linear for speed up
    def compute_source_term_lastentry(self,order,values,last_moment_zero):
        calculates the source term last entry
    def compute_vertical_velocity_profile(self,values):
        reconstruct the vertical velocity profiles from the moment values
    """

    def __init__(self, 
                 initial_condition: str,
                 viscosity: float,
                 slip_length: float,
                 hyperbolic: bool,
                 linear_source: bool):
        """
        Constructs all the necessary attributes for the SWME1D object.

        Parameters
        ----------
        initial_condition : str
            initial condition of the PDE
        viscosity : float
            dynamic viscosity value
        slip_length : float
            slip length value
        hyperbolic : boolean
            true if hyperbolic, false if not hyperbolic
        """
        self.initial_condition = initial_condition
        self.viscosity = viscosity
        self.slip_length = slip_length
        self.hyperbolic = hyperbolic
        self.linear_source = linear_source

    def compute_system_matrix(self,
                              order: int,
                              values: np.array,
                              **kwargs) -> np.array:
        
        g = kwargs["g"] if "g" in kwargs else 1
        A = np.zeros((order+2,order+2)) 
        h = values[0]
        um = values[1]/values[0]
        
        if order == 0:
            A[0][0] = 0
            A[0][1] = 1
            A[1][0] = g*h - um*um
            A[1][1] = 2.*um
        
        elif order == 1:
            alpha1 = values[2]/values[0]

            A[0][0] = 0
            A[0][1] = 1
            A[0][2] = 0
            A[1][0] = g*h - um*um - alpha1*alpha1/3
            A[1][1] = 2*um
            A[1][2] = 2*alpha1/3
            A[2][0] = -2*um*alpha1
            A[2][1] = 2*alpha1
            A[2][2] = um
        
        elif order == 2:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]

            if self.hyperbolic:
                alpha2 = 0

            A[0][0] = 0
            A[0][1] = 1
            A[0][2] = 0
            A[0][3] = 0
            A[1][0] = g*h - um*um - alpha1*alpha1/3 - alpha2*alpha2/5
            A[1][1] = 2*um
            A[1][2] = 2*alpha1/3
            A[1][3] = 2*alpha2/5
            A[2][0] = -2/5*alpha1*(5*um + 2*alpha2)
            A[2][1] = 2*alpha1
            A[2][2] = um + alpha2
            A[2][3] = 3*alpha1/5
            A[3][0] = -2/21*(7*alpha1*alpha1 + 3*alpha2*(7*um + alpha2))
            A[3][1] = 2*alpha2
            A[3][2] = alpha1/3
            A[3][3] = um + 3/7*alpha2
        
        elif order == 3:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0

            A[0][0] = 0
            A[0][1] = 1
            A[0][2] = 0
            A[0][3] = 0
            A[0][4] = 0
            A[1][0] = g*h - um*um - alpha1*alpha1/3 - alpha2*alpha2/5 - alpha3*alpha3/7
            A[1][1] = 2*um
            A[1][2] = 2*alpha1/3
            A[1][3] = 2*alpha2/5
            A[1][4] = 2*alpha3/7
            A[2][0] = -2/35*(7*alpha1*(5*um + 2*alpha2) + 9*alpha2*alpha3)
            A[2][1] = 2*alpha1
            A[2][2] = um + alpha2
            A[2][3] = 3*(alpha1 + alpha3)/5
            A[2][4] = 3*alpha2/7
            A[3][0] = -2/21*(3*alpha2*(7*um + alpha2)+(alpha1 + alpha3)*(7*alpha1 + 2*alpha3))
            A[3][1] = 2*alpha2
            A[3][2] = alpha1/3 + 9*alpha3/7
            A[3][3] = um + 3/7*alpha2
            A[3][4] = 4*alpha1/7 + alpha3/3
            A[4][0] = -2*um*alpha3 - 2*alpha2*(9*alpha1 + 4*alpha3)/15 
            A[4][1] = 2*alpha3
            A[4][2] = 0
            A[4][3] = 2*(alpha1 + alpha3)/5
            A[4][4] = um + alpha2/3
        
        elif order == 4:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0
                alpha4 = 0

            A[0][0] = 0
            A[0][1] = 1
            A[0][2] = 0
            A[0][3] = 0
            A[0][4] = 0
            A[0][5] = 0
            A[1][0] = g*h - um*um - alpha1*alpha1/3. - alpha2*alpha2/5. - \
            alpha3*alpha3/7. - alpha4*alpha4/9.
            A[1][1] = 2*um
            A[1][2] = (2*alpha1)/3.
            A[1][3] = (2*alpha2)/5.
            A[1][4] = (2*alpha3)/7.
            A[1][5] = (2*alpha4)/9.
            A[2][0] = (-2*(21*alpha1*(5*um + 2*alpha2) + alpha3*(27*alpha2 + \
            20*alpha4)))/105.
            A[2][1] = 2*alpha1
            A[2][2] = um + alpha2
            A[2][3] = (3*(alpha1 + alpha3))/5.
            A[2][4] = (3*(alpha2 + alpha4))/7.
            A[2][5] = alpha3/3.
            A[3][0] = (-2*(99*alpha2*alpha2 + 33*(alpha1 + alpha3)*(7*alpha1 + \
            2*alpha3) + 50*alpha4*alpha4 + 99*alpha2*(7*um + 2*alpha4)))/693.
            A[3][1] = 2*alpha2
            A[3][2] = alpha1/3. + (9*alpha3)/7.
            A[3][3] = um + (3*alpha2)/7. + (16*alpha4)/21.
            A[3][4] = (4*alpha1)/7. + alpha3/3.
            A[3][5] = (3*alpha2)/7. + (185*alpha4)/693.
            A[4][0] = alpha1*((-6*alpha2)/5. - (8*alpha4)/9.) - (2*alpha3*(165*um \
            + 44*alpha2 + 30*alpha4))/165.
            A[4][1] = 2*alpha3
            A[4][2] = (14*alpha4)/9.
            A[4][3] = (2*(alpha1 + alpha3))/5.
            A[4][4] = um + (alpha2 + alpha4)/3.
            A[4][5] = (5*alpha1)/9. + (3*alpha3)/11.
            A[5][0] = -2*um*alpha4 - (2*(1287*alpha2*alpha2 + \
            65*alpha3*(44*alpha1 + 9*alpha3) + 1300*alpha2*alpha4 + \
            405*alpha4*alpha4))/5005.
            A[5][1] = 2*alpha4
            A[5][2] = (-2*alpha3)/7.
            A[5][3] = (6*alpha2)/35. + (30*alpha4)/77.
            A[5][4] = (3*alpha1)/7. + (3*alpha3)/11.
            A[5][5] = um + (23*alpha2)/77. + (243*alpha4)/1001.
        
        elif order == 5:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0
                alpha4 = 0
                alpha5 = 0

            A[0][0] = 0
            A[0][1] = 1
            A[0][2] = 0
            A[0][3] = 0
            A[0][4] = 0
            A[0][5] = 0
            A[0][6] = 0
            A[1][0] = g*h - um*um - alpha1*alpha1/3. - alpha2*alpha2/5. - \
            alpha3*alpha3/7. - alpha4*alpha4/9. - alpha5*alpha5/11.
            A[1][1] = 2*um
            A[1][2] = (2*alpha1)/3.
            A[1][3] = (2*alpha2)/5.
            A[1][4] = (2*alpha3)/7.
            A[1][5] = (2*alpha4)/9.
            A[1][6] = (2*alpha5)/11.
            A[2][0] = (-2*alpha1*(5*um + 2*alpha2))/5. - (18*alpha2*alpha3)/35. - \
            (2*alpha4*(44*alpha3 + 35*alpha5))/231.
            A[2][1] = 2*alpha1
            A[2][2] = um + alpha2
            A[2][3] = (3*(alpha1 + alpha3))/5.
            A[2][4] = (3*(alpha2 + alpha4))/7.
            A[2][5] = (alpha3 + alpha5)/3.
            A[2][6] = (3*alpha4)/11.
            A[3][0] = (-2*(429*(alpha1 + alpha3)*(7*alpha1 + 2*alpha3) + \
            13*(99*alpha2*alpha2 + 50*alpha4*alpha4 + 99*alpha2*(7*um + \
            2*alpha4)) + 1950*alpha3*alpha5 + 525*alpha5*alpha5))/9009.
            A[3][1] = 2*alpha2
            A[3][2] = alpha1/3. + (9*alpha3)/7.
            A[3][3] = um + (3*alpha2)/7. + (16*alpha4)/21.
            A[3][4] = (4*alpha1)/7. + alpha3/3. + (125*alpha5)/231.
            A[3][5] = (3*alpha2)/7. + (185*alpha4)/693.
            A[3][6] = (80*alpha3)/231. + (95*alpha5)/429.
            A[4][0] = -2*um*alpha3 + alpha1*((-6*alpha2)/5. - (8*alpha4)/9.) - \
            (4*alpha4*(13*alpha3 + 10*alpha5))/143. - (4*alpha2*(22*alpha3 + \
            25*alpha5))/165.
            A[4][1] = 2*alpha3
            A[4][2] = (14*alpha4)/9.
            A[4][3] = (2*(alpha1 + alpha3))/5. + (10*alpha5)/11.
            A[4][4] = um + (alpha2 + alpha4)/3.
            A[4][5] = (5*alpha1)/9. + (3*(alpha3 + alpha5))/11.
            A[4][6] = (14*(13*alpha2 + 7*alpha4))/429.
            A[5][0] = (-2*(1287*alpha2*alpha2 + 5005*um*alpha4 + \
            1300*alpha2*alpha4 + 405*alpha4*alpha4 + 45*(alpha3 + \
            alpha5)*(13*alpha3 + 7*alpha5) + 65*alpha1*(44*alpha3 + \
            35*alpha5)))/5005.
            A[5][1] = 2*alpha4
            A[5][2] = (-2*alpha3)/7. + (20*alpha5)/11.
            A[5][3] = (6*alpha2)/35. + (30*alpha4)/77.
            A[5][4] = (3*(143*alpha1 + 91*alpha3 + 115*alpha5))/1001.
            A[5][5] = um + (23*alpha2)/77. + (243*alpha4)/1001.
            A[5][6] = (6*(91*alpha1 + 41*alpha3 + 35*alpha5))/1001.
            A[6][0] = (-2*(819*um*alpha5 + 30*alpha2*(13*alpha3 + 7*alpha5) + \
            alpha4*(455*alpha1 + 180*alpha3 + 126*alpha5)))/819.
            A[6][1] = 2*alpha5
            A[6][2] = (-5*alpha4)/9.
            A[6][3] = (5*alpha5)/13.
            A[6][4] = (5*(alpha2 + alpha4))/21.
            A[6][5] = (4*alpha1)/9. + (3*(alpha3 + alpha5))/13.
            A[6][6] = um + (11*alpha2)/39. + (8*alpha4)/39.
        
        elif order == 6:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]
            alpha6 = values[7]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0
                alpha4 = 0
                alpha5 = 0
                alpha6 = 0

            A[0][0] = 0
            A[0][1] = 1
            A[0][2] = 0
            A[0][3] = 0
            A[0][4] = 0
            A[0][5] = 0
            A[0][6] = 0
            A[0][7] = 0
            A[1][0] = g*h - um*um - alpha1*alpha1/3. - alpha2*alpha2/5. - \
            alpha3*alpha3/7. - alpha4*alpha4/9. - alpha5*alpha5/11. - \
            alpha6*alpha6/13.
            A[1][1] = 2*um
            A[1][2] = (2*alpha1)/3.
            A[1][3] = (2*alpha2)/5.
            A[1][4] = (2*alpha3)/7.
            A[1][5] = (2*alpha4)/9.
            A[1][6] = (2*alpha5)/11.
            A[1][7] = (2*alpha5)/13.
            A[2][0] = (-2*alpha1*(5*um + 2*alpha2))/5. - (18*alpha2*alpha3)/35. - \
            (8*alpha3*alpha4)/21. - (2*alpha5*(65*alpha4 + 54*alpha5))/429.
            A[2][1] = 2*alpha1
            A[2][2] = um + alpha2
            A[2][3] = (3*(alpha1 + alpha3))/5.
            A[2][4] = (3*(alpha2 + alpha4))/7.
            A[2][5] = (alpha3 + alpha5)/3.
            A[2][6] = (3*(alpha4 + alpha5))/11.
            A[2][7] = (3*alpha5)/13.
            A[3][0] = (-2*(1287*alpha2*alpha2 + 429*(alpha1 + alpha3)*(7*alpha1 + \
            2*alpha3) + 1287*alpha2*(7*um + 2*alpha4) + 1950*alpha3*alpha5 + \
            525*alpha5*alpha5 + (10*alpha4 + 21*alpha5)*(65*alpha4 + \
            21*alpha5)))/9009.
            A[3][1] = 2*alpha2
            A[3][2] = alpha1/3. + (9*alpha3)/7.
            A[3][3] = um + (3*alpha2)/7. + (16*alpha4)/21.
            A[3][4] = (4*alpha1)/7. + alpha3/3. + (125*alpha5)/231.
            A[3][5] = (3*alpha2)/7. + (185*alpha4)/693. + (60*alpha5)/143.
            A[3][6] = (80*alpha3)/231. + (95*alpha5)/429.
            A[3][7] = (125*alpha4 + 81*alpha5)/429.
            A[4][0] = (-2*(143*alpha1*(27*alpha2 + 20*alpha4) + \
            78*alpha2*(22*alpha3 + 25*alpha5) + 15*alpha5*(60*alpha4 + 49*alpha5) \
            + 15*alpha3*(429*um + 78*alpha4 + 100*alpha5)))/6435.
            A[4][1] = 2*alpha3
            A[4][2] = (14*alpha4)/9.
            A[4][3] = (2*(alpha1 + alpha3))/5. + (10*alpha5)/11.
            A[4][4] = um + alpha2/3. + alpha4/3. + (25*alpha5)/39.
            A[4][5] = (5*alpha1)/9. + (3*(alpha3 + alpha5))/11.
            A[4][6] = (14*(13*alpha2 + 7*(alpha4 + alpha5)))/429.
            A[4][7] = (2*(25*alpha3 + 14*alpha5))/143.
            A[5][0] = (-2*(21879*alpha2*alpha2 + 6885*alpha4*alpha4 + 765*(alpha3 \
            + alpha5)*(13*alpha3 + 7*alpha5) + 1105*alpha1*(44*alpha3 + \
            35*alpha5) + 4410*alpha6*alpha6 + 595*alpha4*(143*um + 20*alpha5) + \
            425*alpha2*(52*alpha4 + 63*alpha5)))/85085.
            A[5][1] = 2*alpha4
            A[5][2] = (-2*alpha3)/7. + (20*alpha5)/11.
            A[5][3] = (6*(143*alpha2 + 325*alpha4 + 875*alpha5))/5005.
            A[5][4] = (3*(143*alpha1 + 91*alpha3 + 115*alpha5))/1001.
            A[5][5] = um + (299*alpha2 + 243*alpha4 + 287*alpha5)/1001.
            A[5][6] = (6*(91*alpha1 + 41*alpha3 + 35*alpha5))/1001.
            A[5][7] = (6*(85*(2*alpha2 + alpha4) + 74*alpha5))/2431.
            A[6][0] = (-2*(510*alpha2*(13*alpha3 + 7*alpha5) + \
            51*alpha3*(60*alpha4 + 49*alpha5) + 119*alpha1*(65*alpha4 + \
            54*alpha5) + 21*alpha5*(663*um + 102*alpha4 + 80*alpha5)))/13923.
            A[6][1] = 2*alpha5
            A[6][2] = (-5*alpha4)/9. + (27*alpha5)/13.
            A[6][3] = (5*alpha5)/13.
            A[6][4] = (5*(alpha2 + alpha4))/21. + (14*alpha5)/39.
            A[6][5] = (4*alpha1)/9. + (3*(alpha3 + alpha5))/13.
            A[6][6] = um + (11*alpha2)/39. + (8*(alpha4 + alpha5))/39.
            A[6][7] = (119*alpha1 + 51*alpha3 + 40*alpha5)/221.
            A[7][0] = (-100*alpha3*alpha3)/231. - (14*alpha3*alpha5)/33. - \
            2*um*alpha5 - (2*alpha2*(25*alpha4 + 14*alpha5))/55. - \
            (4*(19*(85*alpha4*alpha4 + 459*alpha1*alpha5 + 60*alpha5*alpha5) + \
            2394*alpha4*alpha5 + 900*alpha6*alpha6))/31977.
            A[7][1] = 2*alpha5
            A[7][2] = (-9*alpha5)/11.
            A[7][3] = (-5*alpha4)/33. + (21*alpha5)/55.
            A[7][4] = (25*alpha3 + 49*alpha5)/231.
            A[7][5] = (3*alpha2)/11. + (19*alpha4)/99. + (42*alpha5)/187.
            A[7][6] = (255*alpha1 + 119*alpha3 + 104*alpha5)/561.
            A[7][7] = um + (3*alpha2)/11. + (104*alpha4)/561. + \
            (600*alpha5)/3553.
        
        else:
            print("This order is not implemented for the system matrix of the SWME1D-PDE.")

        return A

    def compute_system_matrix_diff(self,
                                   order_low: int,
                                   values: np.array) -> np.array:
        
        g = kwargs["g"] if "g" in kwargs else 1
        A_diff = np.zeros((order_low+2,order_low+2)) 
        h = values[0]
        um = values[1]/values[0]
        
        if order_low == 0:
            alpha1 = values[2]/values[0]
            A_diff[1][0] = -1/3.*alpha1*alpha1
        
        elif order_low == 1:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]

            if self.hyperbolic:
                alpha2 = 0

            A_diff[1][0] = -0.2*alpha2*alpha2
            A_diff[2][0] = (-4*alpha1*alpha2)/5.
            A_diff[2][2] = alpha2

        elif order_low == 2:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0

            A_diff[1][0] = -1./7.*alpha3*alpha3
            A_diff[2][0] = (-18*alpha2*alpha3)/35.
            A_diff[2][3] = (3*alpha3)/5.
            A_diff[3][0] = (-2*alpha3*(9*alpha1 + 2*alpha3))/21.
            A_diff[3][2] = (9*alpha3)/7.
        
        elif order_low == 3:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0
                alpha4 = 0

            A_diff[1][0] = -1./9.*alpha4*alpha4
            A_diff[2][0] = (-8*alpha3*alpha4)/21.
            A_diff[2][4] = (3*alpha4)/7.
            A_diff[3][0] = (-4*alpha4*(99*alpha2 + 25*alpha4))/693.
            A_diff[3][3] = (16*alpha4)/21.
            A_diff[4][0] = (-4*(22*alpha1 + 9*alpha3)*alpha4)/99.
            A_diff[4][2] = (14*alpha4)/9.
            A_diff[4][4] = alpha4/3.
        
        elif order_low == 4:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0
                alpha4 = 0
                alpha5 = 0

            A_diff[1][0] = -1/11.*alpha5*alpha5
            A_diff[2][0] = (-10*alpha4*alpha5)/33.
            A_diff[2][5] = alpha5/3.
            A_diff[3][0] = (-50*alpha5*(26*alpha3 + 7*alpha5))/3003.
            A_diff[3][4] = (125*alpha5)/231.
            A_diff[4][0] = (-20*(13*alpha2 + 6*alpha4)*alpha5)/429.
            A_diff[4][3] = (10*alpha5)/11.
            A_diff[4][5] = (3*alpha5)/11.
            A_diff[5][0] = (-2*alpha5*(455*alpha1 + 180*alpha3 + \
            63*alpha5))/1001.
            A_diff[5][2] = (20*alpha5)/11.
            A_diff[5][4] = (345*alpha5)/1001.
        
        elif order_low == 5:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]
            alpha6 = values[7]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0
                alpha4 = 0
                alpha5 = 0
                alpha6 = 0

            A_diff[1][0] = -1/13.*alpha6*alpha6
            A_diff[2][0] = (-36*alpha5*alpha6)/143.
            A_diff[2][6] = (3*alpha6)/11.
            A_diff[3][0] = (-2*alpha6*(25*alpha4 + 7*alpha6))/143.
            A_diff[3][5] = (60*alpha6)/143.
            A_diff[4][0] = (-2*(100*alpha3 + 49*alpha5)*alpha6)/429.
            A_diff[4][4] = (25*alpha6)/39.
            A_diff[4][6] = (98*alpha6)/429.
            A_diff[5][0] = (-2*alpha6*(765*alpha2 + 340*alpha4 + 126*alpha6))/2431.
            A_diff[5][3] = (150*alpha6)/143.
            A_diff[5][5] = (41*alpha6)/143.
            A_diff[6][0] = (-2*(306*alpha1 + 119*alpha3 + 80*alpha5)*alpha6)/663.
            A_diff[6][2] = (27*alpha6)/13.
            A_diff[6][4] = (14*alpha6)/39.
            A_diff[6][6] = (8*alpha6)/39.
        
        else:
            print("This order is not implemented for the system matrix diff of the SWME1D-PDE.")

        return A_diff
    
    def compute_system_matrix_last_row(self,
                              order: int,
                              values: np.array,
                              **kwargs) -> np.array:

        g = kwargs["g"] if "g" in kwargs else 1
        A_last_row=np.zeros((1,order+2)) 
        h = values[0]
        um = values[1]/values[0]
        
        if order == 0:
            A_last_row[0][0] = 0
            A_last_row[0][1] = 0
        
        elif order == 1:
            alpha1 = values[2]/values[0]

            A_last_row[0][0] = (-2*alpha1*alpha1)/3.
            A_last_row[0][1] = 0
            A_last_row[0][2] = alpha1/3.

        elif order == 2:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]

            if self.hyperbolic:
                alpha2 = 0

            A_last_row[0][0] = (-6*alpha1*alpha2)/5.
            A_last_row[0][1] = 0
            A_last_row[0][2] = 0
            A_last_row[0][3] = (2*alpha1)/5.
        
        elif order == 3:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0

            A_last_row[0][0] = (-2*(1287*alpha2*alpha2 + 65*alpha3*(44*alpha1 + \
            9*alpha3)))/5005.
            A_last_row[0][1] = 0
            A_last_row[0][2] = (-2*alpha3)/7.
            A_last_row[0][3] = (6*alpha2)/35.
            A_last_row[0][4] = (3*alpha1)/7. + (3*alpha3)/11.
        
        elif order == 4:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0
                alpha4 = 0

            A_last_row[0][0] = (-2*(390*alpha2*alpha3 + (455*alpha1 + \
            180*alpha3)*alpha4))/819.
            A_last_row[0][1] = 0
            A_last_row[0][2] = (-5*alpha4)/9.
            A_last_row[0][3] = 0
            A_last_row[0][4] = (5*(alpha2 + alpha4))/21.
            A_last_row[0][5] = (4*alpha1)/9. + (3*alpha3)/13.
        
        elif order == 5:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]

            if self.hyperbolic:
                alpha2 = 0
                alpha3 = 0
                alpha4 = 0
                alpha5 = 0

            A_last_row[0][0] = (-100*alpha3*alpha3)/231. - (10*alpha2*alpha4)/11. \
            - (14*alpha3*alpha5)/33. - (4*(85*alpha4*alpha4 + 459*alpha1*alpha5 + \
            60*alpha5*alpha5))/1683.
            A_last_row[0][1] = 0
            A_last_row[0][2] = (-9*alpha5)/11.
            A_last_row[0][3] = (-5*alpha4)/33.
            A_last_row[0][4] = (25*alpha3 + 49*alpha5)/231.
            A_last_row[0][5] = (3*alpha2)/11. + (19*alpha4)/99.
            A_last_row[0][6] = (255*alpha1 + 119*alpha3 + 104*alpha5)/561.
        
        else:
            print("This order is not implemented for the system matrix last row of the SWME1D-PDE.")

        return A_last_row
    
    def compute_source_term(self,
                            order: int,
                            values: np.array,
                            delta_t: float,
                            **kwargs) -> np.array:
        
        viscosity   = kwargs["viscosity"]   if "viscosity"   in kwargs else self.viscosity
        slip_length = kwargs["slip_length"] if "slip_length" in kwargs else self.slip_length
        g           = kwargs["g"]           if "g"           in kwargs else 1

        if self.linear_source:
            S = self._compute_source_matrix_inverse(order,values,delta_t,**kwargs)
        
        else:
            S = np.zeros(order+2) 
            h = values[0]
            um = values[1]/values[0]
            
            if order == 0:
                S[0] = 0
                S[1] = -viscosity/slip_length*um
            
            elif order == 1:
                alpha1 = values[2]/values[0]

                S[0] = 0
                S[1] = -viscosity/slip_length*(um + alpha1)
                S[2] = -3*viscosity/slip_length*(um + (1 + 4*slip_length/h)*alpha1)

            elif order == 2:
                alpha1 = values[2]/values[0]
                alpha2 = values[3]/values[0]

                S[0] = 0
                S[1] = -viscosity/slip_length*(um + alpha1 + alpha2)
                S[2] = -3*viscosity/slip_length*(um + (1 + 4*slip_length/h)*alpha1 + alpha2)
                S[3] = -5*viscosity/slip_length*(um + alpha1 + (1 + 12*slip_length/h)*alpha2)
            
            elif order == 3:
                alpha1 = values[2]/values[0]
                alpha2 = values[3]/values[0]
                alpha3 = values[4]/values[0]

                S[0] = 0
                S[1] = -viscosity/slip_length*(um + alpha1 + alpha2 + alpha3)
                S[2] = -3*viscosity/slip_length*((h + 4*slip_length)*alpha1 + h*(um + alpha2) + (h + 4*slip_length)*alpha3)/h
                S[3] = -5*viscosity/slip_length*(um + alpha1 + (1 + 12*slip_length/h)*alpha2 + alpha3)
                S[4] = -7*viscosity/slip_length*((h + 4*slip_length)*alpha1 + h*(um + alpha2) + (h + 24*slip_length)*alpha3)/h
            
            elif order == 4:
                alpha1 = values[2]/values[0]
                alpha2 = values[3]/values[0]
                alpha3 = values[4]/values[0]
                alpha4 = values[5]/values[0]

                S[0] = 0
                S[1] = -((viscosity*(um + alpha1 + alpha2 + alpha3 + \
                alpha4))/slip_length)
                S[2] = (-3*viscosity*(um + alpha2 + alpha3 + ((h + \
                4*slip_length)*alpha1 + 4*slip_length*alpha3)/h + \
                alpha4))/slip_length
                S[3] = (-5*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
                12*slip_length)*alpha2 + \
                12*slip_length*alpha4)/h))/slip_length
                S[4] = (-7*viscosity*(um + alpha2 + alpha3 + ((h + \
                4*slip_length)*alpha1 + 24*slip_length*alpha3)/h + \
                alpha4))/slip_length
                S[5] = (-9*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
                12*slip_length)*alpha2 + \
                40*slip_length*alpha4)/h))/slip_length

            elif order == 5:
                alpha1 = values[2]/values[0]
                alpha2 = values[3]/values[0]
                alpha3 = values[4]/values[0]
                alpha4 = values[5]/values[0]
                alpha5 = values[6]/values[0]

                S[0] = 0
                S[1] = -((viscosity*(um + alpha1 + alpha2 + alpha3 + alpha4 + \
                alpha5))/slip_length)
                S[2] = (-3*viscosity*(h*um + (h + 4*slip_length)*alpha1 + \
                h*alpha2 + (h + 4*slip_length)*alpha3 + h*alpha4 + (h + \
                4*slip_length)*alpha5))/(h*slip_length)
                S[3] = (-5*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
                12*slip_length)*alpha2 + 12*slip_length*alpha4)/h + \
                alpha5))/slip_length
                S[4] = (-7*viscosity*(h*um + (h + 4*slip_length)*alpha1 + \
                h*alpha2 + (h + 24*slip_length)*alpha3 + h*alpha4 + (h + \
                24*slip_length)*alpha5))/(h*slip_length)
                S[5] = (-9*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
                12*slip_length)*alpha2 + 40*slip_length*alpha4)/h + \
                alpha5))/slip_length
                S[6] = (-11*viscosity*(h*um + (h + 4*slip_length)*alpha1 + \
                h*alpha2 + (h + 24*slip_length)*alpha3 + h*alpha4 + (h + \
                60*slip_length)*alpha5))/(h*slip_length)

            elif order == 6:
                alpha1 = values[2]/values[0]
                alpha2 = values[3]/values[0]
                alpha3 = values[4]/values[0]
                alpha4 = values[5]/values[0]
                alpha5 = values[6]/values[0]
                alpha6 = values[7]/values[0]

                S[0] = 0
                S[1] = -((viscosity*(um + alpha1 + alpha2 + alpha3 + alpha4 + \
                alpha5 + alpha6))/slip_length)
                S[2] = (-3*viscosity*(um + alpha2 + alpha3 + alpha4 + alpha5 + \
                ((h + 4*slip_length)*alpha1 + 4*slip_length*(alpha3 + \
                alpha5))/h + alpha6))/slip_length
                S[3] = (-5*viscosity*(um + alpha1 + alpha3 + alpha4 + alpha5 + \
                alpha6 + ((h + 12*slip_length)*alpha2 + \
                12*slip_length*(alpha4 + alpha6))/h))/slip_length
                S[4] = (-7*viscosity*(um + alpha2 + alpha3 + alpha4 + alpha5 + \
                ((h + 4*slip_length)*alpha1 + 24*slip_length*(alpha3 + \
                alpha5))/h + alpha6))/slip_length
                S[5] = (-9*viscosity*(um + alpha1 + alpha3 + alpha4 + alpha5 + \
                alpha6 + ((h + 12*slip_length)*alpha2 + \
                40*slip_length*(alpha4 + alpha6))/h))/slip_length
                S[6] = (-11*viscosity*(um + alpha2 + alpha3 + alpha4 + alpha5 + \
                ((h + 4*slip_length)*alpha1 + 12*slip_length*(2*alpha3 + \
                5*alpha5))/h + alpha6))/slip_length
                S[7] = (-13*viscosity*(um + alpha1 + alpha3 + alpha4 + alpha5 + \
                alpha6 + ((h + 12*slip_length)*alpha2 + \
                40*slip_length*alpha4 + \
                84*slip_length*alpha6)/h))/slip_length
            
            else:
                print("This order is not implemented for the source term of the SWME1D-PDE.")

        return S
    
    def _compute_source_matrix_inverse(self,
                                      order: int,
                                      values: np.array,
                                      delta_t: float,
                                      **kwargs) -> np.array:
        
        viscosity   = kwargs["viscosity"]   if "viscosity"   in kwargs else self.viscosity
        slip_length = kwargs["slip_length"] if "slip_length" in kwargs else self.slip_length
        g           = kwargs["g"]           if "g"           in kwargs else 1

        S_inv = np.zeros((order+2,order+2)) 
        h = values[0]
        
        if order == 0:
            S_inv[0][0] = 1
            S_inv[0][1] = 0
            S_inv[1][0] = 0
            S_inv[1][1] = (h*slip_length)/(h*slip_length + \
            delta_t*viscosity)
            
        elif order == 1:
            S_inv[0][0] = 1
            S_inv[0][1] = 0
            S_inv[0][2] = 0
            S_inv[1][0] = 0
            S_inv[1][1] = (h**3*slip_length + 3*h*delta_t*(h + \
            4*slip_length)*viscosity)/(h**3*slip_length + \
            4*h*delta_t*(h + 3*slip_length)*viscosity + \
            12*delta_t**2*viscosity**2)
            S_inv[1][2] = -((h**2*delta_t*viscosity)/(h**3*slip_length \
            + 4*h*delta_t*(h + 3*slip_length)*viscosity + \
            12*delta_t**2*viscosity**2))
            S_inv[2][0] = 0
            S_inv[2][1] = (-3*h**2*delta_t*viscosity)/(h**3*slip_length \
            + 4*h*delta_t*(h + 3*slip_length)*viscosity + \
            12*delta_t**2*viscosity**2)
            S_inv[2][2] = (h**2*(h*slip_length + \
            delta_t*viscosity))/(h**3*slip_length + 4*h*delta_t*(h + \
            3*slip_length)*viscosity + \
            12*delta_t**2*viscosity**2)
        
        elif order == 2:
            S_inv[0][0] = 1
            S_inv[0][1] = 0
            S_inv[0][2] = 0
            S_inv[0][3] = 0
            S_inv[1][0] = 0
            S_inv[1][1] = (h**5*slip_length + 8*h**3*delta_t*(h + \
            9*slip_length)*viscosity + 240*h*delta_t**2*(h + \
            3*slip_length)*viscosity**2)/(h**5*slip_length + \
            9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3)
            S_inv[1][2] = -((h**2*delta_t*viscosity*(h**2 + \
            60*delta_t*viscosity))/(h**5*slip_length + \
            9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3))
            S_inv[1][3] = -((h**2*delta_t*viscosity*(h**2 + \
            12*delta_t*viscosity))/(h**5*slip_length + \
            9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3))
            S_inv[2][0] = 0
            S_inv[2][1] = (-3*h**2*delta_t*viscosity*(h**2 + \
            60*delta_t*viscosity))/(h**5*slip_length + \
            9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3)
            S_inv[2][2] = (h**2*(h**3*slip_length + 6*h*delta_t*(h + \
            10*slip_length)*viscosity + \
            60*delta_t**2*viscosity**2))/(h**5*slip_length + \
            9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3)
            S_inv[2][3] = (-3*h**4*delta_t*viscosity)/(h**5*slip_length \
            + 9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3)
            S_inv[3][0] = 0
            S_inv[3][1] = (-5*h**2*delta_t*viscosity*(h**2 + \
            12*delta_t*viscosity))/(h**5*slip_length + \
            9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3)
            S_inv[3][2] = (-5*h**4*delta_t*viscosity)/(h**5*slip_length \
            + 9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3)
            S_inv[3][3] = (h**2*(h**3*slip_length + 4*h*delta_t*(h + \
            3*slip_length)*viscosity + \
            12*delta_t**2*viscosity**2))/(h**5*slip_length + \
            9*h**3*delta_t*(h + 8*slip_length)*viscosity + \
            24*h*delta_t**2*(13*h + 30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3)

        elif order == 3:
            S_inv[0][0] = 1
            S_inv[0][1] = 0
            S_inv[0][2] = 0
            S_inv[0][3] = 0
            S_inv[0][4] = 0
            S_inv[1][0] = 0
            S_inv[1][1] = (h**7*slip_length + 15*h**5*delta_t*(h + \
            16*slip_length)*viscosity + 960*h**3*delta_t**2*(2*h + \
            13*slip_length)*viscosity**2 + 33600*h*delta_t**3*(h + \
            3*slip_length)*viscosity**3)/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[1][2] = -((h**2*delta_t*viscosity*(h**2 + \
            60*delta_t*viscosity)*(h**2 + \
            140*delta_t*viscosity))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4))
            S_inv[1][3] = -((h**2*delta_t*viscosity*(h**4 + \
            180*h**2*delta_t*viscosity + \
            1680*delta_t**2*viscosity**2))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4))
            S_inv[1][4] = -((h**4*delta_t*viscosity*(h**2 + \
            60*delta_t*viscosity))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4))
            S_inv[2][0] = 0
            S_inv[2][1] = (-3*h**2*delta_t*viscosity*(h**2 + \
            60*delta_t*viscosity)*(h**2 + \
            140*delta_t*viscosity))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[2][2] = (h**2*(h**5*slip_length + h**3*delta_t*(13*h + \
            228*slip_length)*viscosity + 48*h*delta_t**2*(31*h + \
            210*slip_length)*viscosity**2 + \
            10080*delta_t**3*viscosity**3))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[2][3] = (-3*h**4*delta_t*viscosity*(h**2 + \
            140*delta_t*viscosity))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[2][4] = (-3*h**2*delta_t*viscosity*(h**3*(h + \
            4*slip_length) + 12*h*delta_t*(7*h + \
            20*slip_length)*viscosity + \
            240*delta_t**2*viscosity**2))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[3][0] = 0
            S_inv[3][1] = (-5*h**2*delta_t*viscosity*(h**4 + \
            180*h**2*delta_t*viscosity + \
            1680*delta_t**2*viscosity**2))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[3][2] = (-5*h**4*delta_t*viscosity*(h**2 + \
            140*delta_t*viscosity))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[3][3] = (h**2*(h**5*slip_length + h**3*delta_t*(11*h + \
            180*slip_length)*viscosity + 120*h*delta_t**2*(5*h + \
            14*slip_length)*viscosity**2 + \
            1680*delta_t**3*viscosity**3))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[3][4] = (-5*h**6*delta_t*viscosity)/(h**7*slip_length \
            + 16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[4][0] = 0
            S_inv[4][1] = (-7*h**4*delta_t*viscosity*(h**2 + \
            60*delta_t*viscosity))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[4][2] = (-7*h**2*delta_t*viscosity*(h**3*(h + \
            4*slip_length) + 12*h*delta_t*(7*h + \
            20*slip_length)*viscosity + \
            240*delta_t**2*viscosity**2))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[4][3] = (-7*h**6*delta_t*viscosity)/(h**7*slip_length \
            + 16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)
            S_inv[4][4] = (h**2*(h**5*slip_length + 9*h**3*delta_t*(h + \
            8*slip_length)*viscosity + 24*h*delta_t**2*(13*h + \
            30*slip_length)*viscosity**2 + \
            720*delta_t**3*viscosity**3))/(h**7*slip_length + \
            16*h**5*delta_t*(h + 15*slip_length)*viscosity + \
            240*h**3*delta_t**2*(9*h + 52*slip_length)*viscosity**2 + \
            2880*h*delta_t**3*(16*h + 35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4)

        elif order == 4:
            S_inv[0][0] = 1
            S_inv[0][1] = 0
            S_inv[0][2] = 0
            S_inv[0][3] = 0
            S_inv[0][4] = 0
            S_inv[0][5] = 0
            S_inv[1][0] = 0
            S_inv[1][1] = (h**9*slip_length + 24*h**7*delta_t*(h + \
            25*slip_length)*viscosity + 8400*h**5*delta_t**2*(h + \
            11*slip_length)*viscosity**2 + 13440*h**3*delta_t**3*(43*h \
            + 255*slip_length)*viscosity**3 + 8467200*h*delta_t**4*(h + \
            3*slip_length)*viscosity**4)/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[1][2] = -((h**2*delta_t*viscosity*(h**2 + \
            140*delta_t*viscosity)*(h**4 + 420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5))
            S_inv[1][3] = -((h**2*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity)*(h**4 + 180*h**2*delta_t*viscosity + \
            1680*delta_t**2*viscosity**2))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5))
            S_inv[1][4] = -((h**4*delta_t*viscosity*(h**4 + \
            420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5))
            S_inv[1][5] = -((h**4*delta_t*viscosity*(h**4 + \
            180*h**2*delta_t*viscosity + \
            1680*delta_t**2*viscosity**2))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5))
            S_inv[2][0] = 0
            S_inv[2][1] = (-3*h**2*delta_t*viscosity*(h**2 + \
            140*delta_t*viscosity)*(h**4 + 420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[2][2] = (h**2*(h**7*slip_length + 2*h**5*delta_t*(11*h + \
            294*slip_length)*viscosity + 7140*h**3*delta_t**2*(h + \
            12*slip_length)*viscosity**2 + 40320*h*delta_t**3*(10*h + \
            63*slip_length)*viscosity**3 + \
            2540160*delta_t**4*viscosity**4))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[2][3] = (-3*h**4*delta_t*viscosity*(h**2 + \
            140*delta_t*viscosity)*(h**2 + \
            252*delta_t*viscosity))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[2][4] = (3*h**2*delta_t*viscosity*(-(h**5*(h + \
            4*slip_length)) - 240*h**3*delta_t*(2*h + \
            7*slip_length)*viscosity - 1680*h*delta_t**2*(13*h + \
            36*slip_length)*viscosity**2 - \
            60480*delta_t**3*viscosity**3))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[2][5] = (-3*h**6*delta_t*viscosity*(h**2 + \
            140*delta_t*viscosity))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[3][0] = 0
            S_inv[3][1] = (-5*h**2*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity)*(h**4 + 180*h**2*delta_t*viscosity + \
            1680*delta_t**2*viscosity**2))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[3][2] = (-5*h**4*delta_t*viscosity*(h**2 + \
            140*delta_t*viscosity)*(h**2 + \
            252*delta_t*viscosity))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[3][3] = (h**2*(h**7*slip_length + 20*h**5*delta_t*(h + \
            27*slip_length)*viscosity + 60*h**3*delta_t**2*(103*h + \
            1108*slip_length)*viscosity**2 + 2400*h*delta_t**3*(97*h + \
            252*slip_length)*viscosity**3 + \
            604800*delta_t**4*viscosity**4))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[3][4] = (-5*h**6*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[3][5] = (-5*h**2*delta_t*viscosity*(h**5*(h + \
            12*slip_length) + 24*h**3*delta_t*(13*h + \
            90*slip_length)*viscosity + 240*h*delta_t**2*(37*h + \
            84*slip_length)*viscosity**2 + \
            20160*delta_t**3*viscosity**3))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[4][0] = 0
            S_inv[4][1] = (-7*h**4*delta_t*viscosity*(h**4 + \
            420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[4][2] = (7*h**2*delta_t*viscosity*(-(h**5*(h + \
            4*slip_length)) - 240*h**3*delta_t*(2*h + \
            7*slip_length)*viscosity - 1680*h*delta_t**2*(13*h + \
            36*slip_length)*viscosity**2 - \
            60480*delta_t**3*viscosity**3))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[4][3] = (-7*h**6*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[4][4] = (h**2*(h**7*slip_length + 18*h**5*delta_t*(h + \
            24*slip_length)*viscosity + 240*h**3*delta_t**2*(13*h + \
            84*slip_length)*viscosity**2 + 20160*h*delta_t**3*(4*h + \
            9*slip_length)*viscosity**3 + \
            181440*delta_t**4*viscosity**4))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[4][5] = (-7*h**8*delta_t*viscosity)/(h**9*slip_length \
            + 25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[5][0] = 0
            S_inv[5][1] = (-9*h**4*delta_t*viscosity*(h**4 + \
            180*h**2*delta_t*viscosity + \
            1680*delta_t**2*viscosity**2))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[5][2] = (-9*h**6*delta_t*viscosity*(h**2 + \
            140*delta_t*viscosity))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[5][3] = (-9*h**2*delta_t*viscosity*(h**5*(h + \
            12*slip_length) + 24*h**3*delta_t*(13*h + \
            90*slip_length)*viscosity + 240*h*delta_t**2*(37*h + \
            84*slip_length)*viscosity**2 + \
            20160*delta_t**3*viscosity**3))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[5][4] = (-9*h**8*delta_t*viscosity)/(h**9*slip_length \
            + 25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)
            S_inv[5][5] = (h**2*(h**7*slip_length + 16*h**5*delta_t*(h + \
            15*slip_length)*viscosity + 240*h**3*delta_t**2*(9*h + \
            52*slip_length)*viscosity**2 + 2880*h*delta_t**3*(16*h + \
            35*slip_length)*viscosity**3 + \
            100800*delta_t**4*viscosity**4))/(h**9*slip_length + \
            25*h**7*delta_t*(h + 24*slip_length)*viscosity + \
            600*h**5*delta_t**2*(15*h + 154*slip_length)*viscosity**2 + \
            5040*h**3*delta_t**3*(133*h + 680*slip_length)*viscosity**3 \
            + 201600*h*delta_t**4*(59*h + 126*slip_length)*viscosity**4 \
            + 25401600*delta_t**5*viscosity**5)

        elif order == 5:
            S_inv[0][0] = 1
            S_inv[0][1] = 0
            S_inv[0][2] = 0
            S_inv[0][3] = 0
            S_inv[0][4] = 0
            S_inv[0][5] = 0
            S_inv[0][6] = 0
            S_inv[1][0] = 0
            S_inv[1][1] = (h**11*slip_length + 35*h**9*delta_t*(h + \
            36*slip_length)*viscosity + 13440*h**7*delta_t**2*(2*h + \
            33*slip_length)*viscosity**2 + 120960*h**5*delta_t**3*(39*h \
            + 373*slip_length)*viscosity**3 + \
            6773760*h**3*delta_t**4*(37*h + \
            210*slip_length)*viscosity**4 + 3353011200*h*delta_t**5*(h \
            + 3*slip_length)*viscosity**5)/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[1][2] = -((h**2*delta_t*viscosity*(h**4 + \
            420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2)*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6))
            S_inv[1][3] = -((h**2*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity)*(h**6 + 840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6))
            S_inv[1][4] = -((h**4*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**4 + 420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6))
            S_inv[1][5] = -((h**4*delta_t*viscosity*(h**6 + \
            840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6))
            S_inv[1][6] = -((h**6*delta_t*viscosity*(h**4 + \
            420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6))
            S_inv[2][0] = 0
            S_inv[2][1] = (-3*h**2*delta_t*viscosity*(h**4 + \
            420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2)*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[2][2] = (h**2*(h**9*slip_length + 3*h**7*delta_t*(11*h + \
            416*slip_length)*viscosity + 48*h**5*delta_t**2*(509*h + \
            8946*slip_length)*viscosity**2 + \
            30240*h**3*delta_t**3*(127*h + \
            1338*slip_length)*viscosity**3 + 725760*h*delta_t**4*(229*h \
            + 1386*slip_length)*viscosity**4 + \
            1005903360*delta_t**5*viscosity**5))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[2][3] = (-3*h**4*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity)*(h**4 + 756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[2][4] = (-3*h**2*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**5*(h + 4*slip_length) + \
            240*h**3*delta_t*(2*h + 7*slip_length)*viscosity + \
            1680*h*delta_t**2*(13*h + 36*slip_length)*viscosity**2 + \
            60480*delta_t**3*viscosity**3))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[2][5] = (-3*h**6*delta_t*viscosity*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[2][6] = (3*h**4*delta_t*viscosity*(-(h**5*(h + \
            4*slip_length)) - 240*h**3*delta_t*(2*h + \
            7*slip_length)*viscosity - 1680*h*delta_t**2*(13*h + \
            36*slip_length)*viscosity**2 - \
            60480*delta_t**3*viscosity**3))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[3][0] = 0
            S_inv[3][1] = (-5*h**2*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity)*(h**6 + 840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[3][2] = (-5*h**4*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity)*(h**4 + 756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[3][3] = (h**2*(h**9*slip_length + h**7*delta_t*(31*h + \
            1200*slip_length)*viscosity + 240*h**5*delta_t**2*(89*h + \
            1575*slip_length)*viscosity**2 + \
            15120*h**3*delta_t**3*(201*h + \
            1844*slip_length)*viscosity**3 + 604800*h*delta_t**4*(155*h \
            + 396*slip_length)*viscosity**4 + \
            239500800*delta_t**5*viscosity**5))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[3][4] = (-5*h**6*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity)*(h**2 + \
            396*delta_t*viscosity))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[3][5] = (5*h**2*delta_t*viscosity*(-(h**7*(h + \
            12*slip_length)) - 48*h**5*delta_t*(23*h + \
            210*slip_length)*viscosity - 5040*h**3*delta_t**2*(29*h + \
            180*slip_length)*viscosity**2 - 60480*h*delta_t**3*(59*h + \
            132*slip_length)*viscosity**3 - \
            7983360*delta_t**4*viscosity**4))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[3][6] = (-5*h**8*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[4][0] = 0
            S_inv[4][1] = (-7*h**4*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**4 + 420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[4][2] = (-7*h**2*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**5*(h + 4*slip_length) + \
            240*h**3*delta_t*(2*h + 7*slip_length)*viscosity + \
            1680*h*delta_t**2*(13*h + 36*slip_length)*viscosity**2 + \
            60480*delta_t**3*viscosity**3))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[4][3] = (-7*h**6*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity)*(h**2 + \
            396*delta_t*viscosity))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[4][4] = (h**2*(h**9*slip_length + h**7*delta_t*(29*h + \
            1092*slip_length)*viscosity + 336*h**5*delta_t**2*(58*h + \
            907*slip_length)*viscosity**2 + 5040*h**3*delta_t**3*(445*h \
            + 2632*slip_length)*viscosity**3 + \
            282240*h*delta_t**4*(179*h + 396*slip_length)*viscosity**4 \
            + 111767040*delta_t**5*viscosity**5))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[4][5] = (-7*h**8*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[4][6] = (7*h**2*delta_t*viscosity*(-(h**7*(h + \
            24*slip_length)) - 120*h**5*delta_t*(7*h + \
            86*slip_length)*viscosity - 720*h**3*delta_t**2*(117*h + \
            644*slip_length)*viscosity**2 - 20160*h*delta_t**3*(83*h + \
            180*slip_length)*viscosity**3 - \
            3628800*delta_t**4*viscosity**4))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[5][0] = 0
            S_inv[5][1] = (-9*h**4*delta_t*viscosity*(h**6 + \
            840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[5][2] = (-9*h**6*delta_t*viscosity*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[5][3] = (9*h**2*delta_t*viscosity*(-(h**7*(h + \
            12*slip_length)) - 48*h**5*delta_t*(23*h + \
            210*slip_length)*viscosity - 5040*h**3*delta_t**2*(29*h + \
            180*slip_length)*viscosity**2 - 60480*h*delta_t**3*(59*h + \
            132*slip_length)*viscosity**3 - \
            7983360*delta_t**4*viscosity**4))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[5][4] = (-9*h**8*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[5][5] = (h**2*(h**9*slip_length + 9*h**7*delta_t*(3*h + \
            100*slip_length)*viscosity + 600*h**5*delta_t**2*(19*h + \
            210*slip_length)*viscosity**2 + 5040*h**3*delta_t**3*(193*h \
            + 1032*slip_length)*viscosity**3 + \
            362880*h*delta_t**4*(51*h + 110*slip_length)*viscosity**4 + \
            39916800*delta_t**5*viscosity**5))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[5][6] = \
            (-9*h**10*delta_t*viscosity)/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[6][0] = 0
            S_inv[6][1] = (-11*h**6*delta_t*viscosity*(h**4 + \
            420*h**2*delta_t*viscosity + \
            15120*delta_t**2*viscosity**2))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[6][2] = (11*h**4*delta_t*viscosity*(-(h**5*(h + \
            4*slip_length)) - 240*h**3*delta_t*(2*h + \
            7*slip_length)*viscosity - 1680*h*delta_t**2*(13*h + \
            36*slip_length)*viscosity**2 - \
            60480*delta_t**3*viscosity**3))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[6][3] = (-11*h**8*delta_t*viscosity*(h**2 + \
            252*delta_t*viscosity))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[6][4] = (11*h**2*delta_t*viscosity*(-(h**7*(h + \
            24*slip_length)) - 120*h**5*delta_t*(7*h + \
            86*slip_length)*viscosity - 720*h**3*delta_t**2*(117*h + \
            644*slip_length)*viscosity**2 - 20160*h*delta_t**3*(83*h + \
            180*slip_length)*viscosity**3 - \
            3628800*delta_t**4*viscosity**4))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[6][5] = \
            (-11*h**10*delta_t*viscosity)/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)
            S_inv[6][6] = (h**2*(h**9*slip_length + 25*h**7*delta_t*(h + \
            24*slip_length)*viscosity + 600*h**5*delta_t**2*(15*h + \
            154*slip_length)*viscosity**2 + 5040*h**3*delta_t**3*(133*h \
            + 680*slip_length)*viscosity**3 + 201600*h*delta_t**4*(59*h \
            + 126*slip_length)*viscosity**4 + \
            25401600*delta_t**5*viscosity**5))/(h**11*slip_length + \
            36*h**9*delta_t*(h + 35*slip_length)*viscosity + \
            420*h**7*delta_t**2*(67*h + 1056*slip_length)*viscosity**2 \
            + 40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6)

        elif order == 6:
            S_inv[0][0] = 1
            S_inv[0][1] = 0
            S_inv[0][2] = 0
            S_inv[0][3] = 0
            S_inv[0][4] = 0
            S_inv[0][5] = 0
            S_inv[0][6] = 0
            S_inv[0][7] = 0
            S_inv[1][0] = 0
            S_inv[1][1] = (h**13*slip_length + 48*h**11*delta_t*(h + \
            49*slip_length)*viscosity + 70560*h**9*delta_t**2*(h + \
            23*slip_length)*viscosity**2 + 282240*h**7*delta_t**3*(91*h \
            + 1263*slip_length)*viscosity**3 + \
            2661120*h**5*delta_t**4*(1237*h + \
            10983*slip_length)*viscosity**4 + \
            1341204480*h**3*delta_t**5*(113*h + \
            625*slip_length)*viscosity**5 + \
            1917922406400*h*delta_t**6*(h + \
            3*slip_length)*viscosity**6)/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[1][2] = -((h**2*delta_t*viscosity*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2)*(h**6 + \
            1512*h**4*delta_t*viscosity + \
            277200*h**2*delta_t**2*viscosity**2 + \
            8648640*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7))
            S_inv[1][3] = -((h**2*delta_t*viscosity*(h**4 + \
            1188*h**2*delta_t*viscosity + \
            144144*delta_t**2*viscosity**2)*(h**6 + \
            840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7))
            S_inv[1][4] = -((h**4*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**6 + 1512*h**4*delta_t*viscosity \
            + 277200*h**2*delta_t**2*viscosity**2 + \
            8648640*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7))
            S_inv[1][5] = -((h**4*delta_t*viscosity*(h**2 + \
            572*delta_t*viscosity)*(h**6 + 840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7))
            S_inv[1][6] = -((h**6*delta_t*viscosity*(h**6 + \
            1512*h**4*delta_t*viscosity + \
            277200*h**2*delta_t**2*viscosity**2 + \
            8648640*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7))
            S_inv[1][7] = -((h**6*delta_t*viscosity*(h**6 + \
            840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7))
            S_inv[2][0] = 0
            S_inv[2][1] = (-3*h**2*delta_t*viscosity*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2)*(h**6 + \
            1512*h**4*delta_t*viscosity + \
            277200*h**2*delta_t**2*viscosity**2 + \
            8648640*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[2][2] = (h**2*(h**11*slip_length + 2*h**9*delta_t*(23*h + \
            1170*slip_length)*viscosity + 252*h**7*delta_t**2*(261*h + \
            6332*slip_length)*viscosity**2 + \
            24192*h**5*delta_t**3*(929*h + \
            14003*slip_length)*viscosity**3 + \
            2661120*h**3*delta_t**4*(976*h + \
            9621*slip_length)*viscosity**4 + \
            191600640*h*delta_t**5*(509*h + \
            3003*slip_length)*viscosity**5 + \
            575376721920*delta_t**6*viscosity**6))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[2][3] = (-3*h**4*delta_t*viscosity*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2)*(h**4 + \
            1188*h**2*delta_t*viscosity + \
            144144*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[2][4] = (-3*h**2*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**7*(h + 4*slip_length) + \
            56*h**5*delta_t*(29*h + 108*slip_length)*viscosity + \
            25200*h**3*delta_t**2*(13*h + 44*slip_length)*viscosity**2 \
            + 665280*h*delta_t**3*(19*h + 52*slip_length)*viscosity**3 \
            + 34594560*delta_t**4*viscosity**4))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[2][5] = (-3*h**6*delta_t*viscosity*(h**2 + \
            572*delta_t*viscosity)*(h**4 + 756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[2][6] = (3*h**4*delta_t*viscosity*(-(h**7*(h + \
            4*slip_length)) - 56*h**5*delta_t*(29*h + \
            108*slip_length)*viscosity - 25200*h**3*delta_t**2*(13*h + \
            44*slip_length)*viscosity**2 - 665280*h*delta_t**3*(19*h + \
            52*slip_length)*viscosity**3 - \
            34594560*delta_t**4*viscosity**4))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[2][7] = (-3*h**8*delta_t*viscosity*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[3][0] = 0
            S_inv[3][1] = (-5*h**2*delta_t*viscosity*(h**4 + \
            1188*h**2*delta_t*viscosity + \
            144144*delta_t**2*viscosity**2)*(h**6 + \
            840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[3][2] = (-5*h**4*delta_t*viscosity*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2)*(h**4 + \
            1188*h**2*delta_t*viscosity + \
            144144*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[3][3] = (h**2*(h**11*slip_length + 4*h**9*delta_t*(11*h + \
            573*slip_length)*viscosity + 108*h**7*delta_t**2*(569*h + \
            13900*slip_length)*viscosity**2 + \
            12960*h**5*delta_t**3*(1489*h + \
            21868*slip_length)*viscosity**3 + \
            2661120*h**3*delta_t**4*(739*h + \
            6213*slip_length)*viscosity**4 + \
            79833600*h*delta_t**5*(679*h + \
            1716*slip_length)*viscosity**5 + \
            136994457600*delta_t**6*viscosity**6))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[3][4] = (-5*h**6*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**4 + 1188*h**2*delta_t*viscosity \
            + 144144*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[3][5] = (-5*h**2*delta_t*viscosity*(h**2 + \
            572*delta_t*viscosity)*(h**7*(h + 12*slip_length) + \
            48*h**5*delta_t*(23*h + 210*slip_length)*viscosity + \
            5040*h**3*delta_t**2*(29*h + 180*slip_length)*viscosity**2 \
            + 60480*h*delta_t**3*(59*h + 132*slip_length)*viscosity**3 \
            + 7983360*delta_t**4*viscosity**4))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[3][6] = (-5*h**8*delta_t*viscosity*(h**4 + \
            1188*h**2*delta_t*viscosity + \
            144144*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[3][7] = (5*h**4*delta_t*viscosity*(-(h**7*(h + \
            12*slip_length)) - 48*h**5*delta_t*(23*h + \
            210*slip_length)*viscosity - 5040*h**3*delta_t**2*(29*h + \
            180*slip_length)*viscosity**2 - 60480*h*delta_t**3*(59*h + \
            132*slip_length)*viscosity**3 - \
            7983360*delta_t**4*viscosity**4))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[4][0] = 0
            S_inv[4][1] = (-7*h**4*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**6 + 1512*h**4*delta_t*viscosity \
            + 277200*h**2*delta_t**2*viscosity**2 + \
            8648640*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[4][2] = (-7*h**2*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**7*(h + 4*slip_length) + \
            56*h**5*delta_t*(29*h + 108*slip_length)*viscosity + \
            25200*h**3*delta_t**2*(13*h + 44*slip_length)*viscosity**2 \
            + 665280*h*delta_t**3*(19*h + 52*slip_length)*viscosity**3 \
            + 34594560*delta_t**4*viscosity**4))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[4][3] = (-7*h**6*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**4 + 1188*h**2*delta_t*viscosity \
            + 144144*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[4][4] = (h**2*(h**11*slip_length + 42*h**9*delta_t*(h + \
            52*slip_length)*viscosity + 336*h**7*delta_t**2*(162*h + \
            3871*slip_length)*viscosity**2 + \
            1344*h**5*delta_t**3*(12163*h + \
            153351*slip_length)*viscosity**3 + \
            665280*h**3*delta_t**4*(2113*h + \
            11816*slip_length)*viscosity**4 + \
            111767040*h*delta_t**5*(261*h + \
            572*slip_length)*viscosity**5 + \
            63930746880*delta_t**6*viscosity**6))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[4][5] = (-7*h**8*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**2 + \
            572*delta_t*viscosity))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[4][6] = (7*h**2*delta_t*viscosity*(-(h**9*(h + \
            24*slip_length)) - 12*h**7*delta_t*(187*h + \
            3044*slip_length)*viscosity - 1680*h**5*delta_t**2*(403*h + \
            4176*slip_length)*viscosity**2 - \
            241920*h**3*delta_t**3*(216*h + \
            1133*slip_length)*viscosity**3 - \
            7983360*h*delta_t**4*(121*h + 260*slip_length)*viscosity**4 \
            - 2075673600*delta_t**5*viscosity**5))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[4][7] = (-7*h**10*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[5][0] = 0
            S_inv[5][1] = (-9*h**4*delta_t*viscosity*(h**2 + \
            572*delta_t*viscosity)*(h**6 + 840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[5][2] = (-9*h**6*delta_t*viscosity*(h**2 + \
            572*delta_t*viscosity)*(h**4 + 756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[5][3] = (-9*h**2*delta_t*viscosity*(h**2 + \
            572*delta_t*viscosity)*(h**7*(h + 12*slip_length) + \
            48*h**5*delta_t*(23*h + 210*slip_length)*viscosity + \
            5040*h**3*delta_t**2*(29*h + 180*slip_length)*viscosity**2 \
            + 60480*h*delta_t**3*(59*h + 132*slip_length)*viscosity**3 \
            + 7983360*delta_t**4*viscosity**4))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[5][4] = (-9*h**8*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity)*(h**2 + \
            572*delta_t*viscosity))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[5][5] = (h**2*(h**11*slip_length + 8*h**9*delta_t*(5*h + \
            249*slip_length)*viscosity + 48*h**7*delta_t**2*(1063*h + \
            22905*slip_length)*viscosity**2 + \
            4320*h**5*delta_t**3*(3135*h + \
            31234*slip_length)*viscosity**3 + \
            60480*h**3*delta_t**4*(16127*h + \
            82872*slip_length)*viscosity**4 + \
            13063680*h*delta_t**5*(1337*h + \
            2860*slip_length)*viscosity**5 + \
            37362124800*delta_t**6*viscosity**6))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[5][6] = (-9*h**10*delta_t*viscosity*(h**2 + \
            572*delta_t*viscosity))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[5][7] = (9*h**2*delta_t*viscosity*(-(h**9*(h + \
            40*slip_length)) - 60*h**7*delta_t*(31*h + \
            588*slip_length)*viscosity - 1680*h**5*delta_t**2*(277*h + \
            2640*slip_length)*viscosity**2 - \
            80640*h**3*delta_t**3*(382*h + \
            1905*slip_length)*viscosity**3 - \
            3628800*h*delta_t**4*(145*h + 308*slip_length)*viscosity**4 \
            - 1117670400*delta_t**5*viscosity**5))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[6][0] = 0
            S_inv[6][1] = (-11*h**6*delta_t*viscosity*(h**6 + \
            1512*h**4*delta_t*viscosity + \
            277200*h**2*delta_t**2*viscosity**2 + \
            8648640*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[6][2] = (11*h**4*delta_t*viscosity*(-(h**7*(h + \
            4*slip_length)) - 56*h**5*delta_t*(29*h + \
            108*slip_length)*viscosity - 25200*h**3*delta_t**2*(13*h + \
            44*slip_length)*viscosity**2 - 665280*h*delta_t**3*(19*h + \
            52*slip_length)*viscosity**3 - \
            34594560*delta_t**4*viscosity**4))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[6][3] = (-11*h**8*delta_t*viscosity*(h**4 + \
            1188*h**2*delta_t*viscosity + \
            144144*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[6][4] = (11*h**2*delta_t*viscosity*(-(h**9*(h + \
            24*slip_length)) - 12*h**7*delta_t*(187*h + \
            3044*slip_length)*viscosity - 1680*h**5*delta_t**2*(403*h + \
            4176*slip_length)*viscosity**2 - \
            241920*h**3*delta_t**3*(216*h + \
            1133*slip_length)*viscosity**3 - \
            7983360*h*delta_t**4*(121*h + 260*slip_length)*viscosity**4 \
            - 2075673600*delta_t**5*viscosity**5))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[6][5] = (-11*h**10*delta_t*viscosity*(h**2 + \
            572*delta_t*viscosity))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[6][6] = (h**2*(h**11*slip_length + 2*h**9*delta_t*(19*h + \
            846*slip_length)*viscosity + 420*h**7*delta_t**2*(79*h + \
            1312*slip_length)*viscosity**2 + \
            67200*h**5*delta_t**3*(100*h + \
            909*slip_length)*viscosity**3 + \
            120960*h**3*delta_t**4*(3409*h + \
            16720*slip_length)*viscosity**4 + \
            159667200*h*delta_t**5*(43*h + 91*slip_length)*viscosity**5 \
            + 14529715200*delta_t**6*viscosity**6))/(h**13*slip_length \
            + 49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[6][7] = \
            (-11*h**12*delta_t*viscosity)/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[7][0] = 0
            S_inv[7][1] = (-13*h**6*delta_t*viscosity*(h**6 + \
            840*h**4*delta_t*viscosity + \
            75600*h**2*delta_t**2*viscosity**2 + \
            665280*delta_t**3*viscosity**3))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[7][2] = (-13*h**8*delta_t*viscosity*(h**4 + \
            756*h**2*delta_t*viscosity + \
            55440*delta_t**2*viscosity**2))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[7][3] = (13*h**4*delta_t*viscosity*(-(h**7*(h + \
            12*slip_length)) - 48*h**5*delta_t*(23*h + \
            210*slip_length)*viscosity - 5040*h**3*delta_t**2*(29*h + \
            180*slip_length)*viscosity**2 - 60480*h*delta_t**3*(59*h + \
            132*slip_length)*viscosity**3 - \
            7983360*delta_t**4*viscosity**4))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[7][4] = (-13*h**10*delta_t*viscosity*(h**2 + \
            396*delta_t*viscosity))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[7][5] = (13*h**2*delta_t*viscosity*(-(h**9*(h + \
            40*slip_length)) - 60*h**7*delta_t*(31*h + \
            588*slip_length)*viscosity - 1680*h**5*delta_t**2*(277*h + \
            2640*slip_length)*viscosity**2 - \
            80640*h**3*delta_t**3*(382*h + \
            1905*slip_length)*viscosity**3 - \
            3628800*h*delta_t**4*(145*h + 308*slip_length)*viscosity**4 \
            - 1117670400*delta_t**5*viscosity**5))/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[7][6] = \
            (-13*h**12*delta_t*viscosity)/(h**13*slip_length + \
            49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
            S_inv[7][7] = (h**2*(h**11*slip_length + 36*h**9*delta_t*(h + \
            35*slip_length)*viscosity + 420*h**7*delta_t**2*(67*h + \
            1056*slip_length)*viscosity**2 + \
            40320*h**5*delta_t**3*(128*h + \
            1119*slip_length)*viscosity**3 + \
            1814400*h**3*delta_t**4*(163*h + \
            784*slip_length)*viscosity**4 + \
            101606400*h*delta_t**5*(47*h + 99*slip_length)*viscosity**5 \
            + 10059033600*delta_t**6*viscosity**6))/(h**13*slip_length \
            + 49*h**11*delta_t*(h + 48*slip_length)*viscosity + \
            2352*h**9*delta_t**2*(31*h + 690*slip_length)*viscosity**2 \
            + 211680*h**7*delta_t**3*(129*h + \
            1684*slip_length)*viscosity**3 + \
            120960*h**5*delta_t**4*(30161*h + \
            241626*slip_length)*viscosity**4 + \
            279417600*h**3*delta_t**5*(647*h + \
            3000*slip_length)*viscosity**5 + \
            20118067200*h*delta_t**6*(137*h + \
            286*slip_length)*viscosity**6 + \
            5753767219200*delta_t**7*viscosity**7)
        
        else:
            print("This order is not implemented for the speed up of linear_source=True during implicit Euler of the SWME1D-PDE.")
        
        return S_inv
    
    def compute_source_term_lastentry(self,
                            order: int,
                            values: np.array,
                            last_moment_zero: bool,
                            **kwargs) -> np.array:
        
        viscosity   = kwargs["viscosity"]   if "viscosity"   in kwargs else self.viscosity
        slip_length = kwargs["slip_length"] if "slip_length" in kwargs else self.slip_length
        g           = kwargs["g"]           if "g"           in kwargs else 1
    
        value_out = 0
        h = values[0]
        um = values[1]/values[0]

        if order == 1:
            alpha1 = values[2]/values[0]
            if last_moment_zero:
                alpha1 = 0
            value_out = -3*viscosity/slip_length*(um + (1 + 4*slip_length/h)*alpha1)

        elif order == 2:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            
            if last_moment_zero:
                alpha2 = 0

            value_out = -5*viscosity/slip_length*(um + alpha1 + (1 + 12*slip_length/h)*alpha2)
        
        elif order == 3:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            if last_moment_zero:
                alpha3 = 0

            value_out = -7*viscosity/slip_length*((h + 4*slip_length)*alpha1 + h*(um + alpha2) + (h + 24*slip_length)*alpha3)/h
        
        elif order == 4:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            
            if last_moment_zero:
                alpha4 = 0
            
            value_out = (-9*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
            12*slip_length)*alpha2 + \
            40*slip_length*alpha4)/h))/slip_length

        elif order == 5:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]
            
            if last_moment_zero:
                alpha5 = 0
            
            value_out = (-11*viscosity*(h*um + (h + 4*slip_length)*alpha1 + \
            h*alpha2 + (h + 24*slip_length)*alpha3 + h*alpha4 + (h + \
            60*slip_length)*alpha5))/(h*slip_length)

        elif order == 6:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]
            alpha6 = values[7]/values[0]
            
            if last_moment_zero:
                alpha6 = 0
            
            value_out = (-13*viscosity*(um + alpha1 + alpha3 + alpha4 + alpha5 + \
            alpha6 + ((h + 12*slip_length)*alpha2 + \
            40*slip_length*alpha4 + \
            84*slip_length*alpha6)/h))/slip_length
        
        else:
             print("This order is not implemented for the source term last entry of the SWME1D-PDE.")

        return np.abs(value_out)

    def get_initial_values(self,
                           order: int,
                           initial_condition: str,
                           position: float) -> np.array:
        
        initial_values = np.zeros(2+order)
        
        if initial_condition == 'constantHeight_noVelocity':
            initial_values[0] = 1
            initial_values[1] = 0
            if order > 0:
                initial_values[2] = 0 
            if order > 1:
                initial_values[3] = 0 
            if order > 2:
                initial_values[4] = 0 
            if order > 3:
                initial_values[5] = 0 
            if order > 4:
                initial_values[6] = 0 
            if order > 5:
                initial_values[7] = 0
        
        elif initial_condition == 'constantHeight_constantVelocity':
            initial_values[0] = 1
            initial_values[1] = 1*initial_values[0]
            if order > 0:
                initial_values[2] = 0 
            if order > 1:
                initial_values[3] = 0 
            if order > 2:
                initial_values[4] = 0 
            if order > 3:
                initial_values[5] = 0 
            if order > 4:
                initial_values[6] = 0 
            if order > 5:
                initial_values[7] = 0
        
        elif initial_condition == 'linearHeight_noVelocity':
            initial_values[0] = 1 + 0.1*position
            initial_values[1] = 0*initial_values[0]
            if order > 0:
                initial_values[2] = 0 
            if order > 1:
                initial_values[3] = 0 
            if order > 2:
                initial_values[4] = 0 
            if order > 3:
                initial_values[5] = 0 
            if order > 4:
                initial_values[6] = 0 
            if order > 5:
                initial_values[7] = 0
        
        elif initial_condition == 'damBreak_noVelocity':
            x0 = 0
            if position < x0:
                initial_values[0] = 2
                initial_values[1] = 0*initial_values[0]
                if order > 0:
                    initial_values[2] = 0 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
            else:
                initial_values[0] = 1
                initial_values[1] = 0*initial_values[0]
                if order > 0:
                    initial_values[2] = 0 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0 
                if order > 5:
                    initial_values[7] = 0
        
        elif initial_condition == 'damBreak_constantVelocity':
            x0 = 0
            if position < x0:
                initial_values[0] = 3
                initial_values[1] = 0.25*initial_values[0]
                if order > 0:
                    initial_values[2] = 0
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
            else:
                initial_values[0] = 1
                initial_values[1] = 0.25*initial_values[0]
                if order > 0:
                    initial_values[2] = 0
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0

        elif initial_condition == 'lowDamBreak_linearVelocity':
            x0 = 0
            if position < x0:
                initial_values[0] = 1.5
                initial_values[1] = 0.25*initial_values[0]
                if order > 0:
                    initial_values[2] = -0.25*initial_values[0] 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
            else:
                initial_values[0] = 1
                initial_values[1] = 0.25*initial_values[0]
                if order > 0:
                    initial_values[2] = -0.25*initial_values[0] 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0
        
        elif initial_condition == 'highDamBreak_linearVelocity':
            x0 = 0
            if position < x0:
                initial_values[0] = 5
                initial_values[1] = 0.25*initial_values[0]
                if order > 0:
                    initial_values[2] = -0.25*initial_values[0] 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
            else:
                initial_values[0] = 1
                initial_values[1] = 0.25*initial_values[0]
                if order > 0:
                    initial_values[2] = -0.25*initial_values[0] 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
        
        elif initial_condition == 'symmetric_damBreak':
            x0 = -2
            x1 = 2
            if x0 < position < x1:
                initial_values[0] = 2
                initial_values[1] = 0*initial_values[0]
                if order > 0:
                    initial_values[2] = 0 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
            else:
                initial_values[0] = 1
                initial_values[1] = 0*initial_values[0]
                if order > 0:
                    initial_values[2] = 0 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
        
        elif initial_condition == 'linearDamBreak_noVelocity':
            x0 = -4
            x1 = 4
            if x0 < position < x1:
                initial_values[0] = 2 + (position + 4)/8.0
                initial_values[1] = 0*initial_values[0]
                if order > 0:
                    initial_values[2] = 0 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
            else:
                initial_values[0] = 2
                initial_values[1] = 0*initial_values[0]
                if order > 0:
                    initial_values[2] = 0 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0 
                if order > 5:
                    initial_values[7] = 0 
        
        elif initial_condition == 'smoothWave_noVelocity':
            initial_values[0] = 3 + 3*np.exp(-1.5*position**2)
            initial_values[1] = 0*initial_values[0]
            if order > 0:
                initial_values[2] = 0
            if order > 1:
                initial_values[3] = 0 
            if order > 2:
                initial_values[4] = 0 
            if order > 3:
                initial_values[5] = 0 
            if order > 4:
                initial_values[6] = 0
            if order > 5:
                initial_values[7] = 0 
        
        elif initial_condition == 'smoothWave_constantVelocity':
            initial_values[0] = 1 + 0.5*np.exp(-15*position**2)
            initial_values[1] = 0.2*initial_values[0]
            if order > 0:
                initial_values[2] = 0
            if order > 1:
                initial_values[3] = 0 
            if order > 2:
                initial_values[4] = 0 
            if order > 3:
                initial_values[5] = 0 
            if order > 4:
                initial_values[6] = 0
            if order > 5:
                initial_values[7] = 0 
        
        elif initial_condition == 'smoothWave_linearVelocity':
            initial_values[0] = 1 + np.exp(3*np.cos(np.pi*(position + 0.5)))/np.exp(4)
            initial_values[1] = 0.25*initial_values[0]
            if order > 0:
                initial_values[2] = -0.25*initial_values[0]
            if order > 1:
                initial_values[3] = 0 
            if order > 2:
                initial_values[4] = 0 
            if order > 3:
                initial_values[5] = 0 
            if order > 4:
                initial_values[6] = 0
            if order > 5:
                initial_values[7] = 0 
        
        elif initial_condition == 'smooth_plus_damBreak':
            x0 = -7
            x1 = 7
            if position < x0:
                initial_values[0] = 4
                initial_values[1] = 0.05*initial_values[0]
                if order > 0:
                    initial_values[2] = -0.01*initial_values[0] 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0*initial_values[0] 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0
                if order > 5:
                    initial_values[7] = 0 
            else:
                initial_values[0] = 3 + np.exp(-1.5*(position-x1)**2)
                initial_values[1] = 0.05*initial_values[0]
                if order > 0:
                    initial_values[2] = -0.01*initial_values[0] 
                if order > 1:
                    initial_values[3] = 0 
                if order > 2:
                    initial_values[4] = 0*initial_values[0] 
                if order > 3:
                    initial_values[5] = 0 
                if order > 4:
                    initial_values[6] = 0 
                if order > 5:
                    initial_values[7] = 0 
        
        else:
            print("This initial condition is not implemented yet for the SWME1D")
        
        return initial_values
    
    def compute_number_of_variables(self,
                                    order: int) -> int:
        number_of_variables = order + 2
        return int(number_of_variables)
    
    def compute_vertical_velocity_profile(self,
                                          order: int, 
                                          values: np.array,
                                          z_points: np.array) -> np.array:
        """
        reconstructs the vertical velocity profile from the moment values and evaluates the velocity profile pointwise

        Parameters
        ----------
        order: integer
            order of the model
        values: np.array (2D)
            2D numpy array containing the values of the variables in each mesh cell
        z_points: 
            the locations in vertical direction in which the velocity is computed
        
        Returns
        -------
        velocity_profile: numpy 2D array
            lateral velocity evaluated in in each point in z_points in z-direction

        """
        velocity_profile = np.zeros((len(values), len(z_points)))
        if order >= 0:
            for i in range(len(values)):
                velocity_profile[i,:] += values[i,2]*(np.ones(len(z_points)))
        if order >= 1:
            for i in range(len(values)):
                velocity_profile[i,:] += values[i,3]*(np.ones(len(z_points)) - 2*z_points)
        if order >= 2:
            for i in range(len(values)):
                velocity_profile[i,:] += values[i,4]*(np.ones(len(z_points)) - 6*z_points + 6*np.square(z_points))
        if order >= 3:
            for i in range(len(values)):
                velocity_profile[i,:] += values[i,5]*(np.ones(len(z_points)) - 12*z_points + 30*np.square(z_points) - 20*np.power(z_points,3))
        if order >= 4:
            for i in range(len(values)):
                velocity_profile[i,:] += values[i,6]*(np.ones(len(z_points)) - 20*z_points + 90*np.square(z_points) - \
                                                      140*np.power(z_points,3) + 70*np.power(z_points,4))
        if order >= 5:
            for i in range(len(values)):
                velocity_profile[i,:] += values[i,7]*(np.ones(len(z_points)) - 30*z_points + 210*np.square(z_points) - \
                                                      560*np.power(z_points,3) + 630*np.power(z_points,4) - 252*np.power(z_points,4))     
        if order >= 6:
            for i in range(len(values)):
                velocity_profile[i,:] += values[i,8]*(np.ones(len(z_points)) - 42*z_points + 420*np.square(z_points) - \
                                                      1680*np.power(z_points,3) + 3150*np.power(z_points,4) - \
                                                      2772*np.power(z_points,5) + 924*np.power(z_points,6))
                   
        return velocity_profile
    
    def compute_all_breakdown_criteria(self,
                                   values: np.array,
                                   orders: list,
                                   number_of_variables: list,
                                   n,
                                   delta_x,
                                   tolerance_up_height_gradient = 0.3,
                                   tolerance_down_height_gradient = 0.03,
                                   tolerance_up_momentum_gradient = 0.2,
                                   tolerance_down_momentum_gradient = 0.02,
                                   tolerance_up_moment_gradient = 0.1,
                                   tolerance_down_moment_gradient = 0.01,
                                   tolerance_up_last_moment = 0.01,
                                   tolerance_down_last_moment = 0.001,
                                   tolerance_up_source = 0.002,
                                   tolerance_down_source = 0.0002,
                                   **kwargs) -> np.array:

        """
        Compute ALL breakdown criteria for quantifying the required modelling complexity

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
        max_order = max(orders)
        tolerances_up = np.zeros(4+max_order)
        tolerances_down = np.zeros(4+max_order)
        tolerances_up[0] = tolerance_up_last_moment
        tolerances_up[1] = tolerance_up_source
        tolerances_up[2] = tolerance_up_height_gradient
        tolerances_up[3] = tolerance_up_momentum_gradient
        tolerances_down[0] = tolerance_down_last_moment
        tolerances_down[1] = tolerance_down_source
        tolerances_down[2] = tolerance_down_height_gradient
        tolerances_down[3] = tolerance_down_momentum_gradient
        for i in range(max_order):
            tolerances_up[4+i] = tolerance_up_moment_gradient
            tolerances_down[4+i] = tolerance_down_moment_gradient

        breakdown_criterion_flags = np.zeros(n)
        source_term_lastentry = np.zeros(n)
        for i in range(n):
            source_term_lastentry[i] = self.compute_source_term_lastentry(orders[i+1],values[i+1,:number_of_variables[i+1]],False,**kwargs)

        breakdown_estimators = np.zeros((n,max_order+4))
        for i in range(n):
            breakdown_estimators[i,0] = np.abs(self.compute_source_term_lastentry(orders[i+1],values[i+1,:number_of_variables[i+1]],False,**kwargs))
            breakdown_estimators[i,1] = np.abs(values[i+1,number_of_variables[i+1]-1])
            breakdown_estimators[i,2] = np.abs((values[i+1,0] - values[i,0]))/delta_x
            breakdown_estimators[i,3] = np.abs((values[i+1,1] - values[i,1]))/delta_x
            for j in range(orders[i+1]):
                breakdown_estimators[i,4+j] = np.abs((values[i+1,2+j])-values[i,2+j])/delta_x
        breakdown_estimators[0,0] = np.abs(self.compute_source_term_lastentry(orders[1],values[1,:number_of_variables[1]],False,**kwargs))
        breakdown_estimators[0,1] = np.abs(values[1,number_of_variables[1]-1])
        breakdown_estimators[0,2] = np.abs((values[2,0] - values[1,0]))/delta_x
        breakdown_estimators[0,3] = np.abs((values[2,1] - values[1,1]))/delta_x
        for j in range(orders[i+1]):
            breakdown_estimators[0,4+j] = np.abs((values[2,2+j])-values[2,2+j])/delta_x    

        for i in range(n):
            if breakdown_estimators[i,0] > tolerances_up[0]:
                for j in range(2,orders[i+1]+4):
                    if breakdown_estimators[i,j] > tolerances_up[j]:
                        breakdown_criterion_flags[i] = 1
                        break
            elif breakdown_estimators[i,1] > tolerances_up[1]:
                breakdown_criterion_flags[i] = 1
            elif (breakdown_criterion_flags[i] !=1 and\
                  breakdown_estimators[i,0] < tolerances_down[0] and breakdown_estimators[i,1] < tolerances_down[1]):
                breakdown_criterion_flags[i] = -1
                for j in range(2,orders[i+1]+4):
                    if breakdown_estimators[i,j] > tolerances_down[j]:
                        breakdown_criterion_flags[i] = 0
                        break

        # #TODO: this could be made much more efficient
        # #TODO: with break statement, this might be faster. I could stack the tolerance in an array, then for-loop, and then use break.
        # flagged = False
        # if (np.abs(values[1,number_of_variables[1]-1]) > tolerance_up_last_moment\
        #     or np.abs(source_term_lastentry[0] > tolerance_up_source)\
        #     or np.abs((values[2,0] - values[1,0])) > tolerance_up_height_gradient \
        #     or np.abs((values[2,1] - values[1,1])) > tolerance_up_momentum_gradient):
        #     breakdown_criterion_flags[0] = 1
        #     flagged = True
        # elif (not flagged and\
        #     (np.abs(values[1,number_of_variables[1]-1]) < tolerance_down_last_moment\
        #     or np.abs(source_term_lastentry[0] < tolerance_down_source)\
        #     or np.abs((values[2,0] - values[1,0])) < tolerance_down_height_gradient \
        #     or np.abs((values[2,1] - values[1,1])) < tolerance_down_momentum_gradient)):
        #     breakdown_criterion_flags[0] = -1
        #     flagged = True
        # j = 2
        # while j<number_of_variables[1] and not flagged:
        #     if np.abs((values[2,j] - values[1,j])) > tolerance_up_moment_gradient:
        #         breakdown_criterion_flags[0] = 1
        #         flagged = True
        #     elif not flagged and\
        #         np.abs((values[2,j] - values[1,j])) < tolerance_down_moment_gradient:
        #         breakdown_criterion_flags[0] = -1
        #         flagged = True
        #     j += 1
        # if not flagged:
        #     if source_term_lastentry[0] > tolerance_up_source:
        #         breakdown_criterion_flags[0] = 1
        #         flagged = True
        #     elif source_term_lastentry[0] < tolerance_down_source:
        #         breakdown_criterion_flags[0] = -1
        #         flagged = True               

        # flagged = False
        # for i in range(1,n): 
        #     if np.abs(source_term_lastentry[i]) > tolerance_up_source:
        #         breakdown_criterion_flags[i] = 1
        #         flagged = True
        #     elif np.abs(source_term_lastentry[i]) < tolerance_down_source:
        #         breakdown_criterion_flags[i] = -1
        #         flagged = True  
        #     if not flagged and (np.abs((values[i+1,0] - values[i,0])) > tolerance_up_height_gradient \
        #         or np.abs((values[i+1,1] - values[1,1])) > tolerance_up_momentum_gradient \
        #         or np.abs(values[i+1,number_of_variables[i+1]-1]) > tolerance_up_last_moment):
        #         breakdown_criterion_flags[i] = 1
        #         flagged = True
        #     elif (not flagged and\
        #         (np.abs((values[i+1,0] - values[i,0])) < tolerance_down_height_gradient \
        #         or np.abs((values[i+1,1] - values[i,1])) < tolerance_down_momentum_gradient \
        #         or np.abs(values[i+1,number_of_variables[i+1]-1]) < tolerance_down_last_moment)):
        #         breakdown_criterion_flags[i] = -1
        #         flagged = True
        #     j = 2
        #     while j<number_of_variables[i+1] and not flagged:
        #         if np.abs((values[i+1,j] - values[i,j])) > tolerance_up_moment_gradient:
        #             breakdown_criterion_flags[i] = 1
        #             flagged = True
        #         elif not flagged and\
        #             np.abs((values[i+1,j] - values[i,j])) < tolerance_down_moment_gradient:
        #             breakdown_criterion_flags[i] = -1
        #             flagged = True
        #         j += 1

        return breakdown_criterion_flags
    
    def compute_breakdown_criterion(self,
                                   values: np.array,
                                   orders: np.array,
                                   number_of_variables: list,
                                   breakdown_criterion: str,
                                   n,
                                   delta_x) -> np.array:
        breakdown_criterion_values = np.zeros(n)
        
        if breakdown_criterion == 'height_gradient':
            if np.abs(values[1,0]) < 0.001:
                breakdown_criterion_values[0] = np.abs((values[2,0] - values[1,0]))/delta_x
            else:
                breakdown_criterion_values[0] = np.abs((values[2,0] - values[1,0]))/delta_x
            for i in range(2,n): 
                if np.abs(values[i,0]) < 0.001:
                    #breakdown_criterion_values[i] = np.abs((values[i+1,0] - values[i,0])/0.001)
                    breakdown_criterion_values[i] = np.abs((values[i+1,0] - values[i,0]))/delta_x
                else:
                    #breakdown_criterion_values[i] = np.abs((values[i+1,0] - values[i,0])/values[i,0])
                    breakdown_criterion_values[i] = np.abs((values[i+1,0] - values[i,0]))/delta_x
            if np.abs(values[n,0]) < 0.001:
                breakdown_criterion_values[-1] = np.abs((values[n,0] - values[n-1,0]))/delta_x
            else:
                breakdown_criterion_values[-1] = np.abs((values[n,0] - values[n-1,0]))/delta_x
        elif breakdown_criterion == 'momentum_gradient':
            if np.abs(values[1,1]) < 0.001:
                breakdown_criterion_values[0] = np.abs((values[2,1] - values[1,1]))/delta_x
            else:
                breakdown_criterion_values[0] = np.abs((values[2,1] - values[1,1]))/delta_x
            for i in range(2,n): 
                if np.abs(values[i,1]) < 0.001:
                    #breakdown_criterion_values[i] = np.abs((values[i+1,1] - values[i,1])/0.001)
                    breakdown_criterion_values[i] = np.abs((values[i+1,1] - values[i,1]))/delta_x
                else:
                    #breakdown_criterion_values[i] = np.abs((values[i+1,1] - values[i,1])/values[i,1])
                    breakdown_criterion_values[i] = np.abs((values[i+1,1] - values[i,1]))/delta_x
            if np.abs(values[n,1]) < 0.001:
                breakdown_criterion_values[-1] = np.abs((values[n,1] - values[n-1,1]))/delta_x
            else:
                breakdown_criterion_values[-1] = np.abs((values[n,1] - values[n-1,1]))/delta_x
        elif breakdown_criterion == 'last_moment':
            for i in range(1,n+1): 
                breakdown_criterion_values[i-1] = np.abs(values[i,number_of_variables[i]-1]) 
                # If the order is 0, there are no moments and the above value is never used     
        elif breakdown_criterion == 'source_term':
            for i in range(n):
                breakdown_criterion_values[i] = np.abs(self.compute_source_term_lastentry(orders[i+1],values[i+1,:number_of_variables[i+1]],True,**kwargs))
        else:
            print('this criterion is not implemented yet')  

        return breakdown_criterion_values
    

class VegetationSWME1D(SWME1D):
    """
    This class represents the SWME1D with vegetation drag term in the momentum equation.

    ...

    Attributes
    ----------
    initial_condition : str
        initial condition for the SWME1D
    viscosity : float
        value for the dynamic viscosity
    slip_length : float
        value for the slip length
    hyperbolic : boolean
        whether the model is hyperbolic, true (HSWME) or false (SWME)
    diameter : float
        the diameter of the vegetation stem
    CD : float
        the drag force coefficient
    surface_density : float
        number of vegetation elements per squared meter
    
    Methods inherited from interface SWME1D
    ---------------------------------------
    def compute_system_matrix(self,order,values):
        computes the system matrix of the SWME1D evaluated in the given values, for the given order. 
    def get_initial_values(self,order,initial_condition,position):
        calculates the initial values for one specific physical position
    def compute_number_of_variables(self,order):
        computes the number of state variables in the PDE given the order of the moment model
    def compute_breakdown_criterion(self,values,breakdown_criterion,n)
        computes the values of the given breakdown criterion in each mesh cell
    def compute_vertical_velocity_profile(self,values):
        reconstruct the vertical velocity profiles from the moment values

    Implemented methods from interface SWME1D
    ------------------------------------------
    def compute_source_term(self,order,values):
        computes the system matrix of the SWME1D evaluated in the given values, for the given order.

    Instance methods
    ----------------
    def compute_single_legendre_integral(self,order,h_r):
        computes the integrals of the legendre polynomials phi_i on the interval [0,h_r]
    def compute_double_legendre_integral(self,order,h_r):
        computes the integrals of the products of legendre polynomials phi_i*phi_j on the interval [0,h_r]
    def compute_triple_legendre_integral(self,order,h_r):
        computes the integrals of the legendre polynomials phi_i*phi_j*phi_k on the interval [0,h_r]
    def compute_drag_force(self,order,values,h_v,stem_diam,n_stems,drag_coeff):
        computes the drag force caused by vegetation
    """

    def __init__(self, 
                initial_condition: str,
                viscosity: float,
                slip_length: float,
                hyperbolic: bool,
                diameter: float,
                CD: float,
                surface_density: float):
        """
        Constructs all the necessary attributes for the VegetationSWME1D object.

        Parameters
        ----------
        initial_condition : str
            initial condition of the PDE
        viscosity : float
            dynamic viscosity value
        slip_length : float
            slip length value
        hyperbolic : boolean
            true if hyperbolic, false if not hyperbolic
        diameter : float
            diameter of a cylindrical vegetation stem
        CD : float
            drag coefficient
        surface_density : float
            number of vegetation stems per squared meter
        """
        self.initial_condition = initial_condition
        self.viscosity = viscosity
        self.slip_length = slip_length
        self.hyperbolic = hyperbolic
        self.diameter = diameter
        self.CD = CD
        self.surface_density = surface_density

    def compute_single_legendre_integral(self,
                                         order: int,
                                         h_r: float) -> np.array:
        """
        computes the integrals of the legendre polynomials phi_i on the interval [0,h_r]

        Parameters
        ----------
        order : integer
            order of the model
        h_r : float
            relative height of the vegetation
        
        Returns
        -------
        integrals: numpy 1D array
            values of the integral for each Legendre polynomial phi_i, i = 1,...,order

        """

        integrals = np.zeros(order)
        if order >= 1:
            integrals[0] = h_r - h_r**2
        if order >= 2:
            integrals[1] = h_r - 3*h_r**2 + 2*h_r**3
        if order >= 3:
            integrals[2] = h_r - 6*h_r**2 + 10*h_r**3 - 5*h_r**4
        if order >= 4:
            integrals[3] = h_r - 10*h_r**2 + 30*h_r**3 - 35*h_r**4 + 14*h_r**5
        if order >= 5:
            integrals[4] = h_r - 15*h_r**2 + 70*h_r**3 - 140*h_r**4 + 126*h_r**5 \
            - 42*h_r**6
        if order >= 6:
            integrals[5] = h_r - 21*h_r**2 + 140*h_r**3 - 420*h_r**4 + 630*h_r**5 \
            - 462*h_r**6 + 132*h_r**7
        if order >= 7:
            integrals[6] = h_r - 28*h_r**2 + 252*h_r**3 - 1050*h_r**4 + \
            2310*h_r**5 - 2772*h_r**6 + 1716*h_r**7 - 429*h_r**8

        return integrals

    def compute_double_legendre_integral(self,
                                         order: int,
                                         h_r: float) -> np.array:
        """
        computes the integrals of the product of legendre polynomials phi_i*phi_j on the interval [0,h_r]

        Parameters
        ----------
        order : integer
            order of the model
        h_r : float
            relative height of the vegetation
        
        Returns
        -------
        integrals: numpy 2D array
            values of the integral for each product combination of Legendre polynomials phi_i*phi_j, i,j = 1,...,order

        """

        integrals = np.zeros((order,order))

        if order >= 1:
            integrals[0][0] = h_r - 2*h_r**2 + (4*h_r**3)/3.
        if order >= 2:
            integrals[1][0] = integrals[0][1] = h_r - 4*h_r**2 + 6*h_r**3 - \
            3*h_r**4
            integrals[1][1] = h_r - 6*h_r**2 + 16*h_r**3 - 18*h_r**4 + \
            (36*h_r**5)/5.
        if order >= 3:
            integrals[2][0] = integrals[0][2] = h_r - 7*h_r**2 + 18*h_r**3 - \
            20*h_r**4 + 8*h_r**5
            integrals[2][1] = integrals[1][2] = h_r - 9*h_r**2 + 36*h_r**3 - \
            68*h_r**4 + 60*h_r**5 - 20*h_r**6
            integrals[2][2] = h_r - 12*h_r**2 + 68*h_r**3 - 190*h_r**4 + \
            276*h_r**5 - 200*h_r**6 + (400*h_r**7)/7.
        if order >= 4:
            integrals[3][0] = integrals[0][3] = h_r - 11*h_r**2 + (130*h_r**3)/3. \
            - 80*h_r**4 + 70*h_r**5 - (70*h_r**6)/3.
            integrals[3][1] = integrals[1][3] = h_r - 13*h_r**2 + 72*h_r**3 - \
            200*h_r**4 + 290*h_r**5 - 210*h_r**6 + 60*h_r**7
            integrals[3][2] = integrals[2][3] = h_r - 16*h_r**2 + 120*h_r**3 - \
            460*h_r**4 + 970*h_r**5 - 1140*h_r**6 + 700*h_r**7 - 175*h_r**8
            integrals[3][3] = h_r - 20*h_r**2 + (580*h_r**3)/3. - 970*h_r**4 + \
            2768*h_r**5 - (14000*h_r**6)/3. + 4600*h_r**7 - 2450*h_r**8 + \
            (4900*h_r**9)/9.
        if order >= 5:
            integrals[4][0] = integrals[0][4] = h_r - 16*h_r**2 + 90*h_r**3 - \
            245*h_r**4 + 350*h_r**5 - 252*h_r**6 + 72*h_r**7
            integrals[4][1] = integrals[1][4] = h_r - 18*h_r**2 + 132*h_r**3 - \
            500*h_r**4 + 1050*h_r**5 - 1232*h_r**6 + 756*h_r**7 - 189*h_r**8
            integrals[4][2] = integrals[2][4] = h_r - 21*h_r**2 + 200*h_r**3 - \
            1000*h_r**4 + 2850*h_r**5 - 4802*h_r**6 + 4732*h_r**7 - 2520*h_r**8 + \
            560*h_r**9
            integrals[4][3] = integrals[3][4] = h_r - 25*h_r**2 + 300*h_r**3 - \
            1900*h_r**4 + 7000*h_r**5 - 15792*h_r**6 + 22120*h_r**7 - \
            18760*h_r**8 + 8820*h_r**9 - 1764*h_r**10
            integrals[4][4] = h_r - 30*h_r**2 + 440*h_r**3 - 3430*h_r**4 + \
            15792*h_r**5 - 45584*h_r**6 + 84760*h_r**7 - 101430*h_r**8 + \
            75460*h_r**9 - 31752*h_r**10 + (63504*h_r**11)/11.
        if order >= 6:
            integrals[5][0] = integrals[0][5] = h_r - 22*h_r**2 + 168*h_r**3 - \
            630*h_r**4 + 1302*h_r**5 - 1512*h_r**6 + 924*h_r**7 - 231*h_r**8
            integrals[5][1] = integrals[1][5] = h_r - 24*h_r**2 + 226*h_r**3 - \
            1113*h_r**4 + 3150*h_r**5 - 5292*h_r**6 + 5208*h_r**7 - 2772*h_r**8 + \
            616*h_r**9
            integrals[5][2] = integrals[2][5] = h_r - 27*h_r**2 + 318*h_r**3 - \
            2000*h_r**4 + 7350*h_r**5 - 16562*h_r**6 + 23184*h_r**7 - \
            19656*h_r**8 + 9240*h_r**9 - 1848*h_r**10
            integrals[5][3] = integrals[3][5] = h_r - 31*h_r**2 + 450*h_r**3 - \
            3500*h_r**4 + 16100*h_r**5 - 46452*h_r**6 + 86352*h_r**7 - \
            103320*h_r**8 + 76860*h_r**9 - 32340*h_r**10 + 5880*h_r**11
            integrals[5][4] = integrals[4][5] = h_r - 36*h_r**2 + 630*h_r**3 - \
            5915*h_r**4 + 33180*h_r**5 - 118664*h_r**6 + 280224*h_r**7 - \
            442260*h_r**8 + 461580*h_r**9 - 305760*h_r**10 + 116424*h_r**11 - \
            19404*h_r**12
            integrals[5][5] = h_r - 42*h_r**2 + 868*h_r**3 - 9660*h_r**4 + \
            64764*h_r**5 - 280224*h_r**6 + 814728*h_r**7 - 1623762*h_r**8 + \
            2223620*h_r**9 - 2056824*h_r**10 + 1227744*h_r**11 - 426888*h_r**12 + \
            (853776*h_r**13)/13.
        if order >= 7:
            integrals[6][0] = integrals[0][6] = h_r - 29*h_r**2 + (868*h_r**3)/3. \
            - 1428*h_r**4 + 3990*h_r**5 - 6622*h_r**6 + 6468*h_r**7 - 3432*h_r**8 \
            + (2288*h_r**9)/3.
            integrals[6][1] = integrals[1][6] = h_r - 31*h_r**2 + 366*h_r**3 - \
            2268*h_r**4 + (41286*h_r**5)/5. - 18522*h_r**6 + 25872*h_r**7 - \
            21912*h_r**8 + 10296*h_r**9 - (10296*h_r**10)/5.
            integrals[6][2] = integrals[2][6] = h_r - 34*h_r**2 + 486*h_r**3 - \
            3743*h_r**4 + 17150*h_r**5 - 49392*h_r**6 + 91728*h_r**7 - \
            109692*h_r**8 + 81576*h_r**9 - 34320*h_r**10 + 6240*h_r**11
            integrals[6][3] = integrals[3][6] = h_r - 38*h_r**2 + \
            (1966*h_r**3)/3. - 6125*h_r**4 + 34300*h_r**5 - (367696*h_r**6)/3. + \
            289296*h_r**7 - 456444*h_r**8 + 476300*h_r**9 - 315480*h_r**10 + \
            120120*h_r**11 - 20020*h_r**12
            integrals[6][4] = integrals[4][6] = h_r - 43*h_r**2 + 882*h_r**3 - \
            9800*h_r**4 + 65660*h_r**5 - 284004*h_r**6 + 825552*h_r**7 - \
            1645128*h_r**8 + 2252700*h_r**9 - 2083620*h_r**10 + 1243704*h_r**11 - \
            432432*h_r**12 + 66528*h_r**13
            integrals[6][5] = integrals[5][6] = h_r - 49*h_r**2 + 1176*h_r**3 - \
            15288*h_r**4 + 120540*h_r**5 - 619164*h_r**6 + 2165016*h_r**7 - \
            5284344*h_r**8 + 9094932*h_r**9 - 10990980*h_r**10 + 9125424*h_r**11 \
            - 4956336*h_r**12 + 1585584*h_r**13 - 226512*h_r**14
            integrals[6][6] = h_r - 56*h_r**2 + (4648*h_r**3)/3. - 23268*h_r**4 + \
            (1065036*h_r**5)/5. - 1279544*h_r**6 + 5284344*h_r**7 - \
            15439974*h_r**8 + (97219276*h_r**9)/3. - (245144592*h_r**10)/5. + \
            52993584*h_r**11 - 39903864*h_r**12 + 19880784*h_r**13 - \
            5889312*h_r**14 + (3926208*h_r**15)/5.

        return integrals

    def compute_triple_legendre_integral(self,
                                         order: int,
                                         h_r: float) -> np.array:
        """
        computes the integrals of the product of legendre polynomials phi_i*phi_j*phi_k on the interval [0,h_r]

        Parameters
        ----------
        order : integer
            order of the model
        h_r : float
            relative height of the vegetation
        
        Returns
        -------
        integrals: numpy 3D array
            values of the integral for each combination of triple product of 
            Legendre polynomials phi_i*phi_j*phi_k, i,j,k = 1,...,order

        """

        integrals = np.zeros((order,order,order))
        if order >= 1:
            integrals[0,0,0] = h_r - 3*h_r**2 + 4*h_r**3 - 2*h_r**4
        if order >= 2:
            integrals[1,0,0] = integrals[0,1,0] = integrals[0,0,1] = h_r - \
            5*h_r**2 + (34*h_r**3)/3. - 12*h_r**4 + (24*h_r**5)/5.
            integrals[1,1,0] = integrals[1,0,1] = integrals[0,1,1] = h_r - \
            7*h_r**2 + 24*h_r**3 - 42*h_r**4 + 36*h_r**5 - 12*h_r**6
            integrals[1,1,1] = h_r - 9*h_r**2 + 42*h_r**3 - 108*h_r**4 + \
            (756*h_r**5)/5. - 108*h_r**6 + (216*h_r**7)/7.
        if order >= 3:
            integrals[2,0,0] = integrals[0,2,0] = integrals[0,0,2] = h_r - \
            8*h_r**2 + (82*h_r**3)/3. - 47*h_r**4 + 40*h_r**5 - (40*h_r**6)/3.
            integrals[2,1,0] = integrals[2,0,1] = integrals[1,2,0] = \
            integrals[1,0,2] = integrals[0,2,1] = integrals[0,1,2] = h_r - \
            10*h_r**2 + 48*h_r**3 - 122*h_r**4 + (844*h_r**5)/5. - 120*h_r**6 + \
            (240*h_r**7)/7.
            integrals[2,1,1] = integrals[1,2,1] = integrals[1,1,2] = h_r - \
            12*h_r**2 + 74*h_r**3 - 257*h_r**4 + 516*h_r**5 - 592*h_r**6 + \
            360*h_r**7 - 90*h_r**8
            integrals[2,2,0] = integrals[2,0,2] = integrals[0,2,2] = h_r - \
            13*h_r**2 + 84*h_r**3 - 292*h_r**4 + 580*h_r**5 - 660*h_r**6 + \
            400*h_r**7 - 100*h_r**8
            integrals[2,2,1] = integrals[2,1,2] = integrals[1,2,2] = h_r - \
            15*h_r**2 + 118*h_r**3 - 532*h_r**4 + (7164*h_r**5)/5. - 2340*h_r**6 \
            + (15880*h_r**7)/7. - 1200*h_r**8 + (800*h_r**9)/3.
            integrals[2,2,2] = h_r - 18*h_r**2 + 174*h_r**3 - 987*h_r**4 + \
            3420*h_r**5 - 7440*h_r**6 + 10200*h_r**7 - 8550*h_r**8 + 4000*h_r**9 \
            - 800*h_r**10
        if order >= 4:
            integrals[3,0,0] = integrals[0,3,0] = integrals[0,0,3] = h_r - \
            12*h_r**2 + 58*h_r**3 - 145*h_r**4 + 198*h_r**5 - 140*h_r**6 + \
            40*h_r**7
            integrals[3,1,0] = integrals[3,0,1] = integrals[1,3,0] = \
            integrals[1,0,3] = integrals[0,3,1] = integrals[0,1,3] = h_r - \
            14*h_r**2 + (268*h_r**3)/3. - 308*h_r**4 + 610*h_r**5 - \
            (2080*h_r**6)/3. + 420*h_r**7 - 105*h_r**8
            integrals[3,1,1] = integrals[1,3,1] = integrals[1,1,3] = h_r - \
            16*h_r**2 + 126*h_r**3 - 563*h_r**4 + (7546*h_r**5)/5. - 2460*h_r**6 \
            + (16680*h_r**7)/7. - 1260*h_r**8 + 280*h_r**9
            integrals[3,2,0] = integrals[3,0,2] = integrals[2,3,0] = \
            integrals[2,0,3] = integrals[0,3,2] = integrals[0,2,3] = h_r - \
            17*h_r**2 + (424*h_r**3)/3. - 640*h_r**4 + 1706*h_r**5 - \
            (8270*h_r**6)/3. + (18580*h_r**7)/7. - 1400*h_r**8 + (2800*h_r**9)/9.
            integrals[3,2,1] = integrals[3,1,2] = integrals[2,3,1] = \
            integrals[2,1,3] = integrals[1,3,2] = integrals[1,2,3] = h_r - \
            19*h_r**2 + 186*h_r**3 - 1048*h_r**4 + 3610*h_r**5 - 7830*h_r**6 + \
            10720*h_r**7 - 8980*h_r**8 + 4200*h_r**9 - 840*h_r**10
            integrals[3,2,2] = integrals[2,3,2] = integrals[2,2,3] = h_r - \
            22*h_r**2 + 258*h_r**3 - 1785*h_r**4 + 7674*h_r**5 - 21240*h_r**6 + \
            (269280*h_r**7)/7. - 45300*h_r**8 + 33400*h_r**9 - 14000*h_r**10 + \
            (28000*h_r**11)/11.
            integrals[3,3,0] = integrals[3,0,3] = integrals[0,3,3] = h_r - \
            21*h_r**2 + 220*h_r**3 - 1260*h_r**4 + 4320*h_r**5 - 9280*h_r**6 + \
            12600*h_r**7 - 10500*h_r**8 + 4900*h_r**9 - 980*h_r**10
            integrals[3,3,1] = integrals[3,1,3] = integrals[1,3,3] = h_r - \
            23*h_r**2 + (826*h_r**3)/3. - 1900*h_r**4 + 8120*h_r**5 - \
            (67160*h_r**6)/3. + (283240*h_r**7)/7. - 47600*h_r**8 + \
            (315700*h_r**9)/9. - 14700*h_r**10 + (29400*h_r**11)/11.
            integrals[3,3,2] = integrals[3,2,3] = integrals[2,3,3] = h_r - \
            26*h_r**2 + (1090*h_r**3)/3. - 3015*h_r**4 + 15720*h_r**5 - \
            53680*h_r**6 + 123000*h_r**7 - 190350*h_r**8 + (588700*h_r**9)/3. - \
            129080*h_r**10 + 49000*h_r**11 - (24500*h_r**12)/3.
            integrals[3,3,3] = h_r - 30*h_r**2 + 490*h_r**3 - 4805*h_r**4 + \
            29862*h_r**5 - 123000*h_r**6 + (2421600*h_r**7)/7. - 674100*h_r**8 + \
            909300*h_r**9 - 833000*h_r**10 + (5439000*h_r**11)/11. - \
            171500*h_r**12 + (343000*h_r**13)/13.
        if order >= 5:
            integrals[4,0,0] = integrals[0,4,0] = integrals[0,0,4] = h_r - \
            17*h_r**2 + (334*h_r**3)/3. - 380*h_r**4 + 742*h_r**5 - \
            (2506*h_r**6)/3. + 504*h_r**7 - 126*h_r**8
            integrals[4,1,0] = integrals[4,0,1] = integrals[1,4,0] = \
            integrals[1,0,4] = integrals[0,4,1] = integrals[0,1,4] = h_r - \
            19*h_r**2 + 156*h_r**3 - 698*h_r**4 + 1850*h_r**5 - 2982*h_r**6 + \
            2868*h_r**7 - 1512*h_r**8 + 336*h_r**9
            integrals[4,1,1] = integrals[1,4,1] = integrals[1,1,4] = h_r - \
            21*h_r**2 + 206*h_r**3 - 1148*h_r**4 + (19626*h_r**5)/5. - \
            8482*h_r**6 + 11592*h_r**7 - 9702*h_r**8 + 4536*h_r**9 - \
            (4536*h_r**10)/5.
            integrals[4,2,0] = integrals[4,0,2] = integrals[2,4,0] = \
            integrals[2,0,4] = integrals[0,4,2] = integrals[0,2,4] = h_r - \
            22*h_r**2 + 228*h_r**3 - 1300*h_r**4 + 4450*h_r**5 - 9552*h_r**6 + \
            12964*h_r**7 - 10801*h_r**8 + 5040*h_r**9 - 1008*h_r**10
            integrals[4,2,1] = integrals[4,1,2] = integrals[2,4,1] = \
            integrals[2,1,4] = integrals[1,4,2] = integrals[1,2,4] = h_r - \
            24*h_r**2 + 286*h_r**3 - 1963*h_r**4 + 8370*h_r**5 - 23052*h_r**6 + \
            (291496*h_r**7)/7. - 48972*h_r**8 + (108248*h_r**9)/3. - \
            15120*h_r**10 + (30240*h_r**11)/11.
            integrals[4,2,2] = integrals[2,4,2] = integrals[2,2,4] = h_r - \
            27*h_r**2 + 378*h_r**3 - 3120*h_r**4 + 16218*h_r**5 - 55302*h_r**6 + \
            126624*h_r**7 - 195876*h_r**8 + 201880*h_r**9 - 132776*h_r**10 + \
            50400*h_r**11 - 8400*h_r**12
            integrals[4,3,0] = integrals[4,0,3] = integrals[3,4,0] = \
            integrals[3,0,4] = integrals[0,4,3] = integrals[0,3,4] = h_r - \
            26*h_r**2 + (1000*h_r**3)/3. - 2350*h_r**4 + 10040*h_r**5 - \
            (82376*h_r**6)/3. + 49192*h_r**7 - 57470*h_r**8 + (379540*h_r**9)/9. \
            - 17640*h_r**10 + (35280*h_r**11)/11.
            integrals[4,3,1] = integrals[4,1,3] = integrals[3,4,1] = \
            integrals[3,1,4] = integrals[1,4,3] = integrals[1,3,4] = h_r - \
            28*h_r**2 + 402*h_r**3 - 3325*h_r**4 + 17200*h_r**5 - 58392*h_r**6 + \
            133336*h_r**7 - 205954*h_r**8 + 212100*h_r**9 - 139440*h_r**10 + \
            52920*h_r**11 - 8820*h_r**12
            integrals[4,3,2] = integrals[4,2,3] = integrals[3,4,2] = \
            integrals[3,2,4] = integrals[2,4,3] = integrals[2,3,4] = h_r - \
            31*h_r**2 + 510*h_r**3 - 4980*h_r**4 + 30840*h_r**5 - 126792*h_r**6 + \
            (2493864*h_r**7)/7. - 693840*h_r**8 + 935620*h_r**9 - 856940*h_r**10 \
            + (5594680*h_r**11)/11. - 176400*h_r**12 + (352800*h_r**13)/13.
            integrals[4,3,3] = integrals[3,4,3] = integrals[3,3,4] = h_r - \
            35*h_r**2 + (1990*h_r**3)/3. - 7560*h_r**4 + 55014*h_r**5 - \
            268042*h_r**6 + 903840*h_r**7 - 2150820*h_r**8 + (10910620*h_r**9)/3. \
            - 4342268*h_r**10 + 3577000*h_r**11 - (5801600*h_r**12)/3. + \
            617400*h_r**13 - 88200*h_r**14
            integrals[4,4,0] = integrals[4,0,4] = integrals[0,4,4] = h_r - \
            31*h_r**2 + 480*h_r**3 - 4090*h_r**4 + 21280*h_r**5 - 71904*h_r**6 + \
            162904*h_r**7 - 249760*h_r**8 + 255780*h_r**9 - 167580*h_r**10 + \
            63504*h_r**11 - 10584*h_r**12
            integrals[4,4,1] = integrals[4,1,4] = integrals[1,4,4] = h_r - \
            33*h_r**2 + 562*h_r**3 - 5500*h_r**4 + 33840*h_r**5 - 138264*h_r**6 + \
            386872*h_r**7 - 751548*h_r**8 + (3035900*h_r**9)/3. - 926100*h_r**10 \
            + (6043464*h_r**11)/11. - 190512*h_r**12 + (381024*h_r**13)/13.
            integrals[4,4,2] = integrals[4,2,4] = integrals[2,4,4] = h_r - \
            36*h_r**2 + 690*h_r**3 - 7845*h_r**4 + 56880*h_r**5 - 276504*h_r**6 + \
            931224*h_r**7 - 2214450*h_r**8 + 3742900*h_r**9 - 4467680*h_r**10 + \
            3679704*h_r**11 - 1989204*h_r**12 + 635040*h_r**13 - 90720*h_r**14
            integrals[4,4,3] = integrals[4,3,4] = integrals[3,4,4] = h_r - \
            40*h_r**2 + 870*h_r**3 - 11415*h_r**4 + 96126*h_r**5 - 546084*h_r**6 \
            + 2169000*h_r**7 - 6163500*h_r**8 + 12680500*h_r**9 - \
            18914000*h_r**10 + (222670504*h_r**11)/11. - 15143940*h_r**12 + \
            (97707960*h_r**13)/13. - 2222640*h_r**14 + 296352*h_r**15
            integrals[4,4,4] = h_r - 45*h_r**2 + 1110*h_r**3 - 16620*h_r**4 + \
            160398*h_r**5 - 1050126*h_r**6 + 4844280*h_r**7 - 16155090*h_r**8 + \
            39553780*h_r**9 - 71555876*h_r**10 + 95393592*h_r**11 - \
            92468880*h_r**12 + 63345240*h_r**13 - 29053080*h_r**14 + \
            8001504*h_r**15 - 1000188*h_r**16
        if order >= 6:
            integrals[5,0,0] = integrals[0,5,0] = integrals[0,0,5] = h_r - \
            23*h_r**2 + (592*h_r**3)/3. - 882*h_r**4 + 2310*h_r**5 - 3682*h_r**6 \
            + 3516*h_r**7 - 1848*h_r**8 + (1232*h_r**9)/3.
            integrals[5,1,0] = integrals[5,0,1] = integrals[1,5,0] = \
            integrals[1,0,5] = integrals[0,5,1] = integrals[0,1,5] = h_r - \
            25*h_r**2 + 258*h_r**3 - 1452*h_r**4 + (24654*h_r**5)/5. - \
            10542*h_r**6 + 14280*h_r**7 - 11886*h_r**8 + 5544*h_r**9 - \
            (5544*h_r**10)/5.
            integrals[5,1,1] = integrals[1,5,1] = integrals[1,1,5] = h_r - \
            27*h_r**2 + 324*h_r**3 - 2202*h_r**4 + 9306*h_r**5 - 25494*h_r**6 + \
            45924*h_r**7 - 53928*h_r**8 + 39704*h_r**9 - 16632*h_r**10 + 3024*h_r**11
            integrals[5,2,0] = integrals[5,0,2] = integrals[2,5,0] = \
            integrals[2,0,5] = integrals[0,5,2] = integrals[0,2,5] = h_r - \
            28*h_r**2 + 354*h_r**3 - 2477*h_r**4 + 10550*h_r**5 - 28812*h_r**6 + \
            51576*h_r**7 - 60228*h_r**8 + 44184*h_r**9 - 18480*h_r**10 + 3360*h_r**11
            integrals[5,2,1] = integrals[5,1,2] = integrals[2,5,1] = \
            integrals[2,1,5] = integrals[1,5,2] = integrals[1,2,5] = h_r - \
            30*h_r**2 + 428*h_r**3 - 3512*h_r**4 + (90474*h_r**5)/5. - \
            61312*h_r**6 + 139860*h_r**7 - 215901*h_r**8 + 222264*h_r**9 - \
            (730464*h_r**10)/5. + 55440*h_r**11 - 9240*h_r**12
            integrals[5,2,2] = integrals[2,5,2] = integrals[2,2,5] = h_r - \
            33*h_r**2 + 544*h_r**3 - 5272*h_r**4 + 32490*h_r**5 - 133242*h_r**6 + \
            (2617252*h_r**7)/7. - 727608*h_r**8 + (2942072*h_r**9)/3. - \
            897960*h_r**10 + (5861520*h_r**11)/11. - 184800*h_r**12 + \
            (369600*h_r**13)/13.
            integrals[5,3,0] = integrals[5,0,3] = integrals[3,5,0] = \
            integrals[3,0,5] = integrals[0,5,3] = integrals[0,3,5] = h_r - \
            32*h_r**2 + (1474*h_r**3)/3. - 4175*h_r**4 + 21700*h_r**5 - \
            (219856*h_r**6)/3. + 165984*h_r**7 - 254436*h_r**8 + 260540*h_r**9 - \
            170688*h_r**10 + 64680*h_r**11 - 10780*h_r**12
            integrals[5,3,1] = integrals[5,1,3] = integrals[3,5,1] = \
            integrals[3,1,5] = integrals[1,5,3] = integrals[1,3,5] = h_r - \
            34*h_r**2 + 576*h_r**3 - 5618*h_r**4 + 34520*h_r**5 - 140952*h_r**6 + \
            394248*h_r**7 - 765702*h_r**8 + 1030876*h_r**9 - 943320*h_r**10 + \
            (6155520*h_r**11)/11. - 194040*h_r**12 + (388080*h_r**13)/13.
            integrals[5,3,2] = integrals[5,2,3] = integrals[3,5,2] = \
            integrals[3,2,5] = integrals[2,5,3] = integrals[2,3,5] = h_r - \
            37*h_r**2 + 708*h_r**3 - 8020*h_r**4 + 58048*h_r**5 - 281952*h_r**6 + \
            949144*h_r**7 - 2256436*h_r**8 + 3813180*h_r**9 - 4551036*h_r**10 + \
            3748080*h_r**11 - 2026080*h_r**12 + 646800*h_r**13 - 92400*h_r**14
            integrals[5,3,3] = integrals[3,5,3] = integrals[3,3,5] = h_r - \
            41*h_r**2 + (2680*h_r**3)/3. - 11680*h_r**4 + 98150*h_r**5 - \
            (1671026*h_r**6)/3. + 2211172*h_r**7 - 6281240*h_r**8 + \
            (116279380*h_r**9)/9. - 19268340*h_r**10 + (226820160*h_r**11)/11. - \
            15425200*h_r**12 + (99519000*h_r**13)/13. - 2263800*h_r**14 + \
            301840*h_r**15
            integrals[5,4,0] = integrals[5,0,4] = integrals[4,5,0] = \
            integrals[4,0,5] = integrals[0,5,4] = integrals[0,4,5] = h_r - \
            37*h_r**2 + 678*h_r**3 - 6860*h_r**4 + 42644*h_r**5 - 173964*h_r**6 + \
            483648*h_r**7 - 932652*h_r**8 + 1247820*h_r**9 - 1136604*h_r**10 + \
            (7395864*h_r**11)/11. - 232848*h_r**12 + (465696*h_r**13)/13.
            integrals[5,4,1] = integrals[5,1,4] = integrals[4,5,1] = \
            integrals[4,1,5] = integrals[1,5,4] = integrals[1,4,5] = h_r - \
            39*h_r**2 + 776*h_r**3 - 8858*h_r**4 + 63840*h_r**5 - 308224*h_r**6 + \
            1032696*h_r**7 - 2447424*h_r**8 + 4128012*h_r**9 - 4921140*h_r**10 + \
            4050144*h_r**11 - 2188536*h_r**12 + 698544*h_r**13 - 99792*h_r**14
            integrals[5,4,2] = integrals[5,2,4] = integrals[4,5,2] = \
            integrals[4,2,5] = integrals[2,5,4] = integrals[2,4,5] = h_r - \
            42*h_r**2 + 928*h_r**3 - 12130*h_r**4 + 101592*h_r**5 - 575064*h_r**6 \
            + 2279368*h_r**7 - 6469302*h_r**8 + (39899300*h_r**9)/3. - \
            19828200*h_r**10 + (233360064*h_r**11)/11. - 15867768*h_r**12 + \
            (102366096*h_r**13)/13. - 2328480*h_r**14 + 310464*h_r**15
            integrals[5,4,3] = integrals[5,3,4] = integrals[4,5,3] = \
            integrals[4,3,5] = integrals[3,5,4] = integrals[3,4,5] = h_r - \
            46*h_r**2 + 1140*h_r**3 - 17020*h_r**4 + 163870*h_r**5 - \
            1071504*h_r**6 + 4939564*h_r**7 - 16466275*h_r**8 + 40305300*h_r**9 - \
            72902760*h_r**10 + 97177584*h_r**11 - 94190544*h_r**12 + \
            64521240*h_r**13 - 29591520*h_r**14 + 8149680*h_r**15 - 1018710*h_r**16
            integrals[5,4,4] = integrals[4,5,4] = integrals[4,4,5] = h_r - \
            51*h_r**2 + 1420*h_r**3 - 24010*h_r**4 + 262710*h_r**5 - \
            1960266*h_r**6 + 10374076*h_r**7 - 40023480*h_r**8 + \
            (343805420*h_r**9)/3. - 245951580*h_r**10 + (4361927472*h_r**11)/11. \
            - 477534792*h_r**12 + (5496935640*h_r**13)/13. - 267087240*h_r**14 + \
            113848560*h_r**15 - 29338848*h_r**16 + (58677696*h_r**17)/17.
            integrals[5,5,0] = integrals[5,0,5] = integrals[0,5,5] = h_r - \
            43*h_r**2 + 924*h_r**3 - 10962*h_r**4 + 80220*h_r**5 - 388164*h_r**6 \
            + 1295112*h_r**7 - 3049536*h_r**8 + 5110308*h_r**9 - 6059340*h_r**10 \
            + 4967424*h_r**11 - 2677752*h_r**12 + 853776*h_r**13 - 121968*h_r**14
            integrals[5,5,1] = integrals[5,1,5] = integrals[1,5,5] = h_r - \
            45*h_r**2 + 1038*h_r**3 - 13692*h_r**4 + (571284*h_r**5)/5. - \
            642684*h_r**6 + 2533440*h_r**7 - 7162092*h_r**8 + 14685748*h_r**9 - \
            (109292148*h_r**10)/5. + (256990104*h_r**11)/11. - 17463600*h_r**12 + \
            (112620816*h_r**13)/13. - 2561328*h_r**14 + (1707552*h_r**15)/5.
            integrals[5,5,2] = integrals[5,2,5] = integrals[2,5,5] = h_r - \
            48*h_r**2 + 1214*h_r**3 - 18107*h_r**4 + 173460*h_r**5 - \
            1129744*h_r**6 + 5195232*h_r**7 - 17292996*h_r**8 + 42290388*h_r**9 - \
            76448400*h_r**10 + 101863944*h_r**11 - 98706972*h_r**12 + \
            67603536*h_r**13 - 31002048*h_r**14 + 8537760*h_r**15 - 1067220*h_r**16
            integrals[5,5,3] = integrals[5,3,5] = integrals[3,5,5] = h_r - \
            52*h_r**2 + 1458*h_r**3 - 24605*h_r**4 + 268562*h_r**5 - \
            2000964*h_r**6 + 10580808*h_r**7 - 40801572*h_r**8 + 116794300*h_r**9 \
            - 250605600*h_r**10 + (4443832344*h_r**11)/11. - 486452988*h_r**12 + \
            (5599207656*h_r**13)/13. - 272044080*h_r**14 + 115958304*h_r**15 - \
            29882160*h_r**16 + (59764320*h_r**17)/17.
            integrals[5,5,4] = integrals[5,4,5] = integrals[4,5,5] = h_r - \
            57*h_r**2 + 1778*h_r**3 - 33740*h_r**4 + 415506*h_r**5 - \
            3503626*h_r**6 + 21063672*h_r**7 - 92937978*h_r**8 + 306962460*h_r**9 \
            - 768344892*h_r**10 + 1465845192*h_r**11 - 2129799504*h_r**12 + \
            2337712776*h_r**13 - 1904673960*h_r**14 + 1115962848*h_r**15 - \
            444293388*h_r**16 + 107575776*h_r**17 - 11952864*h_r**18
            integrals[5,5,5] = h_r - 63*h_r**2 + 2184*h_r**3 - 46242*h_r**4 + \
            636930*h_r**5 - 6025446*h_r**6 + 40810788*h_r**7 - 203964264*h_r**8 + \
            768344892*h_r**9 - 2212770420*h_r**10 + (54030663792*h_r**11)/11. - \
            8426123496*h_r**12 + (144760659624*h_r**13)/13. - 11224554264*h_r**14 \
            + 8467517520*h_r**15 - 4625758368*h_r**16 + (29368186848*h_r**17)/17. \
            - 394444512*h_r**18 + (788889024*h_r**19)/19.
        if order >= 7:
            integrals[6,0,0] = integrals[0,6,0] = integrals[0,0,6] = h_r - \
            30*h_r**2 + 328*h_r**3 - 1862*h_r**4 + (31374*h_r**5)/5. - \
            13272*h_r**6 + 17820*h_r**7 - 14751*h_r**8 + 6864*h_r**9 - \
            (6864*h_r**10)/5.
            integrals[6,1,0] = integrals[6,0,1] = integrals[1,6,0] = \
            integrals[1,0,6] = integrals[0,6,1] = integrals[0,1,6] = h_r - \
            32*h_r**2 + (1222*h_r**3)/3. - 2817*h_r**4 + 11886*h_r**5 - \
            32284*h_r**6 + 57624*h_r**7 - 67188*h_r**8 + (147752*h_r**9)/3. - \
            20592*h_r**10 + 3744*h_r**11
            integrals[6,1,1] = integrals[1,6,1] = integrals[1,1,6] = h_r - \
            34*h_r**2 + 492*h_r**3 - 4008*h_r**4 + (102306*h_r**5)/5. - \
            68880*h_r**6 + 156516*h_r**7 - 241089*h_r**8 + 247896*h_r**9 - \
            (814176*h_r**10)/5. + 61776*h_r**11 - 10296*h_r**12
            integrals[6,2,0] = integrals[6,0,2] = integrals[2,6,0] = \
            integrals[2,0,6] = integrals[0,6,2] = integrals[0,2,6] = h_r - \
            35*h_r**2 + (1594*h_r**3)/3. - 4472*h_r**4 + (115694*h_r**5)/5. - \
            (233926*h_r**6)/3. + 176400*h_r**7 - 270216*h_r**8 + 276584*h_r**9 - \
            (905784*h_r**10)/5. + 68640*h_r**11 - 11440*h_r**12
            integrals[6,2,1] = integrals[6,1,2] = integrals[2,6,1] = \
            integrals[2,1,6] = integrals[1,6,2] = integrals[1,2,6] = h_r - \
            37*h_r**2 + 624*h_r**3 - 6032*h_r**4 + 36866*h_r**5 - 150114*h_r**6 + \
            419244*h_r**7 - 813528*h_r**8 + 1094664*h_r**9 - 1001352*h_r**10 + \
            593904*h_r**11 - 205920*h_r**12 + 31680*h_r**13
            integrals[6,2,2] = integrals[2,6,2] = integrals[2,2,6] = h_r - \
            40*h_r**2 + 768*h_r**3 - 8632*h_r**4 + (310514*h_r**5)/5. - \
            300612*h_r**6 + (7070260*h_r**7)/7. - 2398531*h_r**8 + 4050504*h_r**9 \
            - (24160704*h_r**10)/5. + 3978480*h_r**11 - 2150280*h_r**12 + \
            686400*h_r**13 - (686400*h_r**14)/7.
            integrals[6,3,0] = integrals[6,0,3] = integrals[3,6,0] = \
            integrals[3,0,6] = integrals[0,6,3] = integrals[0,3,6] = h_r - \
            39*h_r**2 + 706*h_r**3 - 7108*h_r**4 + 44100*h_r**5 - 179732*h_r**6 + \
            499408*h_r**7 - 962712*h_r**8 + 1287756*h_r**9 - 1172820*h_r**10 + \
            693720*h_r**11 - 240240*h_r**12 + 36960*h_r**13
            integrals[6,3,1] = integrals[6,1,3] = integrals[3,6,1] = \
            integrals[3,1,6] = integrals[1,6,3] = integrals[1,3,6] = h_r - \
            41*h_r**2 + (2428*h_r**3)/3. - 9188*h_r**4 + (330296*h_r**5)/5. - \
            (955696*h_r**6)/3. + 1066632*h_r**7 - 2526792*h_r**8 + 4260716*h_r**9 \
            - (25392156*h_r**10)/5. + 4179120*h_r**11 - 2258080*h_r**12 + \
            720720*h_r**13 - 102960*h_r**14
            integrals[6,3,2] = integrals[6,2,3] = integrals[3,6,2] = \
            integrals[3,2,6] = integrals[2,6,3] = integrals[2,3,6] = h_r - \
            44*h_r**2 + (2908*h_r**3)/3. - 12598*h_r**4 + 105200*h_r**5 - \
            (1783856*h_r**6)/3. + 2354968*h_r**7 - 6680522*h_r**8 + \
            (123565324*h_r**9)/9. - 20464320*h_r**10 + (240810960*h_r**11)/11. - \
            16372840*h_r**12 + (105618480*h_r**13)/13. - 2402400*h_r**14 + \
            320320*h_r**15
            integrals[6,3,3] = integrals[3,6,3] = integrals[3,3,6] = h_r - \
            48*h_r**2 + 1192*h_r**3 - 17700*h_r**4 + 169830*h_r**5 - \
            1108492*h_r**6 + 5105076*h_r**7 - 17007909*h_r**8 + 41614900*h_r**9 - \
            75251600*h_r**10 + 100290240*h_r**11 - 97195440*h_r**12 + \
            66574200*h_r**13 - 30531600*h_r**14 + 8408400*h_r**15 - 1051050*h_r**16
            integrals[6,4,0] = integrals[6,0,4] = integrals[4,6,0] = \
            integrals[4,0,6] = integrals[0,6,4] = integrals[0,4,6] = h_r - \
            44*h_r**2 + (2818*h_r**3)/3. - 11123*h_r**4 + 81340*h_r**5 - \
            (1180312*h_r**6)/3. + 1312416*h_r**7 - 3089844*h_r**8 + \
            5177372*h_r**9 - 6138480*h_r**10 + 5032104*h_r**11 - 2712556*h_r**12 \
            + 864864*h_r**13 - 123552*h_r**14
            integrals[6,4,1] = integrals[6,1,4] = integrals[4,6,1] = \
            integrals[4,1,6] = integrals[1,6,4] = integrals[1,4,6] = h_r - \
            46*h_r**2 + 1056*h_r**3 - 13898*h_r**4 + (579376*h_r**5)/5. - \
            651504*h_r**6 + 2567544*h_r**7 - 7257294*h_r**8 + 14879292*h_r**9 - \
            (110724072*h_r**10)/5. + (260343744*h_r**11)/11. - 17690904*h_r**12 + \
            (114084432*h_r**13)/13. - 2594592*h_r**14 + (1729728*h_r**15)/5.
            integrals[6,4,2] = integrals[6,2,4] = integrals[4,6,2] = \
            integrals[4,2,6] = integrals[2,6,4] = integrals[2,4,6] = h_r - \
            49*h_r**2 + 1236*h_r**3 - 18388*h_r**4 + 175960*h_r**5 - \
            1145424*h_r**6 + 5265736*h_r**7 - 17524264*h_r**8 + 42850332*h_r**9 - \
            77453580*h_r**10 + 103196784*h_r**11 - 99994176*h_r**12 + \
            68483184*h_r**13 - 31404912*h_r**14 + 8648640*h_r**15 - 1081080*h_r**16
            integrals[6,4,3] = integrals[6,3,4] = integrals[4,6,3] = \
            integrals[4,3,6] = integrals[3,6,4] = integrals[3,4,6] = h_r - \
            53*h_r**2 + (4456*h_r**3)/3. - 25000*h_r**4 + 272510*h_r**5 - \
            (6087242*h_r**6)/3. + 10725652*h_r**7 - 41350808*h_r**8 + \
            (1065136900*h_r**9)/9. - 253913700*h_r**10 + (4502155584*h_r**11)/11. \
            - 492811312*h_r**12 + (5672184504*h_r**13)/13. - 275583000*h_r**14 + \
            117465040*h_r**15 - 30270240*h_r**16 + (60540480*h_r**17)/17.
            integrals[6,4,4] = integrals[4,6,4] = integrals[4,4,6] = h_r - \
            58*h_r**2 + 1812*h_r**3 - 34300*h_r**4 + 421750*h_r**5 - \
            3553536*h_r**6 + 21354748*h_r**7 - 94197847*h_r**8 + 311069700*h_r**9 \
            - 778531080*h_r**10 + 1485149232*h_r**11 - 2157711696*h_r**12 + \
            2368242744*h_r**13 - 1929487200*h_r**14 + 1130477040*h_r**15 - \
            450066078*h_r**16 + 108972864*h_r**17 - 12108096*h_r**18
            integrals[6,5,0] = integrals[6,0,5] = integrals[5,6,0] = \
            integrals[5,0,6] = integrals[0,6,5] = integrals[0,5,6] = h_r - \
            50*h_r**2 + (3724*h_r**3)/3. - 17052*h_r**4 + (725004*h_r**5)/5. - \
            820064*h_r**6 + 3226440*h_r**7 - 9073122*h_r**8 + \
            (55467964*h_r**9)/3. - (136809288*h_r**10)/5. + 29109024*h_r**11 - \
            21686280*h_r**12 + (139564656*h_r**13)/13. - 3171168*h_r**14 + \
            (2114112*h_r**15)/5.
            integrals[6,5,1] = integrals[6,1,5] = integrals[5,6,1] = \
            integrals[5,1,6] = integrals[1,6,5] = integrals[1,5,6] = h_r - \
            52*h_r**2 + 1374*h_r**3 - 20727*h_r**4 + 198156*h_r**5 - \
            1283016*h_r**6 + 5865888*h_r**7 - 19436916*h_r**8 + 47381508*h_r**9 - \
            85468464*h_r**10 + 113724072*h_r**11 - 110101068*h_r**12 + \
            75365136*h_r**13 - 34550208*h_r**14 + 9513504*h_r**15 - 1189188*h_r**16
            integrals[6,5,2] = integrals[6,2,5] = integrals[5,6,2] = \
            integrals[5,2,6] = integrals[2,6,5] = integrals[2,5,6] = h_r - \
            55*h_r**2 + 1578*h_r**3 - 26612*h_r**4 + (1444324*h_r**5)/5. - \
            2142084*h_r**6 + 11291280*h_r**7 - 43454952*h_r**8 + 124233828*h_r**9 \
            - (1331753628*h_r**10)/5. + (4720423944*h_r**11)/11. - \
            516534480*h_r**12 + (5943910896*h_r**13)/13. - 288742608*h_r**14 + \
            (615317472*h_r**15)/5. - 31711680*h_r**16 + (63423360*h_r**17)/17.
            integrals[6,5,3] = integrals[6,3,5] = integrals[5,6,3] = \
            integrals[5,3,6] = integrals[3,6,5] = integrals[3,5,6] = h_r - \
            59*h_r**2 + (5578*h_r**3)/3. - 35168*h_r**4 + 431410*h_r**5 - \
            (10886722*h_r**6)/3. + 21786576*h_r**7 - 96048024*h_r**8 + \
            317066252*h_r**9 - 793348980*h_r**10 + 1513163064*h_r**11 - \
            2198149856*h_r**12 + 2412421704*h_r**13 - 1965364632*h_r**14 + \
            1151451840*h_r**15 - 458405640*h_r**16 + 110990880*h_r**17 - \
            12332320*h_r**18
            integrals[6,5,4] = integrals[6,4,5] = integrals[5,6,4] = \
            integrals[5,4,6] = integrals[4,6,5] = integrals[4,5,6] = h_r - \
            64*h_r**2 + 2226*h_r**3 - 47033*h_r**4 + 646730*h_r**5 - \
            6112596*h_r**6 + 41380584*h_r**7 - 206750196*h_r**8 + \
            778685868*h_r**9 - 2242242240*h_r**10 + (54744852792*h_r**11)/11. - \
            8536880268*h_r**12 + (146655612936*h_r**13)/13. - 11371043376*h_r**14 \
            + 8577787680*h_r**15 - 4685910768*h_r**16 + (29749747104*h_r**17)/17. \
            - 399567168*h_r**18 + (799134336*h_r**19)/19.
            integrals[6,5,5] = integrals[5,6,5] = integrals[5,5,6] = h_r - \
            70*h_r**2 + 2688*h_r**3 - 63042*h_r**4 + (4820634*h_r**5)/5. - \
            10158456*h_r**6 + 76919940*h_r**7 - 431764713*h_r**8 + \
            1837221372*h_r**9 - (30096583272*h_r**10)/5. + 15337313712*h_r**11 - \
            30546749520*h_r**12 + 47565906696*h_r**13 - 57618912384*h_r**14 + \
            (268403191056*h_r**15)/5. - 37694976858*h_r**16 + 19286799840*h_r**17 \
            - 6782396544*h_r**18 + 1465079616*h_r**19 - (732539808*h_r**20)/5.
            integrals[6,6,0] = integrals[6,0,6] = integrals[0,6,6] = h_r - \
            57*h_r**2 + 1624*h_r**3 - 25592*h_r**4 + 250236*h_r**5 - \
            1634556*h_r**6 + 7477848*h_r**7 - 24687576*h_r**8 + 59855268*h_r**9 - \
            107360484*h_r**10 + 142137072*h_r**11 - 137058768*h_r**12 + \
            93549456*h_r**13 - 42810768*h_r**14 + 11778624*h_r**15 - 1472328*h_r**16
            integrals[6,6,1] = integrals[6,1,6] = integrals[1,6,6] = h_r - \
            59*h_r**2 + (5326*h_r**3)/3. - 30408*h_r**4 + (1651356*h_r**5)/5. - \
            2437652*h_r**6 + 12777744*h_r**7 - 48940728*h_r**8 + \
            (418239676*h_r**9)/3. - (1490677452*h_r**10)/5. + 479510136*h_r**11 - \
            576513168*h_r**12 + (6629104944*h_r**13)/13. - 321873552*h_r**14 + \
            (685727328*h_r**15)/5. - 35335872*h_r**16 + (70671744*h_r**17)/17.
            integrals[6,6,2] = integrals[6,2,6] = integrals[2,6,6] = h_r - \
            62*h_r**2 + (6022*h_r**3)/3. - 38057*h_r**4 + 464716*h_r**5 - \
            (11671408*h_r**6)/3. + 23275728*h_r**7 - 102377916*h_r**8 + \
            337461428*h_r**9 - 843558504*h_r**10 + 1607866392*h_r**11 - \
            2334652628*h_r**12 + 2561405616*h_r**13 - 2086273728*h_r**14 + \
            1222107744*h_r**15 - 486491148*h_r**16 + 117786240*h_r**17 - \
            13087360*h_r**18
            integrals[6,6,3] = integrals[6,3,6] = integrals[3,6,6] = h_r - \
            66*h_r**2 + 2326*h_r**3 - 49063*h_r**4 + (3360546*h_r**5)/5. - \
            6335504*h_r**6 + 42820624*h_r**7 - 213737724*h_r**8 + \
            804506652*h_r**9 - (11578146312*h_r**10)/5. + \
            (56520017064*h_r**11)/11. - 8811828564*h_r**12 + \
            (151355559432*h_r**13)/13. - 11734142112*h_r**14 + \
            (44254951488*h_r**15)/5. - 4834898640*h_r**16 + \
            (30694641120*h_r**17)/17. - 412251840*h_r**18 + (824503680*h_r**19)/19.
            integrals[6,6,4] = integrals[6,4,6] = integrals[4,6,6] = h_r - \
            71*h_r**2 + (8218*h_r**3)/3. - 64148*h_r**4 + (4896626*h_r**5)/5. - \
            (30923578*h_r**6)/3. + 78006096*h_r**7 - 437710536*h_r**8 + \
            1862103548*h_r**9 - (30499437108*h_r**10)/5. + 15540859992*h_r**11 - \
            30949574432*h_r**12 + 48190183272*h_r**13 - 58372437816*h_r**14 + \
            (271904069184*h_r**15)/5. - 38185716888*h_r**16 + 19537566048*h_r**17 \
            - 6870511648*h_r**18 + 1484106624*h_r**19 - (742053312*h_r**20)/5.
            integrals[6,6,5] = integrals[6,5,6] = integrals[5,6,6] = h_r - \
            77*h_r**2 + (9772*h_r**3)/3. - 84252*h_r**4 + 1423506*h_r**5 - \
            16610538*h_r**6 + 139723452*h_r**7 - 874624104*h_r**8 + \
            4170046412*h_r**9 - 15398377596*h_r**10 + (489969867744*h_r**11)/11. \
            - 101632253184*h_r**12 + (2384184343464*h_r**13)/13. - \
            261313784424*h_r**14 + 292085787312*h_r**15 - 252912112992*h_r**16 + \
            (2823570664992*h_r**17)/17. - 79910262432*h_r**18 + \
            (504463063104*h_r**19)/19. - 5441724288*h_r**20 + 518259456*h_r**21
            integrals[6,6,6] = h_r - 84*h_r**2 + 3892*h_r**3 - 110558*h_r**4 + \
            (10272906*h_r**5)/5. - 26418924*h_r**6 + 245506572*h_r**7 - \
            1703222667*h_r**8 + 9035930988*h_r**9 - (186541907328*h_r**10)/5. + \
            121395970080*h_r**11 - 313874746416*h_r**12 + 647745331464*h_r**13 - \
            1067802213456*h_r**14 + (7008603467184*h_r**15)/5. - \
            1454047747614*h_r**16 + 1175840553888*h_r**17 - 725157328896*h_r**18 \
            + 329224319424*h_r**19 - (518324238432*h_r**20)/5. + \
            20212118784*h_r**21 - 1837465344*h_r**22

        return integrals

    def compute_drag_force(self,
                           order: int,
                           values: np.array,
                           h_v: float,
                           stem_diam: float, 
                           n_stems: float, 
                           drag_coeff: float) -> np.array:
        """
        computes the drag force term caused by vegetation

        Parameters
        ----------
        order : integer
            order of the model
        values : numpy array
            values of the variables
        h_v : float
            height of the vegetation
        stem_diam : float
            diameter of a cylindrical plant stem
        n_stems : float
            number of plant stems per squared meter
        drag_coeff : float
            drag force coefficient
        
        Returns
        -------
        drag_values: numpy 1D array
            values of the drag force on each variable

        """

        drag_values = np.zeros(order+2)
        h = values[0]

        single_integrals = self.compute_single_legendre_integral(order,h_v/h)
        double_integrals = self.compute_double_legendre_integral(order,h_v/h)
        triple_integrals = self.compute_triple_legendre_integral(order,h_v/h)

        um = values[1]/values[0]
        
        sum_var = 0

        for i in range(order):
            for j in range(order):
                sum_var += values[i+2]*values[j+2]*double_integrals[i,j]
        drag_values[1] = h_v/h*um*um + 2*um*np.sum(values[2:]*single_integrals) + sum_var

        for k in range(order):
            sum_var = 0
            for i in range(order):
                for j in range(order):
                    sum_var += values[i+2]*values[j+2]*triple_integrals[i,j,k]
            drag_values[k+2] = um*um*single_integrals[k] + 2*um*np.sum(values[2:]*double_integrals[:,k]) + sum_var  
        
        drag_values = -np.sign(um)*drag_values*n_stems*stem_diam*drag_coeff

        return drag_values

    def compute_source_term(self,
                            order: int,
                            values: np.array,
                            **kwargs) -> np.array:
        
        viscosity   = kwargs["viscosity"]   if "viscosity"   in kwargs else self.viscosity
        slip_length = kwargs["slip_length"] if "slip_length" in kwargs else self.slip_length
        h_v         = kwargs["h_v"]         if "h_v"         in kwargs else 0.4
        stem_diam   = kwargs["stem_diam"]   if "stem_diam"   in kwargs else 0.008
        n_stems     = kwargs["n_stems"]     if "n_stems"     in kwargs else 800
        drag_coeff  = kwargs["drag_coeff"]  if "drag_coeff"  in kwargs else 0.97
        
        """
        Computes the source term with a given order of the PDE evaluated in the given values.

        Parameters
        ----------
        order : integer
            order of the model
        values : numpy array
            values of the variables
        h_v : float
            height of the vegetation
        stem_diam : float
            diameter of a cylindrical plant stem
        n_stems : float
            number of plant stems per squared meter
        drag_coeff : float
            drag force coefficient
        
        
        Returns
        -------
        S: numpy 1D array
            source term vector

        """

        
        S = np.zeros(order+2) 
        h = values[0]
        um = values[1]/values[0]
        if order == 0:
            S[0] = 0
            S[1] = -viscosity/slip_length*um
        if order == 1:
            alpha1 = values[2]/values[0]

            S[0] = 0
            S[1] = -viscosity/slip_length*(um + alpha1)
            S[2] = -3*viscosity/slip_length*(um + (1 + 4*slip_length/h)*alpha1)

        if order == 2:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]

            S[0] = 0
            S[1] = -viscosity/slip_length*(um + alpha1 + alpha2)
            S[2] = -3*viscosity/slip_length*(um + (1 + 4*slip_length/h)*alpha1 + alpha2)
            S[3] = -5*viscosity/slip_length*(um + alpha1 + (1 + 12*slip_length/h)*alpha2)
        if order == 3:
            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]

            S[0] = 0
            S[1] = -viscosity/slip_length*(um + alpha1 + alpha2 + alpha3)
            S[2] = -3*viscosity/slip_length*((h + 4*slip_length)*alpha1 + h*(um + alpha2) + (h + 4*slip_length)*alpha3)/h
            S[3] = -5*viscosity/slip_length*(um + alpha1 + (1 + 12*slip_length/h)*alpha2 + alpha3)
            S[4] = -7*viscosity/slip_length*((h + 4*slip_length)*alpha1 + h*(um + alpha2) + (h + 24*slip_length)*alpha3)/h
        if order == 4:

            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]

            S[0] = 0
            S[1] = -((viscosity*(um + alpha1 + alpha2 + alpha3 + \
            alpha4))/slip_length)
            S[2] = (-3*viscosity*(um + alpha2 + alpha3 + ((h + \
            4*slip_length)*alpha1 + 4*slip_length*alpha3)/h + \
            alpha4))/slip_length
            S[3] = (-5*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
            12*slip_length)*alpha2 + \
            12*slip_length*alpha4)/h))/slip_length
            S[4] = (-7*viscosity*(um + alpha2 + alpha3 + ((h + \
            4*slip_length)*alpha1 + 24*slip_length*alpha3)/h + \
            alpha4))/slip_length
            S[5] = (-9*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
            12*slip_length)*alpha2 + \
            40*slip_length*alpha4)/h))/slip_length

        if order == 5:

            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]

            S[0] = 0
            S[1] = -((viscosity*(um + alpha1 + alpha2 + alpha3 + alpha4 + \
            alpha5))/slip_length)
            S[2] = (-3*viscosity*(h*um + (h + 4*slip_length)*alpha1 + \
            h*alpha2 + (h + 4*slip_length)*alpha3 + h*alpha4 + (h + \
            4*slip_length)*alpha5))/(h*slip_length)
            S[3] = (-5*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
            12*slip_length)*alpha2 + 12*slip_length*alpha4)/h + \
            alpha5))/slip_length
            S[4] = (-7*viscosity*(h*um + (h + 4*slip_length)*alpha1 + \
            h*alpha2 + (h + 24*slip_length)*alpha3 + h*alpha4 + (h + \
            24*slip_length)*alpha5))/(h*slip_length)
            S[5] = (-9*viscosity*(um + alpha1 + alpha3 + alpha4 + ((h + \
            12*slip_length)*alpha2 + 40*slip_length*alpha4)/h + \
            alpha5))/slip_length
            S[6] = (-11*viscosity*(h*um + (h + 4*slip_length)*alpha1 + \
            h*alpha2 + (h + 24*slip_length)*alpha3 + h*alpha4 + (h + \
            60*slip_length)*alpha5))/(h*slip_length)

        if order == 6:

            alpha1 = values[2]/values[0]
            alpha2 = values[3]/values[0]
            alpha3 = values[4]/values[0]
            alpha4 = values[5]/values[0]
            alpha5 = values[6]/values[0]
            alpha6 = values[7]/values[0]

            S[0] = 0
            S[1] = -((viscosity*(um + alpha1 + alpha2 + alpha3 + alpha4 + \
            alpha5 + alpha6))/slip_length)
            S[2] = (-3*viscosity*(um + alpha2 + alpha3 + alpha4 + alpha5 + \
            ((h + 4*slip_length)*alpha1 + 4*slip_length*(alpha3 + \
            alpha5))/h + alpha6))/slip_length
            S[3] = (-5*viscosity*(um + alpha1 + alpha3 + alpha4 + alpha5 + \
            alpha6 + ((h + 12*slip_length)*alpha2 + \
            12*slip_length*(alpha4 + alpha6))/h))/slip_length
            S[4] = (-7*viscosity*(um + alpha2 + alpha3 + alpha4 + alpha5 + \
            ((h + 4*slip_length)*alpha1 + 24*slip_length*(alpha3 + \
            alpha5))/h + alpha6))/slip_length
            S[5] = (-9*viscosity*(um + alpha1 + alpha3 + alpha4 + alpha5 + \
            alpha6 + ((h + 12*slip_length)*alpha2 + \
            40*slip_length*(alpha4 + alpha6))/h))/slip_length
            S[6] = (-11*viscosity*(um + alpha2 + alpha3 + alpha4 + alpha5 + \
            ((h + 4*slip_length)*alpha1 + 12*slip_length*(2*alpha3 + \
            5*alpha5))/h + alpha6))/slip_length
            S[7] = (-13*viscosity*(um + alpha1 + alpha3 + alpha4 + alpha5 + \
            alpha6 + ((h + 12*slip_length)*alpha2 + \
            40*slip_length*alpha4 + \
            84*slip_length*alpha6)/h))/slip_length

        S = S + self.compute_drag_force(order, values, h_v, stem_diam, n_stems, drag_coeff)
        return S


class SGSWME1D(PDE):

    """
    This class represents the one-dimensional stochastic Galerkin Shallow Water Moment Equations (SGSWME1D).

    ...

    Attributes
    ----------
    initial_condition : str
        initial condition for the SGSWME1D
    distr : str
        distribution of viscosity parameter
    mu : float
        mean of distribution of viscosity parameter
    sigma : float
        standard deviation of distribution of viscosity parameter
    slip_length : float
        value for the slip length
    hyperbolic : boolean
        whether the model is hyperbolic, true (HSGSWME1D) or false (SGSWME1D)

    
    Implemented methods from interface PDE
    ---------------------------------
    def compute_system_matrix(self, mom_order, SG_order, values):
        computes the system matrix of the SGSWME1D evaluated in the given values, for the given moment order and stochastic Galerkin order. 
    def compute_source_term(self, mom_order, SG_order, values):
        computes the source term of the SWME1D evaluated in the given values, for the given moment order and stochastic Galerkin order.
    def get_initial_values(self, mom_order, SG_order, initial_condition, position):
        calculates the initial values for one specific physical position
    def compute_number_of_variables(self, mom_order, SG_order):
        computes the number of state variables in the PDE given the moment order and stochastic Galerkin order of the model

    Instance methods
    ----------------
    def compute_exp_and_var(self, mom_order, SG_order, values, primitive)
        computes the expectation and variance of h, um and alpha1 or of h, h*um and h*alpha1
    def compute_vertical_velocity_profile(self, mom_order, SG_order, values, z_points):
        reconstruct the expectation and variance of the vertical velocity profiles from the model values
    """

    def __init__(self, 
                initial_condition: str,
                distr: str,
                mu: float,
                sigma: float,
                slip_length: float,
                hyperbolic: bool):
        """
        Constructs all the necessary attributes for the SWME1D object.

        Parameters
        ----------
        initial_condition : str
            initial condition of the PDE
        distr : str
            distribution of viscosity parameter name
        mu : float
            mean of distribution of viscosity parameter value
        sigma : float
            standard deviation of distribution of viscosity parameter value
        slip_length : float
            slip length value
        hyperbolic : boolean
            true if hyperbolic, false if not hyperbolic
        """

        self.initial_condition = initial_condition
        self.distr = distr
        self.mu = mu 
        self.sigma = sigma
        self.slip_length = slip_length
        self.hyperbolic = hyperbolic


    def compute_system_matrix(self,
                              mom_order: int,
                              SG_order: int,
                              values: np.array,
                              **kwargs) -> np.array:
        
        g = kwargs["g"] if "g" in kwargs else 1
        A = np.zeros(((mom_order + 2)*(SG_order + 1), (mom_order + 2)*(SG_order + 1))) 
        
        if mom_order == 0:
            if SG_order == 0:
                h0 = values[0]
                q0 = values[1]

                A[0][0] = 0
                A[0][1] = 1
                A[1][0] = g*h0 - (q0/h0)**2
                A[1][1] = 2.*q0/h0
            
            elif SG_order == 1:
                h0 = values[0]
                h1 = values[1]
                q0 = values[2]
                q1 = values[3]

                if self.distr == "normal" or self.distr == "uniform":
                    A[0][0] = 0
                    A[0][1] = 0
                    A[0][2] = 1
                    A[0][3] = 0
                    A[1][0] = 0
                    A[1][1] = 0
                    A[1][2] = 0
                    A[1][3] = 1
                    A[2][0] = g*h0 - (-((h1*q0)/(h0**2 - h1**2)) + (h0*q1)/(h0**2 - h1**2))**2 - ((h0*q0)/(h0**2 - h1**2) - (h1*q1)/(h0**2 - h1**2))**2
                    A[2][1] = g*h1 - 2*(-((h1*q0)/(h0**2 - h1**2)) + (h0*q1)/(h0**2 - h1**2))*((h0*q0)/(h0**2 - h1**2) - (h1*q1)/(h0**2 - h1**2))
                    A[2][2] = (2*h0*q0)/(h0**2 - h1**2) - (2*h1*q1)/(h0**2 - h1**2)
                    A[2][3] = (-2*h1*q0)/(h0**2 - h1**2) + (2*h0*q1)/(h0**2 - h1**2)
                    A[3][0] = g*h1 - 2*(-((h1*q0)/(h0**2 - h1**2)) + (h0*q1)/(h0**2 - h1**2))*((h0*q0)/(h0**2 - h1**2) - (h1*q1)/(h0**2 - h1**2))
                    A[3][1] = g*h0 - (-((h1*q0)/(h0**2 - h1**2)) + (h0*q1)/(h0**2 - h1**2))**2 - ((h0*q0)/(h0**2 - h1**2) - (h1*q1)/(h0**2 - h1**2))**2
                    A[3][2] = (-2*h1*q0)/(h0**2 - h1**2) + (2*h0*q1)/(h0**2 - h1**2)
                    A[3][3] = (2*h0*q0)/(h0**2 - h1**2) - (2*h1*q1)/(h0**2 -  h1**2)
                
                else:
                    print("This distribution is not implemented yet for mom_order=0 and SG_order=1")

            elif SG_order == 2:
                h0 = values[0]
                h1 = values[1]
                h2 = values[2]
                q0 = values[3]
                q1 = values[4]
                q2 = values[5]

                if self.distr == "normal":
                    denominator = h0**3 - 3*h0*h1**2 + 3*np.sqrt(2)*h0**2*h2 + 3*h0*h2**2 - np.sqrt(2)*h2**3
                    
                    A[0][0] = 0
                    A[0][1] = 0
                    A[0][2] = 0
                    A[0][3] = 1
                    A[0][4] = 0
                    A[0][5] = 0
                    A[1][0] = 0
                    A[1][1] = 0
                    A[1][2] = 0
                    A[1][3] = 0
                    A[1][4] = 1
                    A[1][5] = 0
                    A[2][0] = 0
                    A[2][1] = 0
                    A[2][2] = 0
                    A[2][3] = 0
                    A[2][4] = 0
                    A[2][5] = 1
                    A[3][0] = g*h0                                                                                      \
                            - (((np.sqrt(2)*h1**2 - h0*h2 -   np.sqrt(2)*h2**2)                 *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)**2    \
                            - (((-(h0*h1) - np.sqrt(2)*h1*h2)                                   *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)**2    \
                            - (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator)**2
                    A[3][1] = g*h1                                                                                      \
                            - np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)        *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)       \
                            * (((-(h0*h1) - np.sqrt(2)*h1*h2)                                   *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            - (((-(h0*h1) - np.sqrt(2)*h1*h2)                                   *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator)       \
                            - (((-(h0*h1) - np.sqrt(2)*h1*h2)                                   *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            *(((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                 *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)        *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator))
                    A[3][2] = g*h2                                                                                      \
                            - np.sqrt(2)*(((-(h0*h1) - np.sqrt(2)*h1*h2)                        *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)**2    \
                            - (((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                   *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator)       \
                            - (((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                   *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + 2*np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)      *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator))
                    A[3][3] = (2*(h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)               *q0)/denominator        \
                            + (2*(-(h0*h1) - np.sqrt(2)*h1*h2)                                  *q1)/denominator        \
                            + (2*(np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                  *q2)/denominator
                    A[3][4] = (2*(-(h0*h1) - np.sqrt(2)*h1*h2)                                  *q0)/denominator        \
                            + (2*(h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                           *q1)/denominator        \
                            + (2*(-(np.sqrt(2)*h0*h1) + h1*h2)                                  *q2)/denominator
                    A[3][5] = (2*(np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                  *q0)/denominator        \
                            + (2*(-(np.sqrt(2)*h0*h1) + h1*h2)                                  *q1)/denominator        \
                            + (2*(h0**2 - h1**2 + np.sqrt(2)*h0*h2)                             *q2)/denominator
                    A[4][0] = g*h1                                                                                      \
                            - (((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                   *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)       \
                            * ((np.sqrt(2)*(h0**2 - h1**2 + np.sqrt(2)*h0*h2)                   *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                  *(q0 + np.sqrt(2)*q2))/denominator)       \
                            - (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q1)/denominator        \
                            + (np.sqrt(2)*(np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)         *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                  *(q0 + np.sqrt(2)*q2))/denominator)       \
                            - (((-(h0*h1) - np.sqrt(2)*h1*h2)                                   *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            * ((np.sqrt(2)*(-(np.sqrt(2)*h0*h1) + h1*h2)                        *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)           *(q0 + np.sqrt(2)*q2))/denominator)
                    A[4][1] = g*(h0 + np.sqrt(2)*h2)                                                                    \
                            - np.sqrt(2)*(((-(h0*h1) - np.sqrt(2)*h1*h2)                        *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            *((np.sqrt(2)*(h0**2 - h1**2 + np.sqrt(2)*h0*h2)                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                  *(q0 + np.sqrt(2)*q2))/denominator)       \
                            - (((-(h0*h1) - np.sqrt(2)*h1*h2)                                   *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q1)/denominator        \
                            + (np.sqrt(2)*(np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)         *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                  *(q0 + np.sqrt(2)*q2))/denominator)       \
                            - ((np.sqrt(2)*(-(np.sqrt(2)*h0*h1) + h1*h2)                        *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)           *(q0 + np.sqrt(2)*q2))/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)        *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator))
                    A[4][2] = np.sqrt(2)*g*h1                                                                           \
                            - (((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                   *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q1)/denominator        \
                            + (np.sqrt(2)*(np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)         *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                  *(q0 + np.sqrt(2)*q2))/denominator)       \
                            - np.sqrt(2)*(((-(h0*h1) - np.sqrt(2)*h1*h2)                        *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            * ((np.sqrt(2)*(-(np.sqrt(2)*h0*h1) + h1*h2)                        *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)           *(q0 + np.sqrt(2)*q2))/denominator)       \
                            - ((np.sqrt(2)*(h0**2 - h1**2 + np.sqrt(2)*h0*h2)                   *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                  *(q0 + np.sqrt(2)*q2))/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + 2*np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)      *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator))
                    A[4][3] = ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                 *q1)/denominator        \
                            + (np.sqrt(2)*(np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)         *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                  *(q0 + np.sqrt(2)*q2))/denominator
                    A[4][4] = ((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                 *q0)/denominator        \
                            + (np.sqrt(2)*(-(np.sqrt(2)*h0*h1) + h1*h2)                         *q1)/denominator        \
                            + (2*(-(h0*h1) - np.sqrt(2)*h1*h2)                                  *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)           *(q0 + np.sqrt(2)*q2))/denominator        \
                            + np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)        *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)
                    A[4][5] = (np.sqrt(2)*(h0**2 - h1**2 + np.sqrt(2)*h0*h2)                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                  *(q0 + np.sqrt(2)*q2))/denominator        \
                            + np.sqrt(2)*(((-(h0*h1) - np.sqrt(2)*h1*h2)                        *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)
                    A[5][0] = g*h2                                                                                      \
                            - (((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                   *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)       \
                            * ((np.sqrt(2)*(-(np.sqrt(2)*h0*h1) + h1*h2)                        *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)           *(q0 + 2*np.sqrt(2)*q2))/denominator)       \
                            - (((-(h0*h1) - np.sqrt(2)*h1*h2)                                   *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            * ((np.sqrt(2)*(h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                 *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q2)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                *(q0 + 2*np.sqrt(2)*q2))/denominator)       \
                            - (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator)       \
                            * ((np.sqrt(2)*(-(h0*h1) - np.sqrt(2)*h1*h2)                        *q1)/denominator        \
                            + ((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                 *q2)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)*(q0 + 2*np.sqrt(2)*q2))/denominator)
                    A[5][1] = np.sqrt(2)*g*h1                                                                           \
                            - np.sqrt(2)*(((-(h0*h1) - np.sqrt(2)*h1*h2)                        *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            * ((np.sqrt(2)*(-(np.sqrt(2)*h0*h1) + h1*h2)                        *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)           *(q0 + 2*np.sqrt(2)*q2))/denominator)       \
                            - (((-(h0*h1) - np.sqrt(2)*h1*h2)                                   *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            * ((np.sqrt(2)*(-(h0*h1) - np.sqrt(2)*h1*h2)                        *q1)/denominator        \
                            + ((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                 *q2)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)*(q0 + 2*np.sqrt(2)*q2))/denominator)       \
                            - ((np.sqrt(2)*(h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                 *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q2)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                *(q0 + 2*np.sqrt(2)*q2))/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)        *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator))
                    A[5][2] = g*(h0 + 2*np.sqrt(2)*h2)                                                                  \
                            - np.sqrt(2)*(((-(h0*h1) - np.sqrt(2)*h1*h2)                        *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)       \
                            * ((np.sqrt(2)*(h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                 *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q2)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                *(q0 + 2*np.sqrt(2)*q2))/denominator)       \
                            - (((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                   *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)       \
                            * ((np.sqrt(2)*(-(h0*h1) - np.sqrt(2)*h1*h2)                        *q1)/denominator        \
                            + ((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                 *q2)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)*(q0 + 2*np.sqrt(2)*q2))/denominator)       \
                            - ((np.sqrt(2)*(-(np.sqrt(2)*h0*h1) + h1*h2)                        *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)           *(q0 + 2*np.sqrt(2)*q2))/denominator)       \
                            * (((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                *q0)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q2)/denominator        \
                            + 2*np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)      *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator))
                    A[5][3] = ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                    *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + (np.sqrt(2)*(-(h0*h1) - np.sqrt(2)*h1*h2)                         *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator        \
                            + ((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                 *q2)/denominator        \
                            + ((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)*(q0 + 2*np.sqrt(2)*q2))/denominator
                    A[5][4] = (np.sqrt(2)*(h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                  *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q2)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                *(q0 + 2*np.sqrt(2)*q2))/denominator        \
                            + np.sqrt(2)*(((-(h0*h1) - np.sqrt(2)*h1*h2)                        *q0)/denominator        \
                            + ((h0**2 + 2*np.sqrt(2)*h0*h2 - h2**2)                             *q1)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q2)/denominator)
                    A[5][5] = ((h0**2 - 2*h1**2 + 3*np.sqrt(2)*h0*h2 + 4*h2**2)                 *q0)/denominator        \
                            + (np.sqrt(2)*(-(np.sqrt(2)*h0*h1) + h1*h2)                         *q1)/denominator        \
                            + ((-(h0*h1) - np.sqrt(2)*h1*h2)                                    *q1)/denominator        \
                            + (2*(np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)                  *q2)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)           *(q0 + 2*np.sqrt(2)*q2))/denominator        \
                            + 2*np.sqrt(2)*(((np.sqrt(2)*h1**2 - h0*h2 - np.sqrt(2)*h2**2)      *q0)/denominator        \
                            + ((-(np.sqrt(2)*h0*h1) + h1*h2)                                    *q1)/denominator        \
                            + ((h0**2 - h1**2 + np.sqrt(2)*h0*h2)                               *q2)/denominator)

                elif self.distr == "uniform":
                    denominator = h0**3 - (9*h0*h1**2)/5. + (24*h0**2*h2)/(7.*np.sqrt(5)) + (18*h1**2*h2)/(7.*np.sqrt(5)) - (3*h0*h2**2)/7. - (2*h2**3)/np.sqrt(5)
                    
                    A[0][0] = 0
                    A[0][1] = 0
                    A[0][2] = 0
                    A[0][3] = 1
                    A[0][4] = 0
                    A[0][5] = 0
                    A[1][0] = 0
                    A[1][1] = 0
                    A[1][2] = 0
                    A[1][3] = 0
                    A[1][4] = 1
                    A[1][5] = 0
                    A[2][0] = 0
                    A[2][1] = 0
                    A[2][2] = 0
                    A[2][3] = 0
                    A[2][4] = 0
                    A[2][5] = 1
                    A[3][0] = g*h0                                                                                                              \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)**2            \
                            - (((-h0*h1 + (4*h1*h2)/(7.*np.sqrt(5)))                                            *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)**2            \
                            - ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q0)/denominator                \
                            + ((-h0*h1 + (4*h1*h2)/(7.*np.sqrt(5)))                                             *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*(h2**2)/np.sqrt(5))                           *q2)/denominator)**2
                    A[3][1] = g*h1                                                                                                              \
                            - (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)               \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                \
                            + ((h0**2 +  (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                       *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator))/ np.sqrt(5)  \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)               \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))
                    A[3][2] = g*h2                                                                                                              \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)**2)/np.sqrt(5)\
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)               \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)
                    A[3][3] = (2*(h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)             *q0)/denominator                \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q1)/denominator                \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q2)/denominator                 
                    A[3][4] = (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q0)/denominator                \
                            + (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                      *q1)/denominator                \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *q2)/denominator                 
                    A[3][5] = (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q0)/denominator                \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *q1)/denominator                \
                            + (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                         *q2)/denominator
                    A[4][0] = g*h1                                                                                                              \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)               \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                        *q1)/(np.sqrt(5)*denominator)   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator)               \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q1)/denominator                \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)               \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                               *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator)        
                    A[4][1] = g*(h0 + (2*h2)/np.sqrt(5))                                                                                        \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                        *q1)/(np.sqrt(5)*denominator)   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator))/np.sqrt(5)   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q1)/denominator                \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)               \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                               *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))
                    A[4][2] = (2*g*h1)/np.sqrt(5)                                                                                               \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q1)/denominator                \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)               \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                               *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator))/np.sqrt(5)   \
                            - ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                        *q1)/(np.sqrt(5)*denominator)   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)
                    A[4][3] = ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q1)/denominator                \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q1)/(np.sqrt(5)*denominator)   \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator
                    A[4][4] = ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q0)/denominator                \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *q1)/(np.sqrt(5)*denominator)   \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator                \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5)
                    A[4][5] = (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                         *q1)/(np.sqrt(5)*denominator)   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator                \
                            + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator))/np.sqrt(5)
                    A[5][0] = g*h2                                                                                                              \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)               \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                               *q1)/(np.sqrt(5)*denominator)   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)               \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                     *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)               \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)               \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                        *q1)/(np.sqrt(5)*denominator)   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)
                    A[5][1] = (2*g*h1)/np.sqrt(5)                                                                                               \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                               *q1)/(np.sqrt(5)*denominator)   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator))/np.sqrt(5)   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                        *q1)/(np.sqrt(5)*denominator)   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)               \
                            - ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                     *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))  
                    A[5][2] = g*(h0 + (2*np.sqrt(5)*h2)/7.)                                                                                     \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)               \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                     *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator))/np.sqrt(5)   \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)               \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                        *q1)/(np.sqrt(5)*denominator)   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)               \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                               *q1)/(np.sqrt(5)*denominator)   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)               \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)
                    A[5][3] = (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q1)/(np.sqrt(5)*denominator)   \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator                \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator
                    A[5][4] = (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                      *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator                \
                            + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator))/np.sqrt(5)
                    A[5][5] = ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q0)/denominator                \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *q1)/(np.sqrt(5)*denominator)   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q2)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator                \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.
                
                else:
                    print("This distribution is not implemented yet for mom_order=0 and SG_order=2")
            
            else:
                print("This stochastic Galerkin order is not implemented yet for mom_order=0")

        elif mom_order == 1:
            if SG_order == 0:
                h0 = values[0]
                q0 = values[1]
                r0 = values[2]

                A[0][0] = 0
                A[0][1] = 1
                A[0][2] = 0
                A[1][0] = g*h0 - (q0/h0)**2 - 1/3*(r0/h0)**2
                A[1][1] = 2.*q0/h0
                A[1][2] = 2/3*r0/h0
                A[2][0] = -2.*(q0*r0)/h0**2
                A[2][1] = 2.*r0/h0     
                A[2][2] = q0/h0

            elif SG_order == 1:
                if self.distr == "normal" or self.distr == "uniform":
                    h0 = values[0]
                    h1 = values[1]
                    q0 = values[2]
                    q1 = values[3]
                    r0 = values[4]
                    r1 = values[5]

                    '''
                    B = np.zeros(A.shape)

                    h_plus = h0 + h1
                    h_min  = h0 - h1
                    q_plus = q0 + q1
                    q_min  = q0 - q1
                    r_plus = r0 + r1
                    r_min  = r0 - r1

                    u_plus = q_plus/h_plus
                    u_min  = q_min/h_min
                    a_plus = r_plus/h_plus
                    a_min  = r_min/h_min

                    B[0][0] = 0
                    B[0][1] = 0
                    B[0][2] = 1
                    B[0][3] = 0
                    B[0][4] = 0
                    B[0][5] = 0

                    B[1][0] = 0
                    B[1][1] = 0
                    B[1][2] = 0
                    B[1][3] = 1
                    B[1][4] = 0
                    B[1][5] = 0

                    B[2][0] = -a_min**2/6 + g*h_min/2 - u_min**2/2 - a_plus**2/6 + g*h_plus/2 - u_plus**2/2
                    B[2][1] =  a_min**2/6 - g*h_min/2 + u_min**2/2 - a_plus**2/6 + g*h_plus/2 - u_plus**2/2
                    B[2][2] =  u_min + u_plus
                    B[2][3] = -u_min + u_plus
                    B[2][4] = (a_min  + a_plus)/3
                    B[2][5] = (-a_min + a_plus)/3

                    B[3][0] =  a_min**2/6 - g*h_min/2 + u_min**2/2 - a_plus**2/6 + g*h_plus/2 - u_plus**2/2
                    B[3][1] = -a_min**2/6 + g*h_min/2 - u_min**2/2 - a_plus**2/6 + g*h_plus/2 - u_plus**2/2
                    B[3][2] = -u_min + u_plus
                    B[3][3] =  u_min + u_plus
                    B[3][4] = (-a_min + a_plus)/3
                    B[3][5] = (a_min  + a_plus)/3

                    B[4][0] = -a_min*u_min - a_plus*u_plus
                    B[4][1] =  a_min*u_min - a_plus*u_plus
                    B[4][2] =  a_min + a_plus
                    B[4][3] = -a_min + a_plus
                    B[4][4] = (u_min  + u_plus)/2
                    B[4][5] = (-u_min + u_plus)/2

                    B[5][0] =  a_min*u_min - a_plus*u_plus
                    B[5][1] = -a_min*u_min - a_plus*u_plus
                    B[5][2] = -a_min + a_plus
                    B[5][3] =  a_min + a_plus
                    B[5][4] = (-u_min + u_plus)/2
                    B[5][5] = (u_min + u_plus)/2
                    '''

                    A[0][0] = 0
                    A[0][1] = 0
                    A[0][2] = 1
                    A[0][3] = 0
                    A[0][4] = 0
                    A[0][5] = 0

                    A[1][0] = 0
                    A[1][1] = 0
                    A[1][2] = 0
                    A[1][3] = 1
                    A[1][4] = 0
                    A[1][5] = 0

                    A[2][0] = -1*(-3*g*h0**5 + 6*g*h0**3*h1**2 - h0*h1*(3*g*h1**3 + 12*q0*q1 + 4*r0*r1) + h0**2*(3*q0**2 + 3*q1**2 + r0**2 + r1**2) + h1**2*(3*q0**2 + 3*q1**2 + r0**2 + r1**2))/(3*(h0**2 - h1**2)**2)
                    A[2][1] = (3*g*h0**4*h1 + h1**2*(3*g*h1**3 - 6*q0*q1 - 2*r0*r1) - 2*h0**2*(3*g*h1**3 + 3*q0*q1 + r0*r1) + 2*h0*h1*(3*q0**2 + 3*q1**2 + r0**2 + r1**2))/(3*(h0**2 - h1**2)**2)
                    A[2][2] = (2*(h0*q0 - h1*q1))/(h0**2 - h1**2)
                    A[2][3] = (2*(h0*q1 - h1*q0))/(h0**2 - h1**2)
                    A[2][4] = (2*(h0*r0 - h1*r1))/(3*(h0**2 - h1**2))
                    A[2][5] = (2*(h0*r1 - h1*r0))/(3*(h0**2 - h1**2))

                    A[3][0] = (3*g*h0**4*h1 + h1**2*(3*g*h1**3 - 6*q0*q1 - 2*r0*r1) - 2*h0**2*(3*g*h1**3 + 3*q0*q1 + r0*r1) + 2*h0*h1*(3*q0**2 + 3*q1**2 + r0**2 + r1**2))/(3*(h0**2 - h1**2)**2)
                    A[3][1] = -1*(-3*g*h0**5 + 6*g*h0**3*h1**2 - h0*h1*(3*g*h1**3 + 12*q0*q1 + 4*r0*r1) + h0**2*(3*q0**2 + 3*q1**2 + r0**2 + r1**2) + h1**2*(3*q0**2 + 3*q1**2 + r0**2 + r1**2))/(3*(h0**2 - h1**2)**2)
                    A[3][2] = (2*(h0*q1 - h1*q0))/(h0**2 - h1**2)
                    A[3][3] = (2*(h0*q0 - h1*q1))/(h0**2 - h1**2)
                    A[3][4] = (2*(h0*r1 - h1*r0))/(3*(h0**2 - h1**2))
                    A[3][5] = (2*(h0*r0 - h1*r1))/(3*(h0**2 - h1**2))
                    
                    A[4][0] = (-2*((h0**2 + h1**2)*q0 - 2*h0*h1*q1)*r0 - 2*((h0**2 + h1**2)*q1 - 2*h0*h1*q0)*r1)/((h0**2 - h1**2)**2)
                    A[4][1] = (-2*((h0**2 + h1**2)*q1 - 2*h0*h1*q0)*r0 - 2*((h0**2 + h1**2)*q0 - 2*h0*h1*q1)*r1)/((h0**2 - h1**2)**2)
                    A[4][2] = (2*(h0*r0 - h1*r1))/(h0**2 - h1**2)
                    A[4][3] = (2*(h0*r1 - h1*r0))/(h0**2 - h1**2)
                    A[4][4] = (h0*q0 - h1*q1)/(h0**2 - h1**2)
                    A[4][5] = (h0*q1 - h1*q0)/(h0**2 - h1**2)
                    
                    A[5][0] = (-2*((h0**2 + h1**2)*q1 - 2*h0*h1*q0)*r0 - 2*((h0**2 + h1**2)*q0 - 2*h0*h1*q1)*r1)/((h0**2 - h1**2)**2)
                    A[5][1] = (-2*((h0**2 + h1**2)*q0 - 2*h0*h1*q1)*r0 - 2*((h0**2 + h1**2)*q1 - 2*h0*h1*q0)*r1)/((h0**2 - h1**2)**2)
                    A[5][2] = (2*(h0*r1 - h1*r0))/(h0**2 - h1**2)
                    A[5][3] = (2*(h0*r0 - h1*r1))/(h0**2 - h1**2)
                    A[5][4] = (h0*q1 - h1*q0)/(h0**2 - h1**2)
                    A[5][5] = (h0*q0 - h1*q1)/(h0**2 - h1**2)
                
                else:
                    print("This distribution is not implemented yet for mom_order=1 and SG_order=1")
            
            elif SG_order == 2:
                if self.hyperbolic == False:
                    h0 = values[0]
                    h1 = values[1]
                    h2 = values[2]
                    q0 = values[3]
                    q1 = values[4]
                    q2 = values[5]
                    r0 = values[6]
                    r1 = values[7]
                    r2 = values[8]
                
                else:
                    h0 = values[0]
                    h1 = values[1]
                    h2 = 0
                    q0 = values[3]
                    q1 = values[4]
                    q2 = 0
                    r0 = values[6]
                    r1 = values[7]
                    r2 = 0

                if self.distr == "normal":
                    print("This distribution is not implemented yet for mom_order=1 and SG_order=2")
                
                elif self.distr == "uniform":
                    denominator = h0**3 - (9*h0*h1**2)/5. + (24*h0**2*h2)/(7.*np.sqrt(5)) + (18*h1**2*h2)/(7.*np.sqrt(5)) - (3*h0*h2**2)/7. - (2*h2**3)/np.sqrt(5)
                    
                    A[0][0] = 0
                    A[0][1] = 0
                    A[0][2] = 0
                    A[0][3] = 1
                    A[0][4] = 0
                    A[0][5] = 0
                    A[0][6] = 0
                    A[0][7] = 0
                    A[0][8] = 0
                    A[1][0] = 0
                    A[1][1] = 0
                    A[1][2] = 0
                    A[1][3] = 0
                    A[1][4] = 1
                    A[1][5] = 0
                    A[1][6] = 0
                    A[1][7] = 0
                    A[1][8] = 0
                    A[2][0] = 0
                    A[2][1] = 0
                    A[2][2] = 0
                    A[2][3] = 0
                    A[2][4] = 0
                    A[2][5] = 1
                    A[2][6] = 0
                    A[2][7] = 0
                    A[2][8] = 0
                    A[3][0] = g*h0                                                                                                                  \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)**2                \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)**2                \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)**2                \
                            + (-((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                         *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)**2                \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)**2                \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)**2)/3.
                    A[3][1] = g*h1                                                                                                                  \
                            - (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator))/np.sqrt(5)       \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))      \
                            + ((-2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                      *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)       \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5)))/3.
                    A[3][2] = g*h2                                                                                                                  \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)**2)/np.sqrt(5)    \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)              \
                            + ((-2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                     *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)**2)/np.sqrt(5)    \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)                   \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/7.))/3.
                    A[3][3] = (2*(h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)             *q0)/denominator                    \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q2)/denominator                   
                    A[3][4] = (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q0)/denominator                    \
                            + (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                      *q1)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *q2)/denominator
                    A[3][5] = (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *q1)/denominator                    \
                            + (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                         *q2)/denominator
                    A[3][6] = ((2*(h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)            *r0)/denominator                    \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *r2)/denominator)/3.
                    A[3][7] = ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                        *r0)/denominator                    \
                            + (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                      *r1)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *r2)/denominator)/3.
                    A[3][8] = ((2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                         *r0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *r1)/denominator                    \
                            + (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                         *r2)/denominator)/3.
                    A[4][0] = g*h1                                                                                                                  \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            + (-(((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(r0 + (2*r2)/np.sqrt(5)))/denominator))                  \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(r0 + (2*r2)/np.sqrt(5)))/denominator))/3.
                    A[4][1] = g*(h0 + (2*h2)/np.sqrt(5))                                                                                            \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator))/np.sqrt(5)       \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))      \
                            + ((-2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                     *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(r0 + (2*r2)/np.sqrt(5)))/denominator))/np.sqrt(5)       \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5)))/3.
                    A[4][2] = (2*g*h1)/np.sqrt(5)                                                                                                   \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator))/np.sqrt(5)       \
                            - ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)              \
                            + (-(((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(r0 + (2*r2)/np.sqrt(5)))/denominator))                  \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(r0 + (2*r2)/np.sqrt(5)))/denominator))/np.sqrt(5)       \
                            - ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/7.))/3.
                    A[4][3] = ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator
                    A[4][4] = ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *q1)/(np.sqrt(5)*denominator)                   \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5)
                    A[4][5] = (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                             *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator                    \
                            + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator))/np.sqrt(5)
                    A[4][6] = (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *r1)/(np.sqrt(5)*denominator)                   \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(r0 + (2*r2)/np.sqrt(5)))/denominator)/3.
                    A[4][7] = (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *r1)/(np.sqrt(5)*denominator)                   \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(r0 + (2*r2)/np.sqrt(5)))/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5))/3.
                    A[4][8] = ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(r0 + (2*r2)/np.sqrt(5)))/denominator                    \
                            + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5))/3.
                    A[5][0] = g*h2                                                                                                                  \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            + (-(((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            *((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator))                  \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator))/3.
                    A[5][1] = (2*g*h1)/np.sqrt(5)                                                                                                   \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator))/np.sqrt(5)       \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            - ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))      \
                            + ((-2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                     *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator))/np.sqrt(5)       \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5)))/3.
                    A[5][2] = g*(h0 + (2*np.sqrt(5)*h2)/7.)                                                                                         \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator))/np.sqrt(5)       \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)              \
                            + ((-2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                     *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator))/np.sqrt(5)       \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/7.))/3.
                    A[5][3] = (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                             *q1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator                    \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator
                    A[5][4] = (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                          *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator                    \
                            + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator))/np.sqrt(5)
                    A[5][5] = ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.
                    A[5][6] = ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                             *r1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator                    \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)/3.
                    A[5][7] = ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator                    \
                            + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5))/3.
                    A[5][8] = (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *r2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/7.)/3.
                    A[6][0] = -2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            - 2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                        *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            - 2*(((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)            *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)
                    A[6][1] = (-2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                      *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5)       \
                            - (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)       \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))      \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5))
                    A[6][2] = -((((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)            *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            * ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))                  \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)              \
                            * ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator)                   \
                            - (4*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)       \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)                   \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/7.)
                    A[6][3] = (2*(h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)             *r0)/denominator                    \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *r2)/denominator
                    A[6][4] = (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *r0)/denominator                    \
                            + (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                      *r1)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *r2)/denominator
                    A[6][5] = (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *r0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *r1)/denominator                    \
                            + (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                         *r2)/denominator
                    #A[6][6] = ((2*(h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)            *q0)/denominator                    \
                    #        + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q1)/denominator                    \
                    #        + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q2)/denominator)/2.
                    A[6][6] = ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator
                    #A[6][7] = ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                        *q0)/denominator                    \
                    #        + (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                      *q1)/denominator                    \
                    #        + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *q2)/denominator)/2.
                    A[6][7] = ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator
                    #A[6][8] = ((2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                         *q0)/denominator                    \
                    #        + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                *q1)/denominator                    \
                    #        + (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                         *q2)/denominator)/2.
                    A[6][8] = (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator
                    A[7][0] = -(((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                          *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))                  \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)                   \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(r0 + (2*r2)/np.sqrt(5)))/denominator)
                    A[7][1] = (-2*((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                        *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)       \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(r0 + (2*r2)/np.sqrt(5)))/denominator))/np.sqrt(5)       \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))      \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5)) 
                    A[7][2] = -((((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)            *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))                  \
                            - (2*((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)       \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)              \
                            * ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(r0 + (2*r2)/np.sqrt(5)))/denominator)                   \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(r0 + (2*r2)/np.sqrt(5)))/denominator))/np.sqrt(5)       \
                            - ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/7.)   
                    A[7][3] = ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *r1)/(np.sqrt(5)*denominator)                   \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(r0 + (2*r2)/np.sqrt(5)))/denominator   
                    A[7][4] = ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *r1)/(np.sqrt(5)*denominator)                   \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(r0 + (2*r2)/np.sqrt(5)))/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5)
                    A[7][5] = (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                             *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(r0 + (2*r2)/np.sqrt(5)))/denominator                    \
                            + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)
                    #A[7][6] = (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                    #        + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                    #        + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q1)/denominator                    \
                    #        + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                    #        + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator                    \
                    #        + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator)/2.
                    A[7][6] = ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))              *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                     *(q0 + (2*q2)/np.sqrt(5)))/denominator
                    #A[7][7] = (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                    #        + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *q1)/(np.sqrt(5)*denominator)                   \
                    #        + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                         *q1)/denominator                    \
                    #        + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                    #        + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator                    \
                    #        + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                    #        + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                    #        + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))/2.
                    A[7][7] = (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                  *(q0 + (2*q2)/np.sqrt(5)))/denominator
                    #A[7][8] = ((2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                            *q1)/(np.sqrt(5)*denominator)                   \
                    #        + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                    #        + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator                    \
                    #        + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                    #        + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                    #        + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator))/np.sqrt(5))/2.
                    A[7][8] = (2*(h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                             *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                            *(q0 + (2*q2)/np.sqrt(5)))/denominator
                    A[8][0] = -(((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                 *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))                  \
                            - ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            - ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator)                   \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)
                    A[8][1] = (-2*((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                               *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)       \
                            - ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *q1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator)                   \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator))/np.sqrt(5)       \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/np.sqrt(5))      \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                        *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/np.sqrt(5))
                    A[8][2] = -(((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                          *q1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))                  \
                            - (2*((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                      *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                          *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)       \
                            - (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)              \
                            * ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *r1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator)                   \
                            * ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator))/np.sqrt(5)       \
                            - ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator)                   \
                            * ((2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                            *r1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator)                   \
                            - ((2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                   *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)                   \
                            * (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *r0)/denominator                    \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r2)/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/7.)
                    A[8][3] = (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                             *r1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator                    \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator
                    A[8][4] = (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                          *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator                    \
                            + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *r0)/denominator                    \
                            + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *r1)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r2)/denominator))/np.sqrt(5)
                    A[8][5] = ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *r0)/denominator                    \
                            + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *r1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *r1)/denominator                    \
                            + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *r2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(r0 + (2*np.sqrt(5)*r2)/7.))/denominator                    \
                            + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *r0)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *r1)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *r2)/denominator))/7.
                    #A[8][6] = ((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                           *q0)/denominator                    \
                    #        + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                    #        + (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                             *q1)/(np.sqrt(5)*denominator)                   \
                    #        + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator                    \
                    #        + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                    #        + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator)/2.
                    A[8][6] = (2*(-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                             *q1)/(np.sqrt(5)*denominator)                   \
                            + ((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)               *q2)/denominator                    \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))   *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator
                    #A[8][7] = ((2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                         *q1)/(np.sqrt(5)*denominator)                   \
                    #        + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                    #        + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator                    \
                    #        + (2*(((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                       *q0)/denominator                    \
                    #        + ((h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                                        *q1)/denominator                    \
                    #        + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q2)/denominator))/np.sqrt(5))/2.
                    A[8][7] = (2*(h0**2 + (2*np.sqrt(5)*h0*h2)/7. - h2**2)                          *q1)/(np.sqrt(5)*denominator)                   \
                            + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q2)/denominator                    \
                            + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                         *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator
                    #A[8][8] = (((h0**2 - (4*h1**2)/5. + (24*h0*h2)/(7.*np.sqrt(5)) + (4*h2**2)/7.)              *q0)/denominator                    \
                    #        + (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *q1)/(np.sqrt(5)*denominator)                   \
                    #        + ((-(h0*h1) + (4*h1*h2)/(7.*np.sqrt(5)))                                           *q1)/denominator                    \
                    #        + (2*((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                          *q2)/denominator                    \
                    #        + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator                    \
                    #        + (2*np.sqrt(5)*((((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))             *q0)/denominator                    \
                    #        + (((-2*h0*h1)/np.sqrt(5) + h1*h2)                                                  *q1)/denominator                    \
                    #        + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                                           *q2)/denominator))/7.)/2.
                    A[8][8] = (2*((-2*h0*h1)/np.sqrt(5) + h1*h2)                                    *q1)/(np.sqrt(5)*denominator)                   \
                            + (((2*h1**2)/np.sqrt(5) - h0*h2 - (2*h2**2)/np.sqrt(5))                            *q2)/denominator                    \
                            + ((h0**2 - h1**2 + (2*h0*h2)/np.sqrt(5))                  *(q0 + (2*np.sqrt(5)*q2)/7.))/denominator
            else:
                print("This stochastic Galerkin order is not implemented yet for mom_order=1")
        
        else:
            print("This moment order is not implemented yet for the SGSWME1D")
        
        return A


    def compute_source_term(self,
                            mom_order: int,
                            SG_order: int,
                            values: np.array,
                            delta_t: float,
                            **kwargs) -> np.array:
        
        slip_length = kwargs["slip_length"] if "slip_length" in kwargs else self.slip_length
        g           = kwargs["g"]           if "g"           in kwargs else 1
        mu = self.mu
        sigma = self.sigma
        S = np.zeros((mom_order+2)*(SG_order+1))

        if mom_order == 0:
            if SG_order == 0:
                h0 = values[0]
                q0 = values[1]

                S[0] = 0
                S[1] = -((mu*q0)/(slip_length*h0))
            
            elif SG_order == 1:
                h0 = values[0]
                h1 = values[1]
                q0 = values[2]
                q1 = values[3]

                if self.distr == "normal":
                    S[0] = 0
                    S[1] = 0
                    S[2] = (h1*(sigma*q0 + mu*q1) - h0*(mu*q0 + sigma*q1))/(slip_length*(h0**2 - h1**2))
                    S[3] = (-(h0*(sigma*q0 + mu*q1)) + h1*(mu*q0 + sigma*q1))/(slip_length*(h0**2 - h1**2))

                elif self.distr == "uniform":
                    S[0] = 0
                    S[1] = 0
                    S[2] = (h1*(np.sqrt(3)*sigma*q0 + 3*mu*q1) - h0*(3*mu*q0 + np.sqrt(3)*sigma*q1))/(3.*slip_length*(h0**2 - h1**2))
                    S[3] = (-(h0*(np.sqrt(3)*sigma*q0 + 3*mu*q1)) + h1*(3*mu*q0 + np.sqrt(3)*sigma*q1))/(3.*slip_length*(h0**2 - h1**2))
                
                else:
                    print("This distribution is not implemented yet for mom_order=0 and SG_order=1")
            
            elif SG_order == 2:
                h0 = values[0]
                h1 = values[1]
                h2 = values[2]
                q0 = values[3]
                q1 = values[4]
                q2 = values[5]
                
                if self.distr == "normal":
                    S[0] = 0
                    S[1] = 0
                    S[2] = 0
                    S[3] = (-(h0**2*(mu*q0 + sigma*q1)) + mu*h1**2*(2*q0 - np.sqrt(2)*q2) + h2**2*(-4*mu*q0 + sigma*q1 + np.sqrt(2)*mu*q2) \
                         + h1*h2*(np.sqrt(2)*sigma*q0 + np.sqrt(2)*mu*q1 - sigma*q2) + h0*(h2*(-3*np.sqrt(2)*mu*q0 - 2*np.sqrt(2)*sigma*q1 + mu*q2) \
                         + h1*(sigma*q0 + mu*q1 + np.sqrt(2)*sigma*q2)))/(slip_length*(h0**3 + 3*np.sqrt(2)*h0**2*h2 - np.sqrt(2)*h2**3 + 3*h0*(-h1**2 + h2**2)))
                    S[4] = (-(h0**2*(sigma*q0 + mu*q1 + np.sqrt(2)*sigma*q2)) + h0*(h1*(mu*q0 + 3*sigma*q1 + np.sqrt(2)*mu*q2) \
                         - h2*(2*np.sqrt(2)*sigma*q0 + 2*np.sqrt(2)*mu*q1 + sigma*q2)) + h2*(mu*h1*(np.sqrt(2)*q0 - q2) + h2*(-2*sigma*q0 + mu*q1 \
                         + np.sqrt(2)*sigma*q2)))/(slip_length*(h0**3 + 3*np.sqrt(2)*h0**2*h2 - np.sqrt(2)*h2**3 + 3*h0*(-h1**2 + h2**2)))
                    S[5] = ((2*sigma*h1*h2 + np.sqrt(2)*mu*(-h1**2 + h2**2))*q0 - h0**2*(np.sqrt(2)*sigma*q1 + mu*q2) + (mu*h1 \
                         - np.sqrt(2)*sigma*h2)*(-(h2*q1) + h1*q2) + h0*(h2*(mu*q0 - 4*sigma*q1 - np.sqrt(2)*mu*q2) + h1*(np.sqrt(2)*sigma*q0 + np.sqrt(2)*mu*q1 \
                         + 2*sigma*q2)))/(slip_length*(h0**3 + 3*np.sqrt(2)*h0**2*h2 - np.sqrt(2)*h2**3 + 3*h0*(-h1**2 + h2**2)))
                
                elif self.distr == "uniform":
                    S[0] = 0
                    S[1] = 0
                    S[2] = 0
                    S[3] = (-35*h0**2*(3*mu*q0 + np.sqrt(3)*sigma*q1) + 42*mu*h1**2*(2*q0 - np.sqrt(5)*q2) + h2**2*(-60*mu*q0 + 35*np.sqrt(3)*sigma*q1 + 42*np.sqrt(5)*mu*q2) \
                         - h1*h2*(4*np.sqrt(15)*sigma*q0 + 12*np.sqrt(5)*mu*q1 + 35*np.sqrt(3)*sigma*q2) + h0*(h2*(-72*np.sqrt(5)*mu*q0 - 10*np.sqrt(15)*sigma*q1 + 105*mu*q2) + 7*h1*(5*np.sqrt(3)*sigma*q0 \
                         + 15*mu*q1 + 2*np.sqrt(15)*sigma*q2)))/(3.*slip_length*(35*h0**3 + 24*np.sqrt(5)*h0**2*h2 + 2*np.sqrt(5)*h2*(9*h1**2 - 7*h2**2) - 3*h0*(21*h1**2 + 5*h2**2)))
                    S[4] = (-7*h0**2*(5*np.sqrt(3)*sigma*q0 + 15*mu*q1 + 2*np.sqrt(15)*sigma*q2) + h0*(21*h1*(5*mu*q0 + 3*np.sqrt(3)*sigma*q1 + 2*np.sqrt(5)*mu*q2) + h2*(-10*np.sqrt(15)*sigma*q0 \
                         - 30*np.sqrt(5)*mu*q1 + 7*np.sqrt(3)*sigma*q2)) + h2*(-3*h1*(4*np.sqrt(5)*mu*q0 + 6*np.sqrt(15)*sigma*q1 + 35*mu*q2) + h2*(8*np.sqrt(3)*sigma*q0 + 105*mu*q1 \
                         + 14*np.sqrt(15)*sigma*q2)))/(3.*slip_length*(35*h0**3 + 24*np.sqrt(5)*h0**2*h2 + 2*np.sqrt(5)*h2*(9*h1**2 - 7*h2**2) - 3*h0*(21*h1**2 + 5*h2**2)))
                    S[5] = -(21*mu*(h1*(-2*np.sqrt(5)*h0 + 5*h2)*q1 + h1**2*(2*np.sqrt(5)*q0 - 5*q2) + (5*h0 + 2*np.sqrt(5)*h2)*(-(h2*q0) + h0*q2)) + 2*np.sqrt(3)*sigma*(7*np.sqrt(5)*h0**2*q1 \
                         + h0*(10*h2*q1 - 7*h1*(np.sqrt(5)*q0 + 2*q2)) + h2*(-7*np.sqrt(5)*h2*q1 + h1*(4*q0 + 7*np.sqrt(5)*q2))))/(3.*slip_length*(35*h0**3 + 24*np.sqrt(5)*h0**2*h2 + 2*np.sqrt(5)*h2*(9*h1**2 \
                         - 7*h2**2) - 3*h0*(21*h1**2 + 5*h2**2)))
                
                else:
                    print("This distribution is not implemented yet for mom_order=0 and SG_order=2")
            
            else:
                print("This stochastic Galerkin order is not implemented yet for mom_order=0")
        
        elif mom_order == 1:
            if SG_order == 0:
                h0 = values[0]
                q0 = values[1]
                r0 = values[2]

                S[0] = 0
                S[1] = -((mu*(q0 + r0))/(slip_length*h0))
                S[2] = (-3*mu*(4*slip_length*r0 + h0*(q0 + r0)))/(slip_length*h0**2)

            elif SG_order == 1:
                h0 = values[0]
                h1 = values[1]
                q0 = values[2]
                q1 = values[3]
                r0 = values[4]
                r1 = values[5]

                if self.distr == "normal":
                    S[0] = 0
                    S[1] = 0
                    S[2] = (h1*(sigma*(q0 + r0) + mu*(q1 + r1)) - h0*(mu*(q0 + r0) + sigma*(q1 + r1)))/(slip_length*(h0**2 - h1**2))
                    S[3] = (-(h0*(sigma*(q0 + r0) + mu*(q1 + r1))) + h1*(mu*(q0 + r0) + sigma*(q1 + r1)))/(slip_length*(h0**2 - h1**2))
                    S[4] = (3*(-(h0**3*(mu*(q0 + r0) + sigma*(q1 + r1))) + h0**2*(-4*slip_length*(mu*r0 + sigma*r1) + h1*(sigma*(q0 + r0) \
                         + mu*(q1 + r1))) - h1**2*(4*slip_length*(mu*r0 + sigma*r1) + h1*(sigma*(q0 + r0) + mu*(q1 + r1))) + h0*h1*(8*slip_length*(sigma*r0 \
                         + mu*r1) + h1*(mu*(q0 + r0) + sigma*(q1 + r1)))))/(slip_length*(h0**2 - h1**2)**2)
                    S[5] = (3*(-(h0**3*(sigma*(q0 + r0) + mu*(q1 + r1))) + h0*h1*(8*slip_length*(mu*r0 + sigma*r1) + h1*(sigma*(q0 + r0) \
                         + mu*(q1 + r1))) + h0**2*(-4*slip_length*(sigma*r0 + mu*r1) + h1*(mu*(q0 + r0) + sigma*(q1 + r1))) \
                         - h1**2*(4*slip_length*(sigma*r0 + mu*r1) + h1*(mu*(q0 + r0) + sigma*(q1 + r1)))))/(slip_length*(h0**2 - h1**2)**2)
                
                elif self.distr == "uniform":
                    S[0] = 0
                    S[1] = 0
                    S[2] = (h1*(np.sqrt(3)*sigma*(q0 + r0) + 3*mu*(q1 + r1)) - h0*(3*mu*(q0 + r0) + np.sqrt(3)*sigma*(q1 + r1)))/(3.*slip_length*(h0**2 - h1**2))
                    S[3] = (-(h0*(np.sqrt(3)*sigma*(q0 + r0) + 3*mu*(q1 + r1))) + h1*(3*mu*(q0 + r0) + np.sqrt(3)*sigma*(q1 + r1)))/(3.*slip_length*(h0**2 - h1**2))
                    S[4] = (-(h0**3*(3*mu*(q0 + r0) + np.sqrt(3)*sigma*(q1 + r1))) + h1**2*(-4*slip_length*(3*mu*r0 + np.sqrt(3)*sigma*r1) - h1*(np.sqrt(3)*sigma*(q0 + r0) + 3*mu*(q1 + r1))) \
                         + h0**2*(-4*slip_length*(3*mu*r0 + np.sqrt(3)*sigma*r1) + h1*(np.sqrt(3)*sigma*(q0 + r0) + 3*mu*(q1 + r1))) 
                         + h0*h1*(8*slip_length*(np.sqrt(3)*sigma*r0 + 3*mu*r1) + h1*(3*mu*(q0 + r0) + np.sqrt(3)*sigma*(q1 + r1))))/(slip_length*(h0**2 - h1**2)**2)
                    S[5] = (-(h0**3*(np.sqrt(3)*sigma*(q0 + r0) + 3*mu*(q1 + r1))) + h0*h1*(8*slip_length*(3*mu*r0 + np.sqrt(3)*sigma*r1) + h1*(np.sqrt(3)*sigma*(q0 + r0) + 3*mu*(q1 + r1))) \
                         + h0**2*(-4*slip_length*(np.sqrt(3)*sigma*r0 + 3*mu*r1) + h1*(3*mu*(q0 + r0) + np.sqrt(3)*sigma*(q1 + r1))) - h1**2*(4*slip_length*(np.sqrt(3)*sigma*r0 + 3*mu*r1) + h1*(3*mu*(q0 \
                         + r0) + np.sqrt(3)*sigma*(q1 + r1))))/(slip_length*(h0**2 - h1**2)**2)
                
                else:
                    print("This distribution is not implemented yet for mom_order=1 and SG_order=1")
            
            elif SG_order == 2:
                if self.hyperbolic == False:
                    h0 = values[0]
                    h1 = values[1]
                    h2 = values[2]
                    q0 = values[3]
                    q1 = values[4]
                    q2 = values[5]
                    r0 = values[6]
                    r1 = values[7]
                    r2 = values[8]
                
                else:
                    h0 = values[0]
                    h1 = values[1]
                    h2 = 0
                    q0 = values[3]
                    q1 = values[4]
                    q2 = 0
                    r0 = values[6]
                    r1 = values[7]
                    r2 = 0

                if self.distr == "normal":
                    print("This distribution is not implemented yet for mom_order=1 and SG_order=1")

                elif self.distr == "uniform":
                    S[0] = 0
                    S[1] = 0
                    S[2] = 0
                    S[3] = (-35*h0**2*(3*mu*q0 + np.sqrt(3)*sigma*q1 + 3*mu*r0 + np.sqrt(3)*sigma*r1) + 42*mu*h1**2*(2*q0 - np.sqrt(5)*q2 + 2*r0 - np.sqrt(5)*r2) \
                         + h2**2*(-60*mu*q0 + 35*np.sqrt(3)*sigma*q1 + 42*np.sqrt(5)*mu*q2 - 60*mu*r0 + 35*np.sqrt(3)*sigma*r1 + 42*np.sqrt(5)*mu*r2) - h1*h2*(4*np.sqrt(15)*sigma*q0 \
                         + 12*np.sqrt(5)*mu*q1 + 35*np.sqrt(3)*sigma*q2 + 4*np.sqrt(15)*sigma*r0 + 12*np.sqrt(5)*mu*r1 + 35*np.sqrt(3)*sigma*r2) + h0*(-(h2*(72*np.sqrt(5)*mu*q0 \
                         + 10*np.sqrt(15)*sigma*q1 - 105*mu*q2 + 72*np.sqrt(5)*mu*r0 + 10*np.sqrt(15)*sigma*r1 - 105*mu*r2)) + 7*h1*(5*np.sqrt(3)*sigma*q0 + 15*mu*q1 + 2*np.sqrt(15)*sigma*q2 \
                         + 5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2)))/(3.*slip_length*(35*h0**3 + 24*np.sqrt(5)*h0**2*h2 + 2*np.sqrt(5)*h2*(9*h1**2 - 7*h2**2) - 3*h0*(21*h1**2 + 5*h2**2)))
                    S[4] = (-7*h0**2*(5*np.sqrt(3)*sigma*q0 + 15*mu*q1 + 2*np.sqrt(15)*sigma*q2 + 5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2) + h0*(21*h1*(5*mu*q0 + 3*np.sqrt(3)*sigma*q1 \
                         + 2*np.sqrt(5)*mu*q2 + 5*mu*r0 + 3*np.sqrt(3)*sigma*r1 + 2*np.sqrt(5)*mu*r2) - h2*(10*np.sqrt(15)*sigma*q0 + 30*np.sqrt(5)*mu*q1 - 7*np.sqrt(3)*sigma*q2 + 10*np.sqrt(15)*sigma*r0 \
                         + 30*np.sqrt(5)*mu*r1 - 7*np.sqrt(3)*sigma*r2)) + h2*(-3*h1*(4*np.sqrt(5)*mu*q0 + 6*np.sqrt(15)*sigma*q1 + 35*mu*q2 + 4*np.sqrt(5)*mu*r0 + 6*np.sqrt(15)*sigma*r1 + 35*mu*r2) \
                         + h2*(8*np.sqrt(3)*sigma*q0 + 105*mu*q1 + 14*np.sqrt(15)*sigma*q2 + 8*np.sqrt(3)*sigma*r0 + 105*mu*r1 + 14*np.sqrt(15)*sigma*r2)))/(3.*slip_length*(35*h0**3 \
                         + 24*np.sqrt(5)*h0**2*h2 + 2*np.sqrt(5)*h2*(9*h1**2 - 7*h2**2) - 3*h0*(21*h1**2 + 5*h2**2)))
                    S[5] = (14*np.sqrt(5)*h2**2*(3*mu*q0 + np.sqrt(3)*sigma*q1 + 3*mu*r0 + np.sqrt(3)*sigma*r1) - 21*mu*h1**2*(2*np.sqrt(5)*q0 - 5*q2 + 2*np.sqrt(5)*r0 - 5*r2) - 7*h0**2*(2*np.sqrt(15)*sigma*q1 \
                         + 15*mu*q2 + 2*np.sqrt(15)*sigma*r1 + 15*mu*r2) - h1*h2*(8*np.sqrt(3)*sigma*q0 + 105*mu*q1 + 14*np.sqrt(15)*sigma*q2 + 8*np.sqrt(3)*sigma*r0 + 105*mu*r1 + 14*np.sqrt(15)*sigma*r2) \
                         + h0*(h2*(105*mu*q0 - 20*np.sqrt(3)*sigma*q1 - 42*np.sqrt(5)*mu*q2 + 105*mu*r0 - 20*np.sqrt(3)*sigma*r1 - 42*np.sqrt(5)*mu*r2) + 14*h1*(np.sqrt(15)*sigma*q0 + 3*np.sqrt(5)*mu*q1 \
                         + 2*np.sqrt(3)*sigma*q2 + np.sqrt(15)*sigma*r0 + 3*np.sqrt(5)*mu*r1 + 2*np.sqrt(3)*sigma*r2)))/(3.*slip_length*(35*h0**3 + 24*np.sqrt(5)*h0**2*h2 + 2*np.sqrt(5)*h2*(9*h1**2 - 7*h2**2) \
                         - 3*h0*(21*h1**2 + 5*h2**2)))
                    S[6] = (-1225*h0**5*(3*mu*q0 + np.sqrt(3)*sigma*q1 + 3*mu*r0 + np.sqrt(3)*sigma*r1) + 35*h0**4*(-140*slip_length*(3*mu*r0 + np.sqrt(3)*sigma*r1) - h2*(144*np.sqrt(5)*mu*q0 + 34*np.sqrt(15)*sigma*q1 \
                         - 105*mu*q2 + 144*np.sqrt(5)*mu*r0 + 34*np.sqrt(15)*sigma*r1 - 105*mu*r2) + 7*h1*(5*np.sqrt(3)*sigma*q0 + 15*mu*q1 + 2*np.sqrt(15)*sigma*q2 + 5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2)) \
                         + 5*h0**3*(147*h1**2*(13*mu*q0 + 3*np.sqrt(3)*sigma*q1 - 2*np.sqrt(5)*mu*q2 + 13*mu*r0 + 3*np.sqrt(3)*sigma*r1 - 2*np.sqrt(5)*mu*r2) + h2*(-56*slip_length*(72*np.sqrt(5)*mu*r0 + 10*np.sqrt(15)*sigma*r1 \
                         - 105*mu*r2) + h2*(-1833*mu*q0 + 110*np.sqrt(3)*sigma*q1 + 798*np.sqrt(5)*mu*q2 - 1833*mu*r0 + 110*np.sqrt(3)*sigma*r1 + 798*np.sqrt(5)*mu*r2)) + 7*h1*(h2*(20*np.sqrt(15)*sigma*q0 + 60*np.sqrt(5)*mu*q1 \
                         + 13*np.sqrt(3)*sigma*q2 + 20*np.sqrt(15)*sigma*r0 + 60*np.sqrt(5)*mu*r1 + 13*np.sqrt(3)*sigma*r2) + 56*slip_length*(5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2))) - h0*(2646*mu*h1**4 \
                         * (2*q0 - np.sqrt(5)*q2 + 2*r0 - np.sqrt(5)*r2) - 63*h1**3*h2*(14*np.sqrt(15)*sigma*q0 + 42*np.sqrt(5)*mu*q1 + 55*np.sqrt(3)*sigma*q2 + 14*np.sqrt(15)*sigma*r0 + 42*np.sqrt(5)*mu*r1 + 55*np.sqrt(3)*sigma*r2) \
                         - 5*h2**3*(16*slip_length*(-291*np.sqrt(5)*mu*r0 + 35*np.sqrt(15)*sigma*r1 + 504*mu*r2) + h2*(1188*mu*q0 + 35*np.sqrt(3)*sigma*q1 - 420*np.sqrt(5)*mu*q2 + 1188*mu*r0 + 35*np.sqrt(3)*sigma*r1 - 420*np.sqrt(5)*mu*r2)) \
                         + 9*h1**2*h2*(-112*slip_length*(31*np.sqrt(5)*mu*r0 + 5*np.sqrt(15)*sigma*r1 - 40*mu*r2) + h2*(440*mu*q0 + 345*np.sqrt(3)*sigma*q1 + 14*np.sqrt(5)*mu*q2 + 440*mu*r0 + 345*np.sqrt(3)*sigma*r1 + 14*np.sqrt(5)*mu*r2)) \
                         + 5*h1*h2**2*(h2*(86*np.sqrt(15)*sigma*q0 + 258*np.sqrt(5)*mu*q1 + 91*np.sqrt(3)*sigma*q2 + 86*np.sqrt(15)*sigma*r0 + 258*np.sqrt(5)*mu*r1 + 91*np.sqrt(3)*sigma*r2) + 192*slip_length*(4*np.sqrt(3)*sigma*r0 + 12*mu*r1 \
                         + 7*np.sqrt(15)*sigma*r2))) + 2*(378*mu*h1**4*(h2*(2*np.sqrt(5)*q0 - 5*q2 + 2*np.sqrt(5)*r0 - 5*r2) + 14*slip_length*(-2*r0 + np.sqrt(5)*r2)) + 5*h2**4*(7*h2*(12*np.sqrt(5)*mu*q0 - 7*np.sqrt(15)*sigma*q1 - 42*mu*q2 \
                         + 12*np.sqrt(5)*mu*r0 - 7*np.sqrt(15)*sigma*r1 - 42*mu*r2) - 2*slip_length*(828*mu*r0 + 245*np.sqrt(3)*sigma*r1 - 168*np.sqrt(5)*mu*r2)) - 3*h1**2*h2**2*(h2*(376*np.sqrt(5)*mu*q0 - 105*np.sqrt(15)*sigma*q1 \
                         - 1120*mu*q2 + 376*np.sqrt(5)*mu*r0 - 105*np.sqrt(15)*sigma*r1 - 1120*mu*r2) + 6*slip_length*(-1000*mu*r0 + 145*np.sqrt(3)*sigma*r1 + 434*np.sqrt(5)*mu*r2)) - 9*h1**3*h2*(42*slip_length*(2*np.sqrt(15)*sigma*r0 \
                         + 6*np.sqrt(5)*mu*r1 - 5*np.sqrt(3)*sigma*r2) + 5*h2*(4*np.sqrt(3)*sigma*q0 + 12*mu*q1 + 7*np.sqrt(15)*sigma*q2 + 4*np.sqrt(3)*sigma*r0 + 12*mu*r1 + 7*np.sqrt(15)*sigma*r2)) + 5*h1*h2**3*(2*slip_length*(110*np.sqrt(15)*sigma*r0 \
                         + 330*np.sqrt(5)*mu*r1 + 301*np.sqrt(3)*sigma*r2) + 7*h2*(4*np.sqrt(3)*sigma*q0 + 12*mu*q1 + 7*np.sqrt(15)*sigma*q2 + 4*np.sqrt(3)*sigma*r0 + 12*mu*r1 + 7*np.sqrt(15)*sigma*r2))) - h0**2*(441*h1**3*(5*np.sqrt(3)*sigma*q0 \
                         + 15*mu*q1 + 2*np.sqrt(15)*sigma*q2 + 5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2) - 63*h1**2*(37*mu*h2*(2*np.sqrt(5)*q0 - 5*q2 + 2*np.sqrt(5)*r0 - 5*r2) + 140*slip_length*(mu*r0 - np.sqrt(3)*sigma*r1 \
                         - 2*np.sqrt(5)*mu*r2)) - 5*h2**2*(h2*(222*np.sqrt(5)*mu*q0 + 296*np.sqrt(15)*sigma*q1 + 693*mu*q2 + 222*np.sqrt(5)*mu*r0 + 296*np.sqrt(15)*sigma*r1 + 693*mu*r2) + 12*slip_length*(-1101*mu*r0 + 130*np.sqrt(3)*sigma*r1 \
                         + 462*np.sqrt(5)*mu*r2)) + 15*h1*h2*(-28*slip_length*(4*np.sqrt(15)*sigma*r0 + 12*np.sqrt(5)*mu*r1 - 19*np.sqrt(3)*sigma*r2) + h2*(67*np.sqrt(3)*sigma*q0 + 201*mu*q1 + 70*np.sqrt(15)*sigma*q2 + 67*np.sqrt(3)*sigma*r0 \
                         + 201*mu*r1 + 70*np.sqrt(15)*sigma*r2))))/(slip_length*(35*h0**3 + 24*np.sqrt(5)*h0**2*h2 + 2*np.sqrt(5)*h2*(9*h1**2 - 7*h2**2) - 3*h0*(21*h1**2 + 5*h2**2))**2)
                    S[7] = -((245*h0**5*(5*np.sqrt(3)*sigma*q0 + 15*mu*q1 + 2*np.sqrt(15)*sigma*q2 + 5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2) - 35*h0**4*(21*h1*(5*mu*q0 + 3*np.sqrt(3)*sigma*q1 + 2*np.sqrt(5)*mu*q2 + 5*mu*r0 \
                         + 3*np.sqrt(3)*sigma*r1 + 2*np.sqrt(5)*mu*r2) - h2*(34*np.sqrt(15)*sigma*q0 + 102*np.sqrt(5)*mu*q1 + 41*np.sqrt(3)*sigma*q2 + 34*np.sqrt(15)*sigma*r0 + 102*np.sqrt(5)*mu*r1 + 41*np.sqrt(3)*sigma*r2) \
                         - 28*slip_length*(5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2)) + h0*h2*(-189*h1**3*(14*np.sqrt(5)*mu*q0 + 12*np.sqrt(15)*sigma*q1 + 55*mu*q2 + 14*np.sqrt(5)*mu*r0 + 12*np.sqrt(15)*sigma*r1 + 55*mu*r2) \
                         + 3*h1*h2*(h2*(430*np.sqrt(5)*mu*q0 + 204*np.sqrt(15)*sigma*q1 + 455*mu*q2 + 430*np.sqrt(5)*mu*r0 + 204*np.sqrt(15)*sigma*r1 + 455*mu*r2) + 960*slip_length*(4*mu*r0 + 6*np.sqrt(3)*sigma*r1 + 7*np.sqrt(5)*mu*r2)) \
                         - h2**2*(16*slip_length*(-149*np.sqrt(15)*sigma*r0 + 525*np.sqrt(5)*mu*r1 + 350*np.sqrt(3)*sigma*r2) + h2*(580*np.sqrt(3)*sigma*q0 + 525*mu*q1 - 308*np.sqrt(15)*sigma*q2 + 580*np.sqrt(3)*sigma*r0 + 525*mu*r1 \
                         - 308*np.sqrt(15)*sigma*r2)) + 9*h1**2*(-560*slip_length*(np.sqrt(15)*sigma*r0 + 3*np.sqrt(5)*mu*r1 + 2*np.sqrt(3)*sigma*r2) + 3*h2*(52*np.sqrt(3)*sigma*q0 + 345*mu*q1 + 28*np.sqrt(15)*sigma*q2 \
                         + 52*np.sqrt(3)*sigma*r0 + 345*mu*r1 + 28*np.sqrt(15)*sigma*r2))) + 2*h2*(27*h1**3*(42*slip_length*mu*(2*np.sqrt(5)*r0 - 5*r2) + 5*h2*(4*mu*q0 + 6*np.sqrt(3)*sigma*q1 + 7*np.sqrt(5)*mu*q2 + 4*mu*r0 \
                         + 6*np.sqrt(3)*sigma*r1 + 7*np.sqrt(5)*mu*r2)) - 3*h1*h2**2*(2*slip_length*(550*np.sqrt(5)*mu*r0 + 384*np.sqrt(15)*sigma*r1 + 1505*mu*r2) + 35*h2*(4*mu*q0 + 6*np.sqrt(3)*sigma*q1 + 7*np.sqrt(5)*mu*q2 + 4*mu*r0 \
                         + 6*np.sqrt(3)*sigma*r1 + 7*np.sqrt(5)*mu*r2)) - 9*h1**2*h2*(h2*(8*np.sqrt(15)*sigma*q0 + 105*np.sqrt(5)*mu*q1 + 70*np.sqrt(3)*sigma*q2 + 8*np.sqrt(15)*sigma*r0 + 105*np.sqrt(5)*mu*r1 + 70*np.sqrt(3)*sigma*r2) \
                         + 2*slip_length*(44*np.sqrt(3)*sigma*r0 - 435*mu*r1 - 112*np.sqrt(15)*sigma*r2)) + h2**3*(7*h2*(8*np.sqrt(15)*sigma*q0 + 105*np.sqrt(5)*mu*q1 + 70*np.sqrt(3)*sigma*q2 + 8*np.sqrt(15)*sigma*r0 + 105*np.sqrt(5)*mu*r1 \
                         + 70*np.sqrt(3)*sigma*r2) + 2*slip_length*(820*np.sqrt(3)*sigma*r0 + 3675*mu*r1 + 112*np.sqrt(15)*sigma*r2))) + h0**2*(1323*h1**3*(5*mu*q0 + 3*np.sqrt(3)*sigma*q1 + 2*np.sqrt(5)*mu*q2 + 5*mu*r0 + 3*np.sqrt(3)*sigma*r1 \
                         + 2*np.sqrt(5)*mu*r2) + 9*h1*h2*(28*slip_length*(-20*np.sqrt(5)*mu*r0 + 6*np.sqrt(15)*sigma*r1 + 95*mu*r2) + 5*h2*(67*mu*q0 + 69*np.sqrt(3)*sigma*q1 + 70*np.sqrt(5)*mu*q2 + 67*mu*r0 + 69*np.sqrt(3)*sigma*r1 + 70*np.sqrt(5)*mu*r2)) \
                         + 63*h1**2*(27*np.sqrt(3)*sigma*h2*(q2 + r2) + 28*slip_length*(5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2)) - h2**2*(h2*(832*np.sqrt(15)*sigma*q0 + 4440*np.sqrt(5)*mu*q1 + 2555*np.sqrt(3)*sigma*q2 \
                         + 832*np.sqrt(15)*sigma*r0 + 4440*np.sqrt(5)*mu*r1 + 2555*np.sqrt(3)*sigma*r2) + 12*slip_length*(-295*np.sqrt(3)*sigma*r0 + 1950*mu*r1 + 476*np.sqrt(15)*sigma*r2))) - h0**3*(441*h1**2*(5*np.sqrt(3)*sigma*q0 + 15*mu*q1 \
                         + 2*np.sqrt(15)*sigma*q2 + 5*np.sqrt(3)*sigma*r0 + 15*mu*r1 + 2*np.sqrt(15)*sigma*r2) + 21*h1*(h2*(100*np.sqrt(5)*mu*q0 + 42*np.sqrt(15)*sigma*q1 + 65*mu*q2 + 100*np.sqrt(5)*mu*r0 + 42*np.sqrt(15)*sigma*r1 + 65*mu*r2) \
                         + 280*slip_length*(5*mu*r0 + 3*np.sqrt(3)*sigma*r1 + 2*np.sqrt(5)*mu*r2)) + h2*(-280*slip_length*(10*np.sqrt(15)*sigma*r0 + 30*np.sqrt(5)*mu*r1 - 7*np.sqrt(3)*sigma*r2) + h2*(-395*np.sqrt(3)*sigma*q0 + 1650*mu*q1 \
                         + 868*np.sqrt(15)*sigma*q2 - 395*np.sqrt(3)*sigma*r0 + 1650*mu*r1 + 868*np.sqrt(15)*sigma*r2))))/(slip_length*(35*h0**3 + 24*np.sqrt(5)*h0**2*h2 + 2*np.sqrt(5)*h2*(9*h1**2 - 7*h2**2) - 3*h0*(21*h1**2 + 5*h2**2))**2))
                    S[8] = (-245*h0**5*(2*np.sqrt(15)*sigma*q1 + 15*mu*q2 + 2*np.sqrt(15)*sigma*r1 + 15*mu*r2) + 35*h0**4*(-28*slip_length*(2*np.sqrt(15)*sigma*r1 + 15*mu*r2) + h2*(105*mu*q0 - 68*np.sqrt(3)*sigma*q1 - 114*np.sqrt(5)*mu*q2 \
                         + 105*mu*r0 - 68*np.sqrt(3)*sigma*r1 - 114*np.sqrt(5)*mu*r2) + 14*h1*(np.sqrt(15)*sigma*q0 + 3*np.sqrt(5)*mu*q1 + 2*np.sqrt(3)*sigma*q2 + np.sqrt(15)*sigma*r0 + 3*np.sqrt(5)*mu*r1 + 2*np.sqrt(3)*sigma*r2)) \
                         + h0**3*(-294*h1**2*(5*np.sqrt(5)*mu*q0 - 3*np.sqrt(15)*sigma*q1 - 35*mu*q2 + 5*np.sqrt(5)*mu*r0 - 3*np.sqrt(15)*sigma*r1 - 35*mu*r2) + 5*h2*(h2*(798*np.sqrt(5)*mu*q0 + 44*np.sqrt(15)*sigma*q1 - 693*mu*q2 \
                         + 798*np.sqrt(5)*mu*r0 + 44*np.sqrt(15)*sigma*r1 - 693*mu*r2) + 56*slip_length*(105*mu*r0 - 20*np.sqrt(3)*sigma*r1 - 42*np.sqrt(5)*mu*r2)) + 7*h1*(560*slip_length*(np.sqrt(15)*sigma*r0 + 3*np.sqrt(5)*mu*r1 \
                         + 2*np.sqrt(3)*sigma*r2) + h2*(200*np.sqrt(3)*sigma*q0 + 195*mu*q1 + 26*np.sqrt(15)*sigma*q2 + 200*np.sqrt(3)*sigma*r0 + 195*mu*r1 + 26*np.sqrt(15)*sigma*r2))) - h0**2*(882*h1**3*(np.sqrt(15)*sigma*q0 + 3*np.sqrt(5)*mu*q1 \
                         + 2*np.sqrt(3)*sigma*q2 + np.sqrt(15)*sigma*r0 + 3*np.sqrt(5)*mu*r1 + 2*np.sqrt(3)*sigma*r2) + 63*h1**2*(mu*h2*(185*q0 - 52*np.sqrt(5)*q2 + 185*r0 - 52*np.sqrt(5)*r2) + 56*slip_length*(5*np.sqrt(5)*mu*r0 \
                         + np.sqrt(15)*sigma*r1 - 5*mu*r2)) - 5*h2**2*(12*slip_length*(462*np.sqrt(5)*mu*r0 + 52*np.sqrt(15)*sigma*r1 - 441*mu*r2) + h2*(693*mu*q0 + 592*np.sqrt(3)*sigma*q1 + 420*np.sqrt(5)*mu*q2 + 693*mu*r0 + 592*np.sqrt(3)*sigma*r1 \
                         + 420*np.sqrt(5)*mu*r2)) + 6*h1*h2*(h2*(67*np.sqrt(15)*sigma*q0 + 525*np.sqrt(5)*mu*q1 + 350*np.sqrt(3)*sigma*q2 + 67*np.sqrt(15)*sigma*r0 + 525*np.sqrt(5)*mu*r1 + 350*np.sqrt(3)*sigma*r2) + 14*slip_length*(-40*np.sqrt(3)*sigma*r0 \
                         + 285*mu*r1 + 38*np.sqrt(15)*sigma*r2))) + h0*(1323*mu*h1**4*(2*np.sqrt(5)*q0 - 5*q2 + 2*np.sqrt(5)*r0 - 5*r2) + 63*h1**3*h2*(28*np.sqrt(3)*sigma*q0 + 165*mu*q1 + 22*np.sqrt(15)*sigma*q2 + 28*np.sqrt(3)*sigma*r0 \
                         + 165*mu*r1 + 22*np.sqrt(15)*sigma*r2) - 70*h2**3*(h2*(30*np.sqrt(5)*mu*q0 - np.sqrt(15)*sigma*q1 - 42*mu*q2 + 30*np.sqrt(5)*mu*r0 - np.sqrt(15)*sigma*r1 - 42*mu*r2) - 8*slip_length*(72*mu*r0 + 10*np.sqrt(3)*sigma*r1 \
                         - 21*np.sqrt(5)*mu*r2)) - 9*h1**2*h2*(h2*(14*np.sqrt(5)*mu*q0 + 138*np.sqrt(15)*sigma*q1 + 595*mu*q2 + 14*np.sqrt(5)*mu*r0 + 138*np.sqrt(15)*sigma*r1 + 595*mu*r2) + 560*slip_length*(8*mu*r0 - 2*np.sqrt(3)*sigma*r1 \
                         - 7*np.sqrt(5)*mu*r2)) - h1*h2**2*(192*slip_length*(8*np.sqrt(15)*sigma*r0 + 105*np.sqrt(5)*mu*r1 + 70*np.sqrt(3)*sigma*r2) + h2*(860*np.sqrt(3)*sigma*q0 + 1365*mu*q1 + 182*np.sqrt(15)*sigma*q2 + 860*np.sqrt(3)*sigma*r0 \
                         + 1365*mu*r1 + 182*np.sqrt(15)*sigma*r2))) - 2*(189*mu*h1**4*(-28*np.sqrt(5)*slip_length*r0 + 70*slip_length*r2 + 5*h2*(2*q0 - np.sqrt(5)*q2 + 2*r0 - np.sqrt(5)*r2)) + 70*h2**4*(7*h2*(3*mu*q0 + np.sqrt(3)*sigma*q1 \
                         + 3*mu*r0 + np.sqrt(3)*sigma*r1) + 2*slip_length*(-12*np.sqrt(5)*mu*r0 + 7*np.sqrt(15)*sigma*r1 + 42*mu*r2)) - 3*h1**2*h2**2*(-6*slip_length*(434*np.sqrt(5)*mu*r0 + 58*np.sqrt(15)*sigma*r1 - 245*mu*r2) + 35*h2*(32*mu*q0 \
                         + 6*np.sqrt(3)*sigma*q1 - 7*np.sqrt(5)*mu*q2 + 32*mu*r0 + 6*np.sqrt(3)*sigma*r1 - 7*np.sqrt(5)*mu*r2)) + 9*h1**3*h2*(h2*(8*np.sqrt(15)*sigma*q0 + 105*np.sqrt(5)*mu*q1 + 70*np.sqrt(3)*sigma*q2 + 8*np.sqrt(15)*sigma*r0 \
                         + 105*np.sqrt(5)*mu*r1 + 70*np.sqrt(3)*sigma*r2) + 42*slip_length*(4*np.sqrt(3)*sigma*r0 - 15*mu*r1 - 2*np.sqrt(15)*sigma*r2)) - h1*h2**3*(7*h2*(8*np.sqrt(15)*sigma*q0 + 105*np.sqrt(5)*mu*q1 + 70*np.sqrt(3)*sigma*q2 \
                         + 8*np.sqrt(15)*sigma*r0 + 105*np.sqrt(5)*mu*r1 + 70*np.sqrt(3)*sigma*r2) + 2*slip_length*(1100*np.sqrt(3)*sigma*r0 + 4515*mu*r1 + 602*np.sqrt(15)*sigma*r2))))/(slip_length*(35*h0**3 + 24*np.sqrt(5)*h0**2*h2 \
                         + 2*np.sqrt(5)*h2*(9*h1**2 - 7*h2**2) - 3*h0*(21*h1**2 + 5*h2**2))**2)
                
                else:
                    print("This distribution is not implemented yet for mom_order=1 and SG_order=1")

            else:
                print("This stochastic Galerkin order is not implemented yet for mom_order=1")
        
        else:
            print("This moment order is not implemented yet for the SGSWME1D")
        
        return S
    

    def get_initial_values(self,
                           mom_order: int,
                           SG_order: int,
                           initial_condition: str,
                           position: float) -> np.array:
        
        initial_values = np.zeros((mom_order+2)*(SG_order+1))
        
        if initial_condition == 'constantHeight_noVelocity':
            if mom_order == 0:
                if SG_order == 0:
                    initial_values[0] = 1
                    initial_values[1] = 0*initial_values[0]
                
                elif SG_order == 1:
                    initial_values[0] = 1
                    initial_values[1] = 0
                    initial_values[2] = 0*initial_values[0]
                    initial_values[3] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=0")

            elif mom_order == 1:
                if SG_order == 0:
                    initial_values[0] = 1
                    initial_values[1] = 0*initial_values[0]
                    initial_values[2] = 0*initial_values[0]
                
                elif SG_order == 1:
                    initial_values[0] = 1
                    initial_values[1] = 0
                    initial_values[2] = 0*initial_values[0]
                    initial_values[3] = 0
                    initial_values[4] = 0*initial_values[0]
                    initial_values[5] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                    initial_values[6] = 0*initial_values[0]
                    initial_values[7] = 0
                    initial_values[8] = 0

                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=1")
            
            else:
                print("This moment order is not implemented yet for the SGSWME1D")
        
        elif initial_condition == 'constantHeight_constantVelocity':
            if mom_order == 0:
                if SG_order == 0:
                    initial_values[0] = 1
                    initial_values[1] = 1*initial_values[0]

                elif SG_order == 1:
                    initial_values[0] = 1
                    initial_values[1] = 0
                    initial_values[2] = 1*initial_values[0]
                    initial_values[3] = 0

                elif SG_order == 2:
                    initial_values[0] = 1
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 1*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=0")
            
            elif mom_order == 1:
                if SG_order == 0:
                    initial_values[0] = 1
                    initial_values[1] = 1*initial_values[0]
                    initial_values[2] = 0*initial_values[0]
                
                elif SG_order == 1:
                    initial_values[0] = 1
                    initial_values[1] = 0
                    initial_values[2] = 1*initial_values[0]
                    initial_values[3] = 0
                    initial_values[4] = 0*initial_values[0]
                    initial_values[5] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 1*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                    initial_values[6] = 0*initial_values[0]
                    initial_values[7] = 0
                    initial_values[8] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=1")
            
            else:
                print("This moment order is not implemented yet for the SGSWME1D")
        
        elif initial_condition == 'linearHeight_noVelocity':
            if mom_order == 0:
                if SG_order == 0:
                    initial_values[0] = 1 + 0.1*position
                    initial_values[1] = 0*initial_values[0]
                
                elif SG_order == 1:
                    initial_values[0] = 1 + 0.1*position
                    initial_values[1] = 0
                    initial_values[2] = 0*initial_values[0]
                    initial_values[3] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1 + 0.1*position
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=0")
            
            elif mom_order == 1:
                if SG_order == 0:
                    initial_values[0] = 1 + 0.1*position
                    initial_values[1] = 0*initial_values[0]
                    initial_values[2] = 0*initial_values[0]
                
                elif SG_order == 1:
                    initial_values[0] = 1 + 0.1*position
                    initial_values[1] = 0
                    initial_values[2] = 0*initial_values[0]
                    initial_values[3] = 0
                    initial_values[4] = 0*initial_values[0]
                    initial_values[5] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1 + 0.1*position
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                    initial_values[6] = 0*initial_values[0]
                    initial_values[7] = 0
                    initial_values[8] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=1")   

            else:
                print("This moment order is not implemented yet for the SGSWME1D")         
        
        elif initial_condition == 'damBreak_noVelocity':
            x0 = 0
            if position < x0:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 2
                        initial_values[1] = 0*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 2
                        initial_values[1] = 0*initial_values[0]
                        initial_values[2] = 0*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = 0*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = 0*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")
            
            else:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0*initial_values[0]
                        initial_values[2] = 0*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = 0*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = 0*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")                   
        
        elif initial_condition == 'damBreak_constantVelocity':
            x0 = 0
            if position < x0:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 3
                        initial_values[1] = 0.25*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 3
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 3
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 3
                        initial_values[1] = 0.25*initial_values[0]
                        initial_values[2] = 0*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 3
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = 0*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 3
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = 0*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")
            
            else:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0.25*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0.25*initial_values[0]
                        initial_values[2] = 0*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = 0*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = 0*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")        

        elif initial_condition == 'lowDamBreak_linearVelocity':
            x0 = 0
            if position < x0:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 1.5
                        initial_values[1] = 0.25*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 1.5
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1.5
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 1.5
                        initial_values[1] = 0.25*initial_values[0]
                        initial_values[2] = -0.25*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 1.5
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = -0.25*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1.5
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = -0.25*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")
            
            else:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0.25*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0.25*initial_values[0]
                        initial_values[2] = -0.25*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = -0.25*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = -0.25*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")     
        
        elif initial_condition == 'highDamBreak_linearVelocity':
            x0 = 0
            if position < x0:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 5
                        initial_values[1] = 0.25*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 5
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 5
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 5
                        initial_values[1] = 0.25*initial_values[0]
                        initial_values[2] = -0.25*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 5
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = -0.25*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 5
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = -0.25*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")
            
            else:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0.25*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0.25*initial_values[0]
                        initial_values[2] = -0.25*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0.25*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = -0.25*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.25*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = -0.25*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")     

        elif initial_condition == 'symmetric_damBreak':
            x0 = -2
            x1 = 2
            if x0 < position < x1:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 2
                        initial_values[1] = 0*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 2
                        initial_values[1] = 0*initial_values[0]
                        initial_values[2] = 0*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = 0*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = 0*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")
            
            else:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 1
                        initial_values[1] = 0*initial_values[0]
                        initial_values[2] = 0*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = 0*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 1
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = 0*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")     
        
        elif initial_condition == 'linearDamBreak_noVelocity':
            x0 = -4
            x1 = 4
            if x0 < position < x1:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 2 + (position + 4)/8.0
                        initial_values[1] = 0*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 2 + (position + 4)/8.0
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 2 + (position + 4)/8.0
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 2 + (position + 4)/8.0
                        initial_values[1] = 0*initial_values[0]
                        initial_values[2] = 0*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 2 + (position + 4)/8.0
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = 0*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 2 + (position + 4)/8.0
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = 0*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")
            
            else:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 2
                        initial_values[1] = 0*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 2
                        initial_values[1] = 0*initial_values[0]
                        initial_values[2] = 0*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = 0*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 2
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = 0*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")     

        elif initial_condition == 'smoothWave_noVelocity':
            if mom_order == 0:
                if SG_order == 0:
                    initial_values[0] = 3 + 3*np.exp(-1.5*position**2)
                    initial_values[1] = 0*initial_values[0]

                elif SG_order == 1:
                    initial_values[0] = 3 + 3*np.exp(-1.5*position**2)
                    initial_values[1] = 0
                    initial_values[2] = 0*initial_values[0]
                    initial_values[3] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 3 + 3*np.exp(-1.5*position**2)
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=0")
            
            elif mom_order == 1:
                if SG_order == 0:
                    initial_values[0] = 3 + 3*np.exp(-1.5*position**2)
                    initial_values[1] = 0*initial_values[0]
                    initial_values[2] = 0*initial_values[0]
                
                elif SG_order == 1:
                    initial_values[0] = 3 + 3*np.exp(-1.5*position**2)
                    initial_values[1] = 0
                    initial_values[2] = 0*initial_values[0]
                    initial_values[3] = 0
                    initial_values[4] = 0*initial_values[0]
                    initial_values[5] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 3 + 3*np.exp(-1.5*position**2)
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                    initial_values[6] = 0*initial_values[0]
                    initial_values[7] = 0
                    initial_values[8] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=1")
            
            else:
                print("This moment order is not implemented yet for the SGSWME1D")
        
        elif initial_condition == 'smoothWave_constantVelocity':
            if mom_order == 0:
                if SG_order == 0:
                    initial_values[0] = 1 + 0.5*np.exp(-15*position**2)
                    initial_values[1] = 0.2*initial_values[0]

                elif SG_order == 1:
                    initial_values[0] = 1 + 0.5*np.exp(-15*position**2)
                    initial_values[1] = 0
                    initial_values[2] = 0.2*initial_values[0]
                    initial_values[3] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1 + 0.5*np.exp(-15*position**2)
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0.2*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=0")
            
            elif mom_order == 1:
                if SG_order == 0:
                    initial_values[0] = 1 + 0.5*np.exp(-15*position**2)
                    initial_values[1] = 0.2*initial_values[0]
                    initial_values[2] = 0*initial_values[0]
                
                elif SG_order == 1:
                    initial_values[0] = 1 + 0.5*np.exp(-15*position**2)
                    initial_values[1] = 0
                    initial_values[2] = 0.2*initial_values[0]
                    initial_values[3] = 0
                    initial_values[4] = 0*initial_values[0]
                    initial_values[5] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1 + 0.5*np.exp(-15*position**2)
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0.2*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                    initial_values[6] = 0*initial_values[0]
                    initial_values[7] = 0
                    initial_values[8] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=1")
            
            else:
                print("This moment order is not implemented yet for the SGSWME1D")
        
        elif initial_condition == 'smoothWave_linearVelocity':
            if mom_order == 0:
                if SG_order == 0:
                    initial_values[0] = 1 + np.exp(3*np.cos(np.pi*(position + 0.5)))/np.exp(4)
                    initial_values[1] = 0.25*initial_values[0]

                elif SG_order == 1:
                    initial_values[0] = 1 + np.exp(3*np.cos(np.pi*(position + 0.5)))/np.exp(4)
                    initial_values[1] = 0
                    initial_values[2] = 0.25*initial_values[0]
                    initial_values[3] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1 + np.exp(3*np.cos(np.pi*(position + 0.5)))/np.exp(4)
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0.25*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=0")
            
            elif mom_order == 1:
                if SG_order == 0:
                    initial_values[0] = 1 + np.exp(3*np.cos(np.pi*(position + 0.5)))/np.exp(4)
                    initial_values[1] = 0.25*initial_values[0]
                    initial_values[2] = -0.25*initial_values[0]

                elif SG_order == 1:
                    initial_values[0] = 1 + np.exp(3*np.cos(np.pi*(position + 0.5)))/np.exp(4)
                    initial_values[1] = 0
                    initial_values[2] = 0.25*initial_values[0]
                    initial_values[3] = 0
                    initial_values[4] = -0.25*initial_values[0]
                    initial_values[5] = 0
                
                elif SG_order == 2:
                    initial_values[0] = 1 + np.exp(3*np.cos(np.pi*(position + 0.5)))/np.exp(4)
                    initial_values[1] = 0
                    initial_values[2] = 0
                    initial_values[3] = 0.25*initial_values[0]
                    initial_values[4] = 0
                    initial_values[5] = 0
                    initial_values[6] = -0.25*initial_values[0]
                    initial_values[7] = 0
                    initial_values[8] = 0
                
                else:
                    print("This stochastic Galerkin order is not implemented yet for mom_order=1")
            
            else:
                print("This moment order is not implemented yet for the SGSWME1D") 
        
        elif initial_condition == 'smooth_plus_damBreak':
            x0 = -7
            x1 = 7
            if x0 < position:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 4
                        initial_values[1] = 0.05*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 4
                        initial_values[1] = 0
                        initial_values[2] = 0.05*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 4
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.05*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 4
                        initial_values[1] = 0.05*initial_values[0]
                        initial_values[2] = -0.01*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 4
                        initial_values[1] = 0
                        initial_values[2] = 0.05*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = -0.01*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 4
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.05*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = -0.01*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")
            
            else:
                if mom_order == 0:
                    if SG_order == 0:
                        initial_values[0] = 3 + np.exp(-1.5*(position-x1)**2)
                        initial_values[1] = 0.05*initial_values[0]

                    elif SG_order == 1:
                        initial_values[0] = 3 + np.exp(-1.5*(position-x1)**2)
                        initial_values[1] = 0
                        initial_values[2] = 0.05*initial_values[0]
                        initial_values[3] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 3 + np.exp(-1.5*(position-x1)**2)
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.05*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=0")
                
                elif mom_order == 1:
                    if SG_order == 0:
                        initial_values[0] = 3 + np.exp(-1.5*(position-x1)**2)
                        initial_values[1] = 0.05*initial_values[0]
                        initial_values[2] = -0.01*initial_values[0]
                    
                    elif SG_order == 1:
                        initial_values[0] = 3 + np.exp(-1.5*(position-x1)**2)
                        initial_values[1] = 0
                        initial_values[2] = 0.05*initial_values[0]
                        initial_values[3] = 0
                        initial_values[4] = -0.01*initial_values[0]
                        initial_values[5] = 0
                    
                    elif SG_order == 2:
                        initial_values[0] = 3 + np.exp(-1.5*(position-x1)**2)
                        initial_values[1] = 0
                        initial_values[2] = 0
                        initial_values[3] = 0.05*initial_values[0]
                        initial_values[4] = 0
                        initial_values[5] = 0
                        initial_values[6] = -0.01*initial_values[0]
                        initial_values[7] = 0
                        initial_values[8] = 0
                    
                    else:
                        print("This stochastic Galerkin order is not implemented yet for mom_order=1")
                
                else:
                    print("This moment order is not implemented yet for the SGSWME1D")     
                
        else:
            print("This initial condition is not implemented yet for the SGSWME1D")
        
        return initial_values
    
    def compute_number_of_variables(self,
                                    mom_order: int,
                                    SG_order: int) -> int:
        number_of_variables = (mom_order + 2)*(SG_order + 1)
        return int(number_of_variables)
    
    
    def compute_all_breakdown_criteria(self,
                                   values: np.array,
                                   number_of_variables: list,
                                   n,
                                   delta_x,
                                   tolerance_up_height_gradient,
                                   tolerance_down_height_gradient,
                                   tolerance_up_momentum_gradient,
                                   tolerance_down_momentum_gradient,
                                   tolerance_up_last_moment,
                                   tolerance_down_last_moment) -> np.array:
        pass

    def compute_breakdown_criterion(self,
                                   values: list,
                                   number_of_variables: int,
                                   breakdown_criterion: str,
                                   n: int) -> np.array:
        pass
    
    def compute_exp_and_var(self,
                            mom_order: int,
                            SG_order: int,
                            values: np.array,
                            primitive: bool) -> np.array:
        
        if mom_order == 0:
            if SG_order == 0:
                func1_exp = values[:,0]
                func1_var = np.zeros(len(values[:,0]))
                func2_var = np.zeros(len(values[:,1]))
                
                if primitive == True:
                    func2_exp = np.divide(values[:,1], values[:,0])
                
                else:
                    func2_exp = values[:,1]
            
            elif SG_order == 1:
                func1_exp = values[:,0]
                func1_var = np.square(values[:,1])
                
                if primitive == True:
                    func2_exp = np.zeros(len(values[:,2]))
                    func2_var = np.zeros(len(values[:,3]))
                    
                    if self.distr == "normal":
                        for i in range(len(func2_exp)):
                            func2_exp[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*(values[i,2]*1 + values[i,3]*w)/(values[i,0]*1 + values[i,1]*w), -np.infty, np.infty)[0]
                            func2_var[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*((values[i,2]*1 + values[i,3]*w)/(values[i,0]*1 + values[i,1]*w))**2, -np.infty, np.infty)[0] - func2_exp[i]**2
                    
                    elif self.distr == "uniform":
                        for i in range(len(func2_exp)):
                            func2_exp[i] = 0.5*scipy.integrate.quad(lambda w: (values[i,2]*1 + values[i,3]*np.sqrt(3)*w)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w), -1, 1)[0]
                            func2_var[i] = 0.5*scipy.integrate.quad(lambda w: ((values[i,2]*1 + values[i,3]*np.sqrt(3)*w)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w))**2, -1, 1)[0] - func2_exp[i]**2
                                
                    else:
                        print("This distribution is not implemented yet for mom_order=0 and SG_order=1")
                
                else:
                    func2_exp = values[:,2]
                    func2_var = np.square(values[:,3])
            
            elif SG_order == 2:
                func1_exp = values[:,0]
                func1_var = np.square(values[:,1]) + np.square(values[:,2])
                
                if primitive == True:
                    func2_exp = np.zeros(len(values[:,3]))
                    func2_var = np.zeros(len(values[:,4]))
                    
                    if self.distr == "normal":
                        for i in range(len(func2_exp)):
                            func2_exp[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*(values[i,3]*1 + values[i,4]*w + values[i,5]*(w**2 - 1)/np.sqrt(2))/(values[i,0]*1 + values[i,1]*w + values[i,2]*(w**2 - 1)/np.sqrt(2)), -np.infty, np.infty)[0]
                            func2_var[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*((values[i,3]*1 + values[i,4]*w + values[i,5]*(w**2 - 1)/np.sqrt(2))/(values[i,0]*1+values[i,1]*w + values[i,2]*(w**2 - 1)/np.sqrt(2)))**2, -np.infty, np.infty)[0] - func2_exp[i]**2
                    
                    elif self.distr == "uniform":
                        for i in range(len(func2_exp)):
                            func2_exp[i] = 0.5*scipy.integrate.quad(lambda w: (values[i,3]*1 + values[i,4]*np.sqrt(3)*w + values[i,5]*np.sqrt(5)*(3*w**2 - 1)/2)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w + values[i,2]*np.sqrt(5)*(3*w**2 - 1)/2), -1, 1)[0]
                            func2_var[i] = 0.5*scipy.integrate.quad(lambda w: ((values[i,3]*1 + values[i,4]*np.sqrt(3)*w + values[i,5]*np.sqrt(5)*(3*w**2 - 1)/2)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w + values[i,2]*np.sqrt(5)*(3*w**2 - 1)/2))**2, -1, 1)[0] - func2_exp[i]**2
                    
                    else:
                        print("This distribution is not implemented yet for mom_order=0 and SG_order=2")
                
                else:
                    func2_exp = values[:,3]
                    func2_var = np.square(values[:,4]) + np.square(values[:,5])
            
            else:
                print("This stochastic Galerkin order is not implemented yet for mom_order=0")
            
            func3_exp = np.zeros(len(values[:,0]))
            func3_var = np.zeros(len(values[:,0]))
            
        elif mom_order == 1:
            if SG_order == 0:
                func1_exp = values[:,0]
                func1_var = np.zeros(len(values[:,0]))
                func2_var = np.zeros(len(values[:,1]))
                func3_var = np.zeros(len(values[:,2]))
                
                if primitive == True:
                    func2_exp = np.divide(values[:,1], values[:,0])
                    func3_exp = np.divide(values[:,2], values[:,0])
                
                else:
                    func2_exp = values[:,1]
                    func3_exp = values[:,2]
            
            elif SG_order == 1:
                func1_exp = values[:,0]
                func1_var = np.square(values[:,1])
                
                if primitive == True:
                    func2_exp = np.zeros(len(values[:,2]))
                    func2_var = np.zeros(len(values[:,3]))
                    func3_exp = np.zeros(len(values[:,2]))
                    func3_var = np.zeros(len(values[:,3]))
                    
                    if self.distr == "normal":
                        for i in range(len(func2_exp)):
                            func2_exp[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*(values[i,2]*1 + values[i,3]*w)/(values[i,0]*1 + values[i,1]*w), -np.infty, np.infty)[0]
                            func2_var[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*((values[i,2]*1 + values[i,3]*w)/(values[i,0]*1 + values[i,1]*w))**2, -np.infty, np.infty)[0] - func2_exp[i]**2
                            func3_exp[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*(values[i,4]*1 + values[i,5]*w)/(values[i,0]*1 + values[i,1]*w), -np.infty, np.infty)[0]
                            func2_var[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*((values[i,4]*1 + values[i,5]*w)/(values[i,0]*1 + values[i,1]*w))**2, -np.infty, np.infty)[0] - func3_exp[i]**2                           
                    
                    elif self.distr == "uniform":
                        for i in range(len(func2_exp)):
                            func2_exp[i] = 0.5*scipy.integrate.quad(lambda w: (values[i,2]*1 + values[i,3]*np.sqrt(3)*w)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w), -1, 1)[0]
                            func2_var[i] = 0.5*scipy.integrate.quad(lambda w: ((values[i,2]*1 + values[i,3]*np.sqrt(3)*w)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w))**2, -1, 1)[0] - func2_exp[i]**2
                            func3_exp[i] = 0.5*scipy.integrate.quad(lambda w: (values[i,4]*1 + values[i,5]*np.sqrt(3)*w)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w), -1, 1)[0]
                            func3_var[i] = 0.5*scipy.integrate.quad(lambda w: ((values[i,4]*1 + values[i,5]*np.sqrt(3)*w)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w))**2, -1, 1)[0] - func3_exp[i]**2                            
                    else:
                        print("This distribution is not implemented yet for mom_order=1 and SG_order=1")
                
                else:
                    func2_exp = values[:,2]
                    func2_var = np.square(values[:,3])
                    func3_exp = values[:,4]
                    func3_var = np.square(values[:,5])
            
            elif SG_order == 2:
                func1_exp = values[:,0]
                func1_var = np.square(values[:,1]) + np.square(values[:,2])
                
                if primitive == True:
                    func2_exp = np.zeros(len(values[:,3]))
                    func2_var = np.zeros(len(values[:,4]))
                    func3_exp = np.zeros(len(values[:,3]))
                    func3_var = np.zeros(len(values[:,4]))
                    
                    if self.distr == "normal":
                        for i in range(len(func2_exp)):
                            func2_exp[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*(values[i,3]*1 + values[i,4]*w + values[i,5]*(w**2 - 1)/np.sqrt(2))/(values[i,0]*1 + values[i,1]*w + values[i,2]*(w**2 - 1)/np.sqrt(2)), -np.infty, np.infty)[0]
                            func2_var[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*((values[i,3]*1 + values[i,4]*w + values[i,5]*(w**2 - 1)/np.sqrt(2))/(values[i,0]*1+values[i,1]*w + values[i,2]*(w**2 - 1)/np.sqrt(2)))**2, -np.infty, np.infty)[0] - func2_exp[i]**2
                            func3_exp[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*(values[i,6]*1 + values[i,7]*w + values[i,8]*(w**2 - 1)/np.sqrt(2))/(values[i,0]*1 + values[i,1]*w + values[i,2]*(w**2 - 1)/np.sqrt(2)), -np.infty, np.infty)[0]
                            func3_var[i] = 1/np.sqrt(2*np.pi)*scipy.integrate.quad(lambda w: np.exp(-w**2/2)*((values[i,6]*1 + values[i,7]*w + values[i,8]*(w**2 - 1)/np.sqrt(2))/(values[i,0]*1+values[i,1]*w + values[i,2]*(w**2 - 1)/np.sqrt(2)))**2, -np.infty, np.infty)[0] - func3_exp[i]**2

                    elif self.distr == "uniform":
                        for i in range(len(func2_exp)):
                            func2_exp[i] = 0.5*scipy.integrate.quad(lambda w: (values[i,3]*1 + values[i,4]*np.sqrt(3)*w + values[i,5]*np.sqrt(5)*(3*w**2 - 1)/2)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w + values[i,2]*np.sqrt(5)*(3*w**2 - 1)/2), -1, 1)[0]
                            func2_var[i] = 0.5*scipy.integrate.quad(lambda w: ((values[i,3]*1 + values[i,4]*np.sqrt(3)*w + values[i,5]*np.sqrt(5)*(3*w**2 - 1)/2)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w + values[i,2]*np.sqrt(5)*(3*w**2 - 1)/2))**2, -1, 1)[0] - func2_exp[i]**2
                            func3_exp[i] = 0.5*scipy.integrate.quad(lambda w: (values[i,6]*1 + values[i,7]*np.sqrt(3)*w + values[i,8]*np.sqrt(5)*(3*w**2 - 1)/2)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w + values[i,2]*np.sqrt(5)*(3*w**2 - 1)/2), -1, 1)[0]
                            func3_var[i] = 0.5*scipy.integrate.quad(lambda w: ((values[i,6]*1 + values[i,7]*np.sqrt(3)*w + values[i,8]*np.sqrt(5)*(3*w**2 - 1)/2)/(values[i,0]*1 + values[i,1]*np.sqrt(3)*w + values[i,2]*np.sqrt(5)*(3*w**2 - 1)/2))**2, -1, 1)[0] - func3_exp[i]**2
                    
                    else:
                        print("This distribution is not implemented yet for mom_order=1 and SG_order=2")
                
                else:
                    func2_exp = values[:,3]
                    func2_var = np.square(values[:,4]) + np.square(values[:,5])
                    func3_exp = values[:,6]
                    func3_var = np.square(values[:,7]) + np.square(values[:,8])
            
            else:    
                print("This stochastic Galerkin order is not implemented yet for mom_order=1")
        
        else:
            print("This moment order is not implemented yet for the SGSWME1D")     
        
        return func1_exp, func1_var, func2_exp, func2_var, func3_exp, func3_var
        
    
    def compute_vertical_velocity_profile(self,
                                          mom_order: int,
                                          SG_order: int,
                                          values: np.array,
                                          z_points: np.array) -> np.array:
        """
        reconstructs the vertical velocity profile from the solution values and evaluates the velocity profile pointwise

        Parameters
        ----------
        mom_order: integer
            moment order of the model
        SG_order: integer
            stochastic Galerkin order of the model
        values: np.array (2D)
            2D numpy array containing the values of the variables in each mesh cell
        z_points: 
            the locations in vertical direction in which the velocity is computed
        
        Returns
        -------
        velocity_profile_exp: numpy 2D array
            expected lateral velocity evaluated in in each point in z_points in z-direction
        velocity_profile_var: numpy 2D array
            variance of lateral velocity evaluated in in each point in z_points in z-direction
        """
        
        velocity_profile_exp = np.zeros((len(values), len(z_points)))
        velocity_profile_var = np.zeros((len(values), len(z_points)))
        if mom_order == 0:
            _, _, um_exp, um_var, _, _ = self.compute_exp_and_var(0, SG_order, values, True)
            for i in range(len(um_exp)):
                velocity_profile_exp[i,:] = um_exp[i]*(np.ones(len(z_points)))
                velocity_profile_var[i,:] = um_var[i]*(np.ones(len(z_points)))
        elif mom_order == 1:
            _, _, um_exp, um_var, alpha1_exp, alpha1_var = self.compute_exp_and_var(1, SG_order, values, True)
            for i in range(len(um_exp)):
                velocity_profile_exp[i,:] = um_exp[i]*(np.ones(len(z_points))) + alpha1_exp[i]*(np.ones(len(z_points)) - 2*z_points)
                velocity_profile_var[i,:] = um_var[i]*(np.ones(len(z_points))) + alpha1_var[i]*(np.ones(len(z_points)) - 2*z_points)**2
   
        return velocity_profile_exp, velocity_profile_var