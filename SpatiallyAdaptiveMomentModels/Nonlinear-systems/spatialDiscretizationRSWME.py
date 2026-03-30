from abc import ABC, abstractmethod
import numpy as np
from collections.abc import Callable

#TODO: use duck typing
class SpatialDiscretization(ABC):

    """
    This interface represents a spatial discretization.

    ...

    Attributes
    ----------
    None

    
    Abstract methods
    -------
    def compute_fluctuation(self):
        computes a fluctuation between two cells
    """

    @abstractmethod
    def compute_fluctuation(self):

        """
        Documented in the child classes

        """
        pass

class PVM(SpatialDiscretization):

    """
    This abstract class represents a polynomial viscosity method.

    ...

    Attributes
    ----------
    None

    Class methods
    -------------
    def compute_fluctuation(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        computes the fluctuation between two cells with values value_left and value_right

            
    Abstract methods    
    -----------------
    def compute_viscosity(self,roe_matrix,delta_t,delta_x):
        computes the numerical viscosity matrix
    """

    @abstractmethod
    def __init__(self):
        pass

    def compute_fluctuation(self,
                            value_left: np.array,
                            value_right: np.array,
                            system_matrix: Callable[...,np.array],
                            direction: str,
                            delta_t: float,
                            delta_x: float) -> np.array:
        
        """
        Computes the fluctuations between two cells containing the values value_left and value_right

        Parameters
        ----------
        value_left: np.array
            value of the cell left of the boundary
        value_right: np.array
            value of the cell right of the boundary
        system_matrix: function
            function that computes the system matrix along the path between value_left and value_right
        direction: string 
            direction of the fluctuation: 'positive' if from left to right, 'negative' if from right to left
        delta_t: float
            time step size
        delta_x: float
            spatial discretization step size
        
        Returns
        -------
        fluctuation: np.array
            fluctuation between two cells containing the values value_left and value_right

        """
        
        generalized_roe = system_matrix((value_left+value_right)/2)
        viscosity = self.compute_viscosity(generalized_roe,delta_t,delta_x)
        if direction == 'negative':
            viscosity *= -1
        fluctuation = (np.dot(generalized_roe,value_right-value_left) + np.dot(viscosity,value_right-value_left))/2

        return fluctuation

    @abstractmethod
    def compute_viscosity(self,
                          roe_matrix: np.array,
                          delta_t: float,
                          delta_x: float):
        
        """
        Computes the numerical viscosity matrix

        Parameters
        ----------
        roe_matrix: np.array
            roe matrix at the interface
        delta_t: float
            time step size
        delta_x: float
            spatial discretization step size
        
        Returns
        -------
        viscosity: np.array
            numerical viscosity matrix

        """

        pass

class PRICE(PVM):

    def __init__(self):
        pass

    def compute_viscosity(self,
                          roe_matrix: np.array,
                          delta_t: float,
                          delta_x: float):
        #viscosity = delta_x/(2*delta_t)*np.identity(np.shape(roe_matrix)[0])+delta_t/(2*delta_x)*roe_matrix 
        viscosity = 0.5*delta_x/delta_t*np.identity(np.shape(roe_matrix)[0])+0.5*delta_t/delta_x*roe_matrix@roe_matrix 
        return viscosity
    
class LF(PVM):

    def __init__(self):
        pass

    def compute_viscosity(self,
                          roe_matrix: np.array,
                          delta_t: float,
                          delta_x: float):
        viscosity = delta_x/delta_t*np.identity(np.shape(roe_matrix)[0])
        return viscosity
