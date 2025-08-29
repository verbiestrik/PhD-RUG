from abc import ABC, abstractmethod
import numpy as np
from collections.abc import Callable

class SpatialDiscretization(ABC):

    """
    This interface represents a spatial discretization.

    ...

    Attributes
    ----------
    non_conservative : boolean
        True if the spatial discretization is of the conservative type, false if non-conservative

    
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
    def compute_generalized_roe_and_viscosity(self,value_left,value_right,system_matrix,direction,delta_t,delta_x)
        compute the generalized roe matrix and the viscosity matrix between two cells with values value_left and value_right
            
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
                            delta_x: float,
                            **kwargs) -> np.array:
        
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
        # Nodes on [0, 1]
        quadrature_nodes = [
            (1 - (1/3) * np.sqrt((5 + 2*np.sqrt(10/7))/3)) / 2,
            (1 - (1/3) * np.sqrt((5 - 2*np.sqrt(10/7))/3)) / 2,
            1/2,
            (1 + (1/3) * np.sqrt((5 - 2*np.sqrt(10/7))/3)) / 2,
            (1 + (1/3) * np.sqrt((5 + 2*np.sqrt(10/7))/3)) / 2
        ]

        # Weights on [0, 1]
        quadrature_weights = [
            (322 - 13*np.sqrt(70)) / 1800,
            (322 + 13*np.sqrt(70)) / 1800,
            128 / 450,
            (322 + 13*np.sqrt(70)) / 1800,
            (322 - 13*np.sqrt(70)) / 1800
        ]

        # Nodes on [0, 1]
        quadrature_nodes = [1/2]

        # Weights on [0, 1]
        quadrature_weights = [1]

        generalized_roe = 0
        for i in range(len(quadrature_nodes)):
            generalized_roe += quadrature_weights[i]*(system_matrix((1-quadrature_nodes[i])*value_left+(quadrature_nodes[i])*value_right,**kwargs))
        viscosity = self.compute_viscosity(generalized_roe,delta_t,delta_x)
        if direction == 'negative':
            viscosity *= -1
        fluctuation = (np.dot(generalized_roe,value_right-value_left) + np.dot(viscosity,value_right-value_left))/2

        return fluctuation
    
    def compute_generalized_roe_and_viscosity(self,value_left,value_right,system_matrix,direction,delta_t,delta_x,**kwargs):
        # Nodes on [0, 1]
        quadrature_nodes = [
            (1 - (1/3) * np.sqrt((5 + 2*np.sqrt(10/7))/3)) / 2,
            (1 - (1/3) * np.sqrt((5 - 2*np.sqrt(10/7))/3)) / 2,
            1/2,
            (1 + (1/3) * np.sqrt((5 - 2*np.sqrt(10/7))/3)) / 2,
            (1 + (1/3) * np.sqrt((5 + 2*np.sqrt(10/7))/3)) / 2
        ]

        # Weights on [0, 1]
        quadrature_weights = [
            (322 - 13*np.sqrt(70)) / 1800,
            (322 + 13*np.sqrt(70)) / 1800,
            128 / 450,
            (322 + 13*np.sqrt(70)) / 1800,
            (322 - 13*np.sqrt(70)) / 1800
        ]

        generalized_roe = 0
        for i in range(len(quadrature_nodes)):
            generalized_roe += quadrature_weights[i]*(system_matrix((1-quadrature_nodes[i])*value_left+(quadrature_nodes[i])*value_right,**kwargs))
        viscosity = self.compute_viscosity(generalized_roe,delta_t,delta_x)
        if direction == 'negative':
            viscosity *= -1
        return (generalized_roe + viscosity)/2.

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