from abc import ABC, abstractmethod
import numpy as np
from collections.abc import Callable
import mpmath as mp
from numpy.polynomial.legendre import leggauss
from functools import lru_cache

class SpatialDiscretization(ABC):

    """
    This abstract class represents a spatial discretization.

    ...

    Attributes
    ----------
    None
    

    Abstract methods
    ----------------
    def compute_fluctuation(self):
        computes a fluctuation between two cells
    def __init__(self):
        initializes the SpatialDiscretization object

    """

    @abstractmethod
    def __init__(self):
        """
        Documented in the child classes

        """
        pass

    @abstractmethod
    def compute_fluctuation(self):

        """
        Documented in the child classes

        """
        pass

class PVM(SpatialDiscretization,ABC):

    """
    This abstract class represents a polynomial viscosity method.

    ...

    Attributes
    ----------
    quadrature_nodes: np.ndarray
        the nodes of the Gauss-Legendre quadrature rule on the interval [0,1]
    quadrature_weights: np.ndarray
        the weights of the Gauss-Legendre quadrature rule on the interval [0,1]

    Methods implemented from the abstract class SpatialDiscretization
    -----------------------------------------------------------------
    def __init__(self,nr_of_quadrature_points):
        initializes the PVM object
    def compute_fluctuation(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        computes the fluctuation between two cells with values value_left and value_right
    
    Abstract methods    
    -----------------
    def compute_viscosity(self,roe_matrix,delta_t,delta_x):
        computes the numerical viscosity matrix


    Instance methods
    -------------
    def compute_generalized_roe_and_viscosity(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        compute the generalized roe matrix and the viscosity matrix between two cells with values value_left and value_right
    def _compute_quadrature_points(self,nr_of_quadrature_points):
        computes the nodes and the weights for a Gauss-Legendre rule on the interval [0,1]
        with nr_of_quadrature_points quadrature points
    """

    def __init__(self,
                 nr_of_quadrature_points: int):
        """
        Initializes the PVM object and sets the number of quadrature points.

        Parameters
        ----------
        nr_of_quadrature_points: int
            the number of quadrature points used in the numerical integration of the integral appearing in
            the generalized Roe linearization
        """
        self.quadrature_nodes, self.quadrature_weights = self._compute_quadrature_points(
            nr_of_quadrature_points)

    def compute_fluctuation(self,
                            value_left: np.ndarray,
                            value_right: np.ndarray,
                            system_matrix: Callable[...,np.ndarray],
                            delta_t: float,
                            delta_x: float) -> tuple[np.ndarray,np.ndarray]:
        
        """
        Computes the fluctuations between two cells containing the values value_left and value_right

        Parameters
        ----------
        value_left: np.ndarray
            value of the cell left of the boundary
        value_right: np.ndarray
            value of the cell right of the boundary
        system_matrix: function
            function that computes the system matrix along the path between value_left and value_right
        delta_t: float
            time step size
        delta_x: float
            spatial discretization step size
        
        Returns
        -------
        fluctuation_min: np.ndarray
            fluctuation between two cells containing the values value_left and value_right in negative direction
        fluctuation_plus: np.ndarray
            fluctuation between two cells containing the values value_left and value_right in positive direction

        """

        generalized_roe = 0
        for i in range(len(self.quadrature_nodes)):
            generalized_roe += self.quadrature_weights[i]*\
                (system_matrix((1-self.quadrature_nodes[i])*value_left+(self.quadrature_nodes[i])*value_right))
        viscosity = np.dot(self.compute_viscosity(generalized_roe,delta_t,delta_x),value_right-value_left)
        generalized_roe = np.dot(generalized_roe,value_right-value_left)
        fluctuation_min = (generalized_roe - viscosity)/2
        fluctuation_plus = (generalized_roe + viscosity)/2

        return fluctuation_min, fluctuation_plus
    
    def compute_generalized_roe_and_viscosity(self,
                                              value_left: np.ndarray,
                                              value_right: np.ndarray,
                                              system_matrix: np.ndarray,
                                              delta_t: float,
                                              delta_x: float) -> tuple[np.ndarray,np.ndarray]:

        """
        Computes the generalized roe matrix and the viscosity matrix between two cells with values value_left and value_right.
        The fluctuations in the PVM scheme can be written in the form D^+- = A^+- . (value_right - value_left).
        This function returns the matrices A^+ and A^-.

        Parameters
        ----------
        value_left: np.ndarray
            value of the cell left of the boundary
        value_right: np.ndarray
            value of the cell right of the boundary
        system_matrix: function
            function that computes the system matrix along the path between value_left and value_right
        delta_t: float
            time step size
        delta_x: float
            spatial discretization step size
        
        Returns
        -------
        fluct_matrix_min: np.ndarray
            fluctuation matrix in the negative direction
        fluct_matrix_plus: np.ndarray
            fluctuation matrix in the positive direction

        """

        generalized_roe = 0
        for i in range(len(self.quadrature_nodes)):
            generalized_roe += self.quadrature_weights[i]*\
                (system_matrix((1-self.quadrature_nodes[i])*value_left+(self.quadrature_nodes[i])*value_right))
        viscosity = self.compute_viscosity(generalized_roe,delta_t,delta_x)
        fluct_matrix_min = (generalized_roe - viscosity)/2
        fluct_matrix_plus = (generalized_roe + viscosity)/2

        return fluct_matrix_min, fluct_matrix_plus

    def _compute_quadrature_points(self,
                                   nr_of_quadrature_points: int):
        """
        Computes the nodes and the weights for a Gauss-Legendre rule on the interval [0,1]
        with nr_of_quadrature_points quadrature points.

        Parameters
        ----------
        nr_of_quadrature_points: int
            the number of quadrature points
        
        Returns
        -------
        quadrature_nodes: np.ndarray
            the nodes of the quadrature rule
        quadrature_weights: np.ndarray
            the weights of the quadrature rule

        """
        # Gauss-Legendre on [-1, 1]
        x, w = leggauss(nr_of_quadrature_points)

        # Map to [0, 1]
        quadrature_nodes = 0.5 * (x + 1.0)
        quadrature_weights = 0.5 * w

        return quadrature_nodes, quadrature_weights

    @abstractmethod
    def compute_viscosity(self,
                          roe_matrix: np.ndarray,
                          delta_t: float,
                          delta_x: float):
        
        """
        Computes the numerical viscosity matrix

        Parameters
        ----------
        roe_matrix: np.ndarray
            roe matrix at the interface
        delta_t: float
            time step size
        delta_x: float
            spatial discretization step size
        
        Returns
        -------
        viscosity: np.ndarray
            numerical viscosity matrix

        """

        pass

class PRICE(PVM):
    """
    This class represents the PRICE scheme, a PVM scheme with viscosity function Q(A) = delta_x/(2*delta_t)*I+delta_t/(2*delta_x)*A^2.
    This method is the arithmetic average of the Lax-Friedrichs method and the Lax-Wendroff method. 

    ...

    Attributes
    ----------
    quadrature_nodes: np.ndarray
        the nodes of the Gauss-Legendre quadrature rule on the interval [0,1]
    quadrature_weights: np.ndarray
        the weights of the Gauss-Legendre quadrature rule on the interval [0,1]

    
    Methods inherited from abtract parent class PVM
    ------------------------------------------------
    def __init__(self,nr_of_quadrature_points):
        initializes the PVM object
    def compute_fluctuation(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        computes the fluctuation between two cells with values value_left and value_right
    def compute_generalized_roe_and_viscosity(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        compute the generalized roe matrix and the viscosity matrix between two cells with values value_left and value_right
    def _compute_quadrature_points(self,nr_of_quadrature_points):
        computes the nodes and the weights for a Gauss-Legendre rule on the interval [0,1]
        with nr_of_quadrature_points quadrature points
    
    Methods implemented from abstract parent class PVM
    def compute_viscosity(self,roe_matrix,delta_t,delta_x):
        computes the viscosity matrix for the PRICE scheme

    """

    def compute_viscosity(self,
                          roe_matrix: np.ndarray,
                          delta_t: float,
                          delta_x: float):
        
        viscosity = 0.5*delta_x/delta_t*np.identity(np.shape(roe_matrix)[0])+0.5*delta_t/delta_x*roe_matrix@roe_matrix 
        return viscosity
    
class LF(PVM):
    """
    This class represents the Lax-Friedrichs scheme, a PVM scheme with viscosity function Q(A) = delta_x/delta_t*I. 

    ...

    Attributes
    ----------
    quadrature_nodes: np.ndarray
        the nodes of the Gauss-Legendre quadrature rule on the interval [0,1]
    quadrature_weights: np.ndarray
        the weights of the Gauss-Legendre quadrature rule on the interval [0,1]

    
    Methods inherited from abtract parent class PVM
    ------------------------------------------------
    def __init__(self,nr_of_quadrature_points):
        initializes the PVM object
    def compute_fluctuation(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        computes the fluctuation between two cells with values value_left and value_right
    def compute_generalized_roe_and_viscosity(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        compute the generalized roe matrix and the viscosity matrix between two cells with values value_left and value_right
    def _compute_quadrature_points(self,nr_of_quadrature_points):
        computes the nodes and the weights for a Gauss-Legendre rule on the interval [0,1]
        with nr_of_quadrature_points quadrature points

    Methods implemented from abstract parent class PVM
    def compute_viscosity(self,roe_matrix,delta_t,delta_x):
        computes the viscosity matrix for the Lax-Friedrichs scheme

    """

    def compute_viscosity(self,
                          roe_matrix: np.ndarray,
                          delta_t: float,
                          delta_x: float):
        
        viscosity = delta_x/delta_t*np.identity(np.shape(roe_matrix)[0])
        return viscosity

class Roe(PVM):
    """
    This class represents the Roe scheme, a PVM scheme with viscosity function Q(A) = |A|. 

    ...

    Attributes
    ----------
    quadrature_nodes: np.ndarray
        the nodes of the Gauss-Legendre quadrature rule on the interval [0,1]
    quadrature_weights: np.ndarray
        the weights of the Gauss-Legendre quadrature rule on the interval [0,1]

    
    Methods inherited from abtract parent class PVM
    ------------------------------------------------
    def __init__(self,nr_of_quadrature_points):
        initializes the PVM object
    def compute_fluctuation(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        computes the fluctuation between two cells with values value_left and value_right
    def compute_generalized_roe_and_viscosity(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        compute the generalized roe matrix and the viscosity matrix between two cells with values value_left and value_right
    def _compute_quadrature_points(self,nr_of_quadrature_points):
        computes the nodes and the weights for a Gauss-Legendre rule on the interval [0,1]
        with nr_of_quadrature_points quadrature points

    Methods implemented from abstract parent class PVM
    def compute_viscosity(self,roe_matrix,delta_t,delta_x):
        computes the viscosity matrix for the Roe scheme

    """

    def compute_viscosity(self,
                          roe_matrix: np.ndarray,
                          delta_t: float,
                          delta_x: float):

        # Eigen-decomposition: A = R D R^-1
        eigenvalues, R = np.linalg.eig(roe_matrix)
        
        # Construct |D|
        D_abs = np.diag(np.abs(eigenvalues))
        
        # Compute inverse of R
        R_inv = np.linalg.inv(R)
        
        # Return B = R |D| R^-1
        viscosity = R @ D_abs @ R_inv

        return viscosity
    
class Osher(PVM):
    """
    This class represents the method of Osher and Solomon, a PVM scheme with viscosity function 
    Q(A) = sum_{i=1}^{nr_of_quadrature_points}weight_i*|A(node_i)|. 

    ...

    Attributes
    ----------
    quadrature_nodes: np.ndarray
        the nodes of the Gauss-Legendre quadrature rule on the interval [0,1]
    quadrature_weights: np.ndarray
        the weights of the Gauss-Legendre quadrature rule on the interval [0,1]
    eigenstructure_available: boolean
        whether the eigenvalues and eigenvectors are given analytically or not
    compute_eigenvalues_and_eigenvectors: function
        the function that evaluates the analytical eigenvalues and the eigenvectors
    

    Methods inherited from abtract parent class PVM
    ------------------------------------------------
    def __init__(self,nr_of_quadrature_points):
        initializes the PVM object
    def compute_fluctuation(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        computes the fluctuation between two cells with values value_left and value_right
    def compute_generalized_roe_and_viscosity(self,value_left,value_right,system_matrix,direction,delta_t,delta_x):
        compute the generalized roe matrix and the viscosity matrix between two cells with values value_left and value_right
    def _compute_quadrature_points(self,nr_of_quadrature_points):
        computes the nodes and the weights for a Gauss-Legendre rule on the interval [0,1]
        with nr_of_quadrature_points quadrature points

    Methods implemented from abstract parent class PVM
    def compute_viscosity(self,roe_matrix,delta_t,delta_x):
        computes the viscosity matrix for the Roe scheme

    """

    def __init__(self,
                 nr_of_quadrature_points: int,
                 eigenstructure_available: bool,
                 compute_eigenvalues_and_eigenvectors: Callable[...,tuple[np.ndarray,np.ndarray]]):
        super().__init__(nr_of_quadrature_points)
        self.eigenstructure_available = eigenstructure_available
        self.compute_eigenvalues_and_eigenvectors = compute_eigenvalues_and_eigenvectors

    def compute_fluctuation(self,
                            value_left: np.ndarray,
                            value_right: np.ndarray,
                            system_matrix: Callable[...,np.ndarray],
                            delta_t: float,
                            delta_x: float) -> tuple[np.ndarray,np.ndarray]:

        generalized_roe = 0
        viscosity = 0
        for i in range(len(self.quadrature_nodes)):
            quadrature_point = (1-self.quadrature_nodes[i])*value_left+(self.quadrature_nodes[i])*value_right
            generalized_roe_point = system_matrix(quadrature_point)
            generalized_roe += self.quadrature_weights[i]*generalized_roe_point
            viscosity_point = self.compute_viscosity(generalized_roe_point,quadrature_point,delta_t,delta_x)
            viscosity += self.quadrature_weights[i]*viscosity_point
        viscosity = np.dot(viscosity,value_right-value_left)
        generalized_roe = np.dot(generalized_roe,value_right-value_left)
        fluctuation_min = (generalized_roe - viscosity)/2
        fluctuation_plus = (generalized_roe + viscosity)/2

        return fluctuation_min, fluctuation_plus
    
    def compute_generalized_roe_and_viscosity(self,
                                              value_left: np.ndarray,
                                              value_right: np.ndarray,
                                              system_matrix: np.ndarray,
                                              delta_t: float,
                                              delta_x: float):

        generalized_roe = 0
        viscosity = 0
        for i in range(len(self.quadrature_nodes)):
            quadrature_point = (1-self.quadrature_nodes[i])*value_left+(self.quadrature_nodes[i])*value_right
            generalized_roe_point = system_matrix(quadrature_point)
            generalized_roe += self.quadrature_weights[i]*generalized_roe_point
            viscosity_point = self.compute_viscosity(generalized_roe_point,quadrature_point,delta_t,delta_x)
            viscosity += self.quadrature_weights[i]*viscosity_point
        fluct_matrix_min = (generalized_roe - viscosity)/2
        fluct_matrix_plus = (generalized_roe + viscosity)/2

        return fluct_matrix_min, fluct_matrix_plus

    def compute_viscosity(self,
                          roe_matrix: np.ndarray,
                          values: np.ndarray,
                          delta_t: float,
                          delta_x: float):

        if self.eigenstructure_available:
            eigenvalue_diag_entries, right_eigenvectors_matrix = self.compute_eigenvalues_and_eigenvectors(values)

            D_abs = np.array(np.abs(eigenvalue_diag_entries))
            left_eigenvectors_matrix = np.linalg.inv(right_eigenvectors_matrix)
            viscosity = (right_eigenvectors_matrix*D_abs[None,:])@left_eigenvectors_matrix

        else:
            # Eigen-decomposition: A = R D R^-1
            eigenvalues, R = np.linalg.eig(roe_matrix)
            
            # Construct |D|
            D_abs = np.array(np.abs(eigenvalues))
            
            # Compute inverse of R
            R_inv = np.linalg.inv(R)
            
            # Return B = R |D| R^-1
            viscosity = (R*D_abs[None,:]) @ R_inv

        return viscosity