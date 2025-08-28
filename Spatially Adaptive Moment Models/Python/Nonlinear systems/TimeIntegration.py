from abc import ABC, abstractmethod
import numpy as np
import scipy.optimize as spopt
from collections.abc import Callable

class TimeIntegration(ABC):

    """
    This interface represents a time integrator.

    ...

    Attributes
    ----------
    None

    
    Abstract methods
    -------
    def integrate(self):
        integrates the equation in time
    """

    @abstractmethod
    def integrate(self):

        """
        Documented in the child classes

        """
        pass

class Implicit(TimeIntegration):

    """
    This abstract class represents an implicit time integrator.

    ...

    Attributes
    ----------
    linear : boolean
        True if the source term can be written in linear form, false if the source term can not be written in linear form

    Implemented methods from interface TimeIntegration
    -------------
    def integrate(self,initial_value,rhs_f,delta_t):
        integrates the equation dw/dt = rhs_f(w) with a time step of delta_t

            
    Instance methods    
    -----------------
    def _compute_residual(self,initial_value,rhs_f,delta_t):
        construct the residual function residual(x), the roots of which will be computed numerically
    """

    @abstractmethod
    def __init__(self,linear):
        self.linear = linear

    def integrate(self,
                  initial_value: np.array,
                  rhs_f: Callable[...,np.array],
                  delta_t: float) -> np.array:
        
        """
        integrates the equation dw/dt = rhs_f(w) starting from initial_value with a time step of delta_t

        Parameters
        ----------
        initial_value: np.array
            initial value w0
        rhs_f: function
            the right-hand side function describing the time evolution of w
        delta_t: float
            time step size
        
        Returns
        -------
        end_values: np.array
            final values 

        """
        if self.linear:
            end_values = rhs_f(initial_value,delta_t)@initial_value
        else: 
            end_values = spopt.newton(self._compute_residual(initial_value,rhs_f,delta_t),initial_value,maxiter=100)

        return end_values
    
    @abstractmethod
    def _compute_residual(self,
                  initial_value: np.array,
                  rhs_f: Callable[...,np.array],
                  delta_t: float) -> np.array:
        
        """
        construct the residual function residual(x), the roots of which will be computed numerically

        Parameters
        ----------
        initial_value: np.array
            initial value w0
        rhs_f: function
            the right-hand side function describing the time evolution of w
        delta_t: float
            time step size
        
        Returns
        -------
        residual: function
            function of which the roots will be computed (residual(x) == 0, where x is the numerical solution of the implicit solver) 

        """
        
        pass

class Explicit(TimeIntegration):

    """
    This abstract class represents an explicit time integrator.

    ...

    Attributes
    ----------
    None

    Implemented methods from interface TimeIntegration
    -------------
    def integrate(self,initial_value,rhs_f,delta_t):
        integrates the equation dw/dt = rhs_f(w) with a time step of delta_t

    """

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod    
    def integrate(self,
                  initial_value: np.array,
                  rhs_f: Callable[...,np.array],
                  delta_t: float) -> np.array:
        
        """
        integrates the equation dw/dt = rhs_f(w) starting from initial_value with a time step of delta_t

        Parameters
        ----------
        initial_value: np.array
            initial value w0
        rhs_f: function
            the right-hand side function describing the time evolution of w
        delta_t: float
            time step size
        
        Returns
        -------
        end_values: np.array
            final values 

        """
        
        pass

class ImplicitEuler(Implicit):

    def __init__(self,linear_source):
        super().__init__(linear_source)

    def _compute_residual(self,
                  initial_value: np.array,
                  rhs_f: Callable[...,np.array],
                  delta_t: float,
                  **kwargs) -> np.array:
        
        residual = lambda end_value : end_value - delta_t*rhs_f(end_value,delta_t,**kwargs) - initial_value

        return residual
    
class ExplicitEuler(Explicit):

    def __init__(self):
        pass

    def integrate(self,
                  initial_value: np.array,
                  rhs_f: Callable[...,np.array],
                  delta_t: float,
                  **kwargs) -> np.array:
        
        """
        integrates the equation dw/dt = rhs_f(w) starting from initial_value with a time step of delta_t

        Parameters
        ----------
        initial_value: np.array
            initial value w0
        rhs_f: function
            the right-hand side function describing the time evolution of w
        delta_t: float
            time step size
        
        Returns
        -------
        end_values: np.array
            final values 

        """
        
        end_values = initial_value + delta_t*rhs_f(initial_value,delta_t,**kwargs)

        return end_values