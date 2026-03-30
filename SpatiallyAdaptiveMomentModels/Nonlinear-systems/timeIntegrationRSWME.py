from abc import ABC, abstractmethod
import numpy as np
import scipy.optimize as spopt
from collections.abc import Callable

class IterationCounter:
    def __init__(self):
        self.count = 0

    def __call__(self, func):
        def wrapped_func(*args, **kwargs):
            self.count += 1
            return func(*args, **kwargs)
        return wrapped_func

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
    None

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
    def __init__(self):
        pass

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
        # Initialize the iteration counter
        counter = IterationCounter()

# Wrap the function to count iterations
        wrapped_function = counter(self._compute_residual(initial_value,rhs_f,delta_t))

        end_values = spopt.newton(wrapped_function,initial_value,maxiter=350)
        #print("Number of iterations:", counter.count)
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

    def __init__(self):
        pass

    def _compute_residual(self,
                  initial_value: np.array,
                  rhs_f: Callable[...,np.array],
                  delta_t: float) -> np.array:
        
        residual = lambda end_value : end_value - delta_t*rhs_f(end_value) - initial_value

        return residual
    
class ExplicitEuler(Explicit):

    def __init__(self):
        pass

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
        
        end_values = initial_value + delta_t*rhs_f(initial_value)

        return end_values

class ImplicitDirect(Explicit):
    def __init__(self,viscosity: float,
                slip_length: float, N:int):
        self.viscosity = viscosity
        self.slip_length = slip_length
        self.N=N
        pass

    def integrate(self,
                  initial_value: np.array,
                  rhs_f: Callable[...,np.array],
                  delta_t: float) -> np.array:
        
        """
        integrates the equation dw/dt = rhs_f(w), for f linear.

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
        end_values = np.zeros(initial_value.shape)      
        end_values[0] = initial_value[0]
        II = np.identity(self.N+1)
        D1 = (delta_t*self.viscosity)/(end_values[0]*self.slip_length)
        D2 = (delta_t*self.viscosity)/(end_values[0]*end_values[0])
        S1 = np.tile((D1 * (1 + 2 * np.arange(self.N + 1))), (self.N + 1, 1)).T
        i_vals = np.arange(0, self.N + 1)[:, None]  
        j_vals = np.arange(0, self.N + 1)
        min_ij = np.minimum(i_vals, j_vals)
        S2 = np.where(((i_vals - j_vals) % 2 == 0) ,
                       D2*(2*i_vals+1)*2 * min_ij * (min_ij + 1),
                      0.0)
        print(II+S1+S2)
        print(initial_value[1:])
        print(end_values[1:])
        end_values[1:]= np.linalg.solve((II+S1+S2),initial_value[1:])
        return end_values

class ImplicitFast(Explicit):
    def __init__(self,viscosity: float,
                slip_length: float, N:int):
        self.viscosity = viscosity
        self.slip_length = slip_length
        self.N=N
        pass

    def integrate(self,
                  initial_value: np.array,
                  rhs_f: Callable[...,np.array],
                  delta_t: float) -> np.array:
        
        """
        integrates the equation dw/dt = rhs_f(w), for f linear.

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
        end_values = np.zeros(initial_value.shape)      
        h = initial_value[0]
        n = self.viscosity
        l = self.slip_length
        dt = delta_t
        if self.N == 1:
            det = 12 * dt**2 * n**2 + 4 * dt * h**2 * n + 12 * dt * h * l * n + h**3 * l
            S_inv = (1/det)*np.array([[h * (3 * dt * h * n + 12 * dt * l * n + h**2 * l),-dt * h**2 * n],
            [-3 * dt * h**2 * n , h**2 * (dt * n + h * l)]])
        if self.N==2:
            det= 720*dt**3*n**3 + 312*dt**2*h**2*n**2 + 720*dt**2*h*l*n**2 + 9*dt*h**4*n + 72*dt*h**3*l*n + h**5*l
            S_inv = (1/det)*np.array( [[h*(240*dt**2*h*n**2 + 720*dt**2*l*n**2 + 8*dt*h**3*n + 72*dt*h**2*l*n + h**4*l), dt*h**2*n*(-60*dt*n - h**2), dt*h**2*n*(-12*dt*n - h**2)], [3*dt*h**2*n*(-60*dt*n - h**2), h**2*(60*dt**2*n**2 + 6*dt*h**2*n + 60*dt*h*l*n + h**3*l), -3*dt*h**4*n], [5*dt*h**2*n*(-12*dt*n - h**2), -5*dt*h**4*n, h**2*(12*dt**2*n**2 + 4*dt*h**2*n + 12*dt*h*l*n + h**3*l)]])
        if self.N==3:
            det = (
                100800 * dt**4 * n**4 + 46080 * dt**3 * h**2 * n**3 + 100800 * dt**3 * h * l * n**3 +2160 * dt**2 * h**4 * n**2 +
                12480 * dt**2 * h**3 * l * n**2 + 16 * dt * h**6 * n + 240 * dt * h**5 * l * n + h**7 * l)
            S_inv =S_inv = (1 / det) * np.array([
                [
                    h * (33600 * dt**3 * h * n**3 + 100800 * dt**3 * l * n**3 + 1920 * dt**2 * h**3 * n**2 + 12480 * dt**2 * h**2 * l * n**2 + 15 * dt * h**5 * n + 240 * dt * h**4 * l * n + h**6 * l),
                    dt * h**2 * n * (-8400 * dt**2 * n**2 - 200 * dt * h**2 * n - h**4),
                    dt * h**2 * n * (-1680 * dt**2 * n**2 - 180 * dt * h**2 * n - h**4),
                    dt * h**4 * n * (-60 * dt * n - h**2)
                ],
                [
                    3 * dt * h**2 * n * (-8400 * dt**2 * n**2 - 200 * dt * h**2 * n - h**4),
                    h**2 * (10080 * dt**3 * n**3 + 1488 * dt**2 * h**2 * n**2 + 10080 * dt**2 * h * l * n**2 + 13 * dt * h**4 * n + 228 * dt * h**3 * l * n + h**5 * l),
                    3 * dt * h**4 * n * (-140 * dt * n - h**2),
                    3 * dt * h**2 * n * (-240 * dt**2 * n**2 - 84 * dt * h**2 * n - 240 * dt * h * l * n - h**4 - 4 * h**3 * l)
                ],
                [
                    5 * dt * h**2 * n * (-1680 * dt**2 * n**2 - 180 * dt * h**2 * n - h**4),
                    5 * dt * h**4 * n * (-140 * dt * n - h**2),
                    h**2 * (1680 * dt**3 * n**3 + 600 * dt**2 * h**2 * n**2 + 1680 * dt**2 * h * l * n**2 + 11 * dt * h**4 * n + 180 * dt * h**3 * l * n + h**5 * l),
                    - 5 * dt * h**6 * n
                ],
                [
                    7 * dt * h**4 * n * (-60 * dt * n - h**2),
                    7 * dt * h**2 * n * (-240 * dt**2 * n**2 - 84 * dt * h**2 * n - 240 * dt * h * l * n - h**4 - 4 * h**3 * l),
                    -7 * dt * h**6 * n,
                    h**2 * (720 * dt**3 * n**3 + 312 * dt**2 * h**2 * n**2 + 720 * dt**2 * h * l * n**2 + 9 * dt * h**4 * n + 72 * dt * h**3 * l * n + h**5 * l)
                ]])
        if self.N==4:
            det = (
                25401600 * dt**5 * n**5 + 11894400 * dt**4 * h**2 * n**4 + 25401600 * dt**4 * h * l * n**4 + 670320 * dt**3 * h**4 * n**3 +
                3427200 * dt**3 * h**3 * l * n**3 + 9000 * dt**2 * h**6 * n**2 + 92400 * dt**2 * h**5 * l * n**2 + 25 * dt * h**8 * n +
                600 * dt * h**7 * l * n + h**9 * l)
            S_inv = (1 / det) * np.array([
                [
                    h * (8467200 * dt**4 * h * n**4 + 25401600 * dt**4 * l * n**4 + 577920 * dt**3 * h**3 * n**3 + 3427200 * dt**3 * h**2 * l * n**3 + 8400 * dt**2 * h**5 * n**2 + 92400 * dt**2 * h**4 * l * n**2 + 24 * dt * h**7 * n + 600 * dt * h**6 * l * n + h**8 * l),
                    dt * h**2 * n * (-2116800 * dt**3 * n**3 - 73920 * dt**2 * h**2 * n**2 - 560 * dt * h**4 * n - h**6),
                    dt * h**2 * n * (-423360 * dt**3 * n**3 - 47040 * dt**2 * h**2 * n**2 - 432 * dt * h**4 * n - h**6),
                    dt * h**4 * n * (-15120 * dt**2 * n**2 - 420 * dt * h**2 * n - h**4),
                    dt * h**4 * n * (-1680 * dt**2 * n**2 - 180 * dt * h**2 * n - h**4)
                ],
                [
                    3 * dt * h**2 * n * (-2116800 * dt**3 * n**3 - 73920 * dt**2 * h**2 * n**2 - 560 * dt * h**4 * n - h**6),
                    h**2 * (2540160 * dt**4 * n**4 + 403200 * dt**3 * h**2 * n**3 + 2540160 * dt**3 * h * l * n**3 + 7140 * dt**2 * h**4 * n**2 + 85680 * dt**2 * h**3 * l * n**2 + 22 * dt * h**6 * n + 588 * dt * h**5 * l * n + h**7 * l),
                    3 * dt * h**4 * n * (-35280 * dt**2 * n**2 - 392 * dt * h**2 * n - h**4),
                    3 * dt * h**2 * n * (-60480 * dt**3 * n**3 - 21840 * dt**2 * h**2 * n**2 - 60480 * dt**2 * h * l * n**2 - 480 * dt * h**4 * n - 1680 * dt * h**3 * l * n - h**6 - 4 * h**5 * l),
                    3 * dt * h**6 * n * (-140 * dt * n - h**2)
                ],
                [
                    5 * dt * h**2 * n * (-423360 * dt**3 * n**3 - 47040 * dt**2 * h**2 * n**2 - 432 * dt * h**4 * n - h**6),
                    5 * dt * h**4 * n * (-35280 * dt**2 * n**2 - 392 * dt * h**2 * n - h**4),
                    h**2 * (604800 * dt**4 * n**4 + 232800 * dt**3 * h**2 * n**3 + 604800 * dt**3 * h * l * n**3 + 6180 * dt**2 * h**4 * n**2 + 66480 * dt**2 * h**3 * l * n**2 + 20 * dt * h**6 * n + 540 * dt * h**5 * l * n + h**7 * l),
                    -5 * dt * h**6 * n * (252 * dt * n + h**2),
                    5 * dt * h**2 * n * (-20160 * dt**3 * n**3 - 8880 * dt**2 * h**2 * n**2 - 20160 * dt**2 * h * l * n**2 - 312 * dt * h**4 * n - 2160 * dt * h**3 * l * n - h**6 - 12 * h**5 * l)
                ],
                [
                    7 * dt * h**4 * n * (-15120 * dt**2 * n**2 - 420 * dt * h**2 * n - h**4),
                    7 * dt * h**2 * n * (-60480 * dt**3 * n**3 - 21840 * dt**2 * h**2 * n**2 - 60480 * dt**2 * h * l * n**2 - 480 * dt * h**4 * n - 1680 * dt * h**3 * l * n - h**6 - 4 * h**5 * l),
                    -7 * dt * h**6 * n * (252 * dt * n + h**2),
                    h**2 * (181440 * dt**4 * n**4 + 80640 * dt**3 * h**2 * n**3 + 181440 * dt**3 * h * l * n**3 + 3120 * dt**2 * h**4 * n**2 + 20160 * dt**2 * h**3 * l * n**2 + 18 * dt * h**6 * n + 432 * dt * h**5 * l * n + h**7 * l),
                    -7 * dt * h**8 * n
                ],
                [
                    9 * dt * h**4 * n * (-1680 * dt**2 * n**2 - 180 * dt * h**2 * n - h**4),
                    9 * dt * h**6 * n * (-140 * dt * n - h**2),
                    9 * dt * h**2 * n * (-20160 * dt**3 * n**3 - 8880 * dt**2 * h**2 * n**2 - 20160 * dt**2 * h * l * n**2 - 312 * dt * h**4 * n - 2160 * dt * h**3 * l * n - h**6 - 12 * h**5 * l),
                    -9 * dt * h**8 * n,
                    h**2 * (100800 * dt**4 * n**4 + 46080 * dt**3 * h**2 * n**3 + 100800 * dt**3 * h * l * n**3 + 2160 * dt**2 * h**4 * n**2 + 12480 * dt**2 * h**3 * l * n**2 + 16 * dt * h**6 * n + 240 * dt * h**5 * l * n + h**7 * l)
                ]])
        if self.N==5:
                det = (
                    10059033600 * dt**6 * n**6 + 4775500800 * dt**5 * h**2 * n**5 + 10059033600 * dt**5 * h * l * n**5 + 295747200 * dt**4 * h**4 * n**4 +
                    1422489600 * dt**4 * h**3 * l * n**4 + 5160960 * dt**3 * h**6 * n**3 + 45118080 * dt**3 * h**5 * l * n**3 + 28140 * dt**2 * h**8 * n**2 +
                    443520 * dt**2 * h**7 * l * n**2 + 36 * dt * h**10 * n + 1260 * dt * h**9 * l * n + h**11 * l)
                S_inv = (1 / denominator) * np.array([
                    [
                        h * (3353011200 * dt**5 * h * n**5 + 10059033600 * dt**5 * l * n**5 + 250629120 * dt**4 * h**3 * n**4 + 1422489600 * dt**4 * h**2 * l * n**4 + 4717440 * dt**3 * h**5 * n**3 + 45118080 * dt**3 * h**4 * l * n**3 + 26880 * dt**2 * h**7 * n**2 + 443520 * dt**2 * h**6 * l * n**2 + 35 * dt * h**9 * n + 1260 * dt * h**8 * l * n + h**10 * l),
                        dt * h**2 * n * (-838252800 * dt**4 * n**4 - 34715520 * dt**3 * h**2 * n**3 - 388080 * dt**2 * h**4 * n**2 - 1176 * dt * h**6 * n - h**8),
                        dt * h**2 * n * (-167650560 * dt**4 * n**4 - 19716480 * dt**3 * h**2 * n**3 - 287280 * dt**2 * h**4 * n**2 - 1092 * dt * h**6 * n - h**8),
                        dt * h**4 * n * (-5987520 * dt**3 * n**3 - 181440 * dt**2 * h**2 * n**2 - 816 * dt * h**4 * n - h**6),
                        dt * h**4 * n * (-665280 * dt**3 * n**3 - 75600 * dt**2 * h**2 * n**2 - 840 * dt * h**4 * n - h**6),
                        dt * h**6 * n * (-15120 * dt**2 * n**2 - 420 * dt * h**2 * n - h**4)
                    ],
                    [
                        3 * dt * h**2 * n * (-838252800 * dt**4 * n**4 - 34715520 * dt**3 * h**2 * n**3 - 388080 * dt**2 * h**4 * n**2 - 1176 * dt * h**6 * n - h**8),
                        h**2 * (1005903360 * dt**5 * n**5 + 166199040 * dt**4 * h**2 * n**4 + 1005903360 * dt**4 * h * l * n**4 + 3840480 * dt**3 * h**4 * n**3 + 40461120 * dt**3 * h**3 * l * n**3 + 24432 * dt**2 * h**6 * n**2 + 429408 * dt**2 * h**5 * l * n**2 + 33 * dt * h**8 * n + 1248 * dt * h**7 * l * n + h**9 * l),
                        3 * dt * h**4 * n * (-13970880 * dt**3 * n**3 - 245952 * dt**2 * h**2 * n**2 - 1008 * dt * h**4 * n - h**6),
                        3 * dt * h**2 * n * (-23950080 * dt**4 * n**4 - 8709120 * dt**3 * h**2 * n**3 - 23950080 * dt**3 * h * l * n**3 - 211920 * dt**2 * h**4 * n**2 - 725760 * dt**2 * h**3 * l * n**2 - 876 * dt * h**6 * n - 3264 * dt * h**5 * l * n - h**8 - 4 * h**7 * l),
                        3 * dt * h**6 * n * (-55440 * dt**2 * n**2 - 756 * dt * h**2 * n - h**4),
                        3 * dt * h**4 * n * (-60480 * dt**3 * n**3 - 21840 * dt**2 * h**2 * n**2 - 60480 * dt**2 * h * l * n**2 - 480 * dt * h**4 * n - 1680 * dt * h**3 * l * n - h**6 - 4 * h**5 * l)
                    ],
                    [
                        5 * dt * h**2 * n * (-167650560 * dt**4 * n**4 - 19716480 * dt**3 * h**2 * n**3 - 287280 * dt**2 * h**4 * n**2 - 1092 * dt * h**6 * n - h**8),
                        5 * dt * h**4 * n * (-13970880 * dt**3 * n**3 - 245952 * dt**2 * h**2 * n**2 - 1008 * dt * h**4 * n - h**6),
                        h**2 * (239500800 * dt**5 * n**5 + 93744000 * dt**4 * h**2 * n**4 + 239500800 * dt**4 * h * l * n**4 + 3039120 * dt**3 * h**4 * n**3 + 27881280 * dt**3 * h**3 * l * n**3 + 21360 * dt**2 * h**6 * n**2 + 378000 * dt**2 * h**5 * l * n**2 + 31 * dt * h**8 * n + 1200 * dt * h**7 * l * n + h**9 * l),
                        -5 * dt * h**6 * n * (99792 * dt**2 * n**2 + 648 * dt * h**2 * n + h**4),
                        5 * dt * h**2 * n * (-7983360 * dt**4 * n**4 - 3568320 * dt**3 * h**2 * n**3 - 7983360 * dt**3 * h * l * n**3 - 146160 * dt**2 * h**4 * n**2 - 907200 * dt**2 * h**3 * l * n**2 - 1104 * dt * h**6 * n - 10080 * dt * h**5 * l * n - h**8 - 12 * h**7 * l),
                        5 * dt * h**8 * n * (-252 * dt * n - h**2)
                    ],
                    [
                        7 * dt * h**4 * n * (-5987520 * dt**3 * n**3 - 181440 * dt**2 * h**2 * n**2 - 816 * dt * h**4 * n - h**6),
                        7 * dt * h**2 * n * (-23950080 * dt**4 * n**4 - 8709120 * dt**3 * h**2 * n**3 - 23950080 * dt**3 * h * l * n**3 - 211920 * dt**2 * h**4 * n**2 - 725760 * dt**2 * h**3 * l * n**2 - 876 * dt * h**6 * n - 3264 * dt * h**5 * l * n - h**8 - 4 * h**7 * l),
                        -7 * dt * h**6 * n * (99792 * dt**2 * n**2 + 648 * dt * h**2 * n + h**4),
                        h**2 * (111767040 * dt**5 * n**5 + 50520960 * dt**4 * h**2 * n**4 + 111767040 * dt**4 * h * l * n**4 + 2242800 * dt**3 * h**4 * n**3 + 13265280 * dt**3 * h**3 * l * n**3 + 19488 * dt**2 * h**6 * n**2 + 304752 * dt**2 * h**5 * l * n**2 + 29 * dt * h**8 * n + 1092 * dt * h**7 * l * n + h**9 * l),
                        -7 * dt * h**8 * n * (396 * dt * n + h**2),
                        7 * dt * h**2 * n * (-3628800 * dt**4 * n**4 - 1673280 * dt**3 * h**2 * n**3 - 3628800 * dt**3 * h * l * n**3 - 84240 * dt**2 * h**4 * n**2 - 463680 * dt**2 * h**3 * l * n**2 - 840 * dt * h**6 * n - 10320 * dt * h**5 * l * n - h**8 - 24 * h**7 * l)
                    ],
                    [
                        9 * dt * h**4 * n * (-665280 * dt**3 * n**3 - 75600 * dt**2 * h**2 * n**2 - 840 * dt * h**4 * n - h**6),
                        9 * dt * h**6 * n * (-55440 * dt**2 * n**2 - 756 * dt * h**2 * n - h**4),
                        9 * dt * h**2 * n * (-7983360 * dt**4 * n**4 - 3568320 * dt**3 * h**2 * n**3 - 7983360 * dt**3 * h * l * n**3 - 146160 * dt**2 * h**4 * n**2 - 907200 * dt**2 * h**3 * l * n**2 - 1104 * dt * h**6 * n - 10080 * dt * h**5 * l * n - h**8 - 12 * h**7 * l),
                        -9 * dt * h**8 * n * (396 * dt * n + h**2),
                        h**2 * (39916800 * dt**5 * n**5 + 18506880 * dt**4 * h**2 * n**4 + 39916800 * dt**4 * h * l * n**4 + 972720 * dt**3 * h**4 * n**3 + 5201280 * dt**3 * h**3 * l * n**3 + 11400 * dt**2 * h**6 * n**2 + 126000 * dt**2 * h**5 * l * n**2 + 27 * dt * h**8 * n + 900 * dt * h**7 * l * n + h**9 * l),
                        -9 * dt * h**10 * n
                    ],
                    [
                        11 * dt * h**6 * n * (-15120 * dt**2 * n**2 - 420 * dt * h**2 * n - h**4),
                        11 * dt * h**4 * n * (-60480 * dt**3 * n**3 - 21840 * dt**2 * h**2 * n**2 - 60480 * dt**2 * h * l * n**2 - 480 * dt * h**4 * n - 1680 * dt * h**3 * l * n - h**6 - 4 * h**5 * l),
                        11 * dt * h**8 * n * (-252 * dt * n - h**2),
                        11 * dt * h**2 * n * (-3628800 * dt**4 * n**4 - 1673280 * dt**3 * h**2 * n**3 - 3628800 * dt**3 * h * l * n**3 - 84240 * dt**2 * h**4 * n**2 - 463680 * dt**2 * h**3 * l * n**2 - 840 * dt * h**6 * n - 10320 * dt * h**5 * l * n - h**8 - 24 * h**7 * l),
                        -11 * dt * h**10 * n,
                        h**2 * (25401600 * dt**5 * n**5 + 11894400 * dt**4 * h**2 * n**4 + 25401600 * dt**4 * h * l * n**4 + 670320 * dt**3 * h**4 * n**3 + 3427200 * dt**3 * h**3 * l * n**3 + 9000 * dt**2 * h**6 * n**2 + 92400 * dt**2 * h**5 * l * n**2 + 25 * dt * h**8 * n + 600 * dt * h**7 * l * n + h**9 * l)                                                                                                                                     
                    ]])
        if self.N==6:
            denominator=(
                5753767219200 * dt**7 * n**7 + 2756175206400 * dt**6 * h**2 * n**6 + 5753767219200 * dt**6 * h * l * n**6 + 180783187200 * dt**5 * h**4 * n**5 
                + 838252800000 * dt**5 * h**3 * l * n**5 + 3648274560 * dt**4 * h**6 * n**4 + 29227080960 * dt**4 * h**5 * l * n**4 + 27306720 * dt**3 * h**8 * n**3 + 356469120 * dt**3 * h**7 * l * n**3 + 72912 * dt**2 * h**10 * n**2 + 1622880 * dt**2 * h**9 * l * n**2 + 49 * dt * h**12 * n + 2352 * dt * h**11 * l * n + h**13 * l)
            S_inv=(1/denominator) * np.array([[h*(1917922406400*dt**6*h*n**6 + 5753767219200*dt**6*l*n**6 + 151556106240*dt**5*h**3*n**5 + 838252800000*dt**5*h**2*l*n**5 + 3291805440*dt**4*h**5*n**4 + 29227080960*dt**4*h**4*l*n**4 + 25683840*dt**3*h**7*n**3 + 356469120*dt**3*h**6*l*n**3 + 70560*dt**2*h**9*n**2 + 1622880*dt**2*h**8*l*n**2 + 48*dt*h**11*n + 2352*dt*h**10*l*n + h**12*l), dt*h**2*n*(-479480601600*dt**5*n**5 - 21906339840*dt**4*h**2*n**4 - 302037120*dt**3*h**4*n**3 - 1475712*dt**2*h**6*n**2 - 2268*dt*h**8*n - h**10), dt*h**2*n*(-95896120320*dt**5*n**5 - 11687639040*dt**4*h**2*n**4 - 211559040*dt**3*h**4*n**3 - 1217664*dt**2*h**6*n**2 - 2028*dt*h**8*n - h**10), dt*h**4*n*(-3424861440*dt**4*n**4 - 118419840*dt**3*h**2*n**3 - 875952*dt**2*h**4*n**2 - 1908*dt*h**6*n - h**8), dt*h**4*n*(-380540160*dt**4*n**4 - 43908480*dt**3*h**2*n**3 - 556080*dt**2*h**4*n**2 - 1412*dt*h**6*n - h**8), dt*h**6*n*(-8648640*dt**3*n**3 - 277200*dt**2*h**2*n**2 - 1512*dt*h**4*n - h**6), dt*h**6*n*(-665280*dt**3*n**3 - 75600*dt**2*h**2*n**2 - 840*dt*h**4*n - h**6)], [3*dt*h**2*n*(-479480601600*dt**5*n**5 - 21906339840*dt**4*h**2*n**4 - 302037120*dt**3*h**4*n**3 - 1475712*dt**2*h**6*n**2 - 2268*dt*h**8*n - h**10), h**2*(575376721920*dt**6*n**6 + 97524725760*dt**5*h**2*n**5 + 575376721920*dt**5*h*l*n**5 + 2597253120*dt**4*h**4*n**4 + 25602635520*dt**4*h**3*l*n**4 + 22474368*dt**3*h**6*n**3 + 338760576*dt**3*h**5*l*n**3 + 65772*dt**2*h**8*n**2 + 1595664*dt**2*h**7*l*n**2 + 46*dt*h**10*n + 2340*dt*h**9*l*n + h**11*l), 3*dt*h**4*n*(-7991343360*dt**4*n**4 - 174835584*dt**3*h**2*n**3 - 1097712*dt**2*h**4*n**2 - 1944*dt*h**6*n - h**8), 3*dt*h**2*n*(-13699445760*dt**5*n**5 - 5040161280*dt**4*h**2*n**4 - 13699445760*dt**4*h*l*n**4 - 142369920*dt**3*h**4*n**3 - 473679360*dt**3*h**3*l*n**3 - 970704*dt**2*h**6*n**2 - 3503808*dt**2*h**5*l*n**2 - 2020*dt*h**8*n - 7632*dt*h**7*l*n - h**10 - 4*h**9*l), 3*dt*h**6*n*(-31711680*dt**3*n**3 - 487872*dt**2*h**2*n**2 - 1328*dt*h**4*n - h**6), 3*dt*h**4*n*(-34594560*dt**4*n**4 - 12640320*dt**3*h**2*n**3 - 34594560*dt**3*h*l*n**3 - 327600*dt**2*h**4*n**2 - 1108800*dt**2*h**3*l*n**2 - 1624*dt*h**6*n - 6048*dt*h**5*l*n - h**8 - 4*h**7*l), 3*dt*h**8*n*(-55440*dt**2*n**2 - 756*dt*h**2*n - h**4)], [5*dt*h**2*n*(-95896120320*dt**5*n**5 - 11687639040*dt**4*h**2*n**4 - 211559040*dt**3*h**4*n**3 - 1217664*dt**2*h**6*n**2 - 2028*dt*h**8*n - h**10), 5*dt*h**4*n*(-7991343360*dt**4*n**4 - 174835584*dt**3*h**2*n**3 - 1097712*dt**2*h**4*n**2 - 1944*dt*h**6*n - h**8), h**2*(136994457600*dt**6*n**6 + 54207014400*dt**5*h**2*n**5 + 136994457600*dt**5*h*l*n**5 + 1966567680*dt**4*h**4*n**4 + 16533538560*dt**4*h**3*l*n**4 + 19297440*dt**3*h**6*n**3 + 283409280*dt**3*h**5*l*n**3 + 61452*dt**2*h**8*n**2 + 1501200*dt**2*h**7*l*n**2 + 44*dt*h**10*n + 2292*dt*h**9*l*n + h**11*l), 5*dt*h**6*n*(-57081024*dt**3*n**3 - 614592*dt**2*h**2*n**2 - 1584*dt*h**4*n - h**6), 5*dt*h**2*n*(-4566481920*dt**5*n**5 - 2049062400*dt**4*h**2*n**4 - 4566481920*dt**4*h*l*n**4 - 87171840*dt**3*h**4*n**3 - 526901760*dt**3*h**3*l*n**3 - 777648*dt**2*h**6*n**2 - 6672960*dt**2*h**5*l*n**2 - 1676*dt*h**8*n - 16944*dt*h**7*l*n - h**10 - 12*h**9*l), 5*dt*h**8*n*(-144144*dt**2*n**2 - 1188*dt*h**2*n - h**4), 5*dt*h**4*n*(-7983360*dt**4*n**4 - 3568320*dt**3*h**2*n**3 - 7983360*dt**3*h*l*n**3 - 146160*dt**2*h**4*n**2 - 907200*dt**2*h**3*l*n**2 - 1104*dt*h**6*n - 10080*dt*h**5*l*n - h**8 - 12*h**7*l)], [7*dt*h**4*n*(-3424861440*dt**4*n**4 - 118419840*dt**3*h**2*n**3 - 875952*dt**2*h**4*n**2 - 1908*dt*h**6*n - h**8), 7*dt*h**2*n*(-13699445760*dt**5*n**5 - 5040161280*dt**4*h**2*n**4 - 13699445760*dt**4*h*l*n**4 - 142369920*dt**3*h**4*n**3 - 473679360*dt**3*h**3*l*n**3 - 970704*dt**2*h**6*n**2 - 3503808*dt**2*h**5*l*n**2 - 2020*dt*h**8*n - 7632*dt*h**7*l*n - h**10 - 4*h**9*l), 7*dt*h**6*n*(-57081024*dt**3*n**3 - 614592*dt**2*h**2*n**2 - 1584*dt*h**4*n - h**6), h**2*(63930746880*dt**6*n**6 + 29171197440*dt**5*h**2*n**5 + 63930746880*dt**5*h*l*n**5 + 1405736640*dt**4*h**4*n**4 + 7860948480*dt**4*h**3*l*n**4 + 16347072*dt**3*h**6*n**3 + 206103744*dt**3*h**5*l*n**3 + 54432*dt**2*h**8*n**2 + 1300656*dt**2*h**7*l*n**2 + 42*dt*h**10*n + 2184*dt*h**9*l*n + h**11*l), 7*dt*h**8*n*(-226512*dt**2*n**2 - 968*dt*h**2*n - h**4), 7*dt*h**2*n*(-2075673600*dt**5*n**5 - 965986560*dt**4*h**2*n**4 - 2075673600*dt**4*h*l*n**4 - 52254720*dt**3*h**4*n**3 - 274095360*dt**3*h**3*l*n**3 - 677040*dt**2*h**6*n**2 - 7015680*dt**2*h**5*l*n**2 - 2244*dt*h**8*n - 36528*dt*h**7*l*n - h**10 - 24*h**9*l), 7*dt*h**10*n*(-396*dt*n - h**2)], [9*dt*h**4*n*(-380540160*dt**4*n**4 - 43908480*dt**3*h**2*n**3 - 556080*dt**2*h**4*n**2 - 1412*dt*h**6*n - h**8), 9*dt*h**6*n*(-31711680*dt**3*n**3 - 487872*dt**2*h**2*n**2 - 1328*dt*h**4*n - h**6), 9*dt*h**2*n*(-4566481920*dt**5*n**5 - 2049062400*dt**4*h**2*n**4 - 4566481920*dt**4*h*l*n**4 - 87171840*dt**3*h**4*n**3 - 526901760*dt**3*h**3*l*n**3 - 777648*dt**2*h**6*n**2 - 6672960*dt**2*h**5*l*n**2 - 1676*dt*h**8*n - 16944*dt*h**7*l*n - h**10 - 12*h**9*l), 9*dt*h**8*n*(-226512*dt**2*n**2 - 968*dt*h**2*n - h**4), h**2*(37362124800*dt**6*n**6 + 17466140160*dt**5*h**2*n**5 + 37362124800*dt**5*h*l*n**5 + 975360960*dt**4*h**4*n**4 + 5012098560*dt**4*h**3*l*n**4 + 13543200*dt**3*h**6*n**3 + 134930880*dt**3*h**5*l*n**3 + 51024*dt**2*h**8*n**2 + 1099440*dt**2*h**7*l*n**2 + 40*dt*h**10*n + 1992*dt*h**9*l*n + h**11*l), 9*dt*h**10*n*(-572*dt*n - h**2), 9*dt*h**2*n*(-1117670400*dt**5*n**5 - 526176000*dt**4*h**2*n**4 - 1117670400*dt**4*h*l*n**4 - 30804480*dt**3*h**4*n**3 - 153619200*dt**3*h**3*l*n**3 - 465360*dt**2*h**6*n**2 - 4435200*dt**2*h**5*l*n**2 - 1860*dt*h**8*n - 35280*dt*h**7*l*n - h**10 - 40*h**9*l)], [11*dt*h**6*n*(-8648640*dt**3*n**3 - 277200*dt**2*h**2*n**2 - 1512*dt*h**4*n - h**6), 11*dt*h**4*n*(-34594560*dt**4*n**4 - 12640320*dt**3*h**2*n**3 - 34594560*dt**3*h*l*n**3 - 327600*dt**2*h**4*n**2 - 1108800*dt**2*h**3*l*n**2 - 1624*dt*h**6*n - 6048*dt*h**5*l*n - h**8 - 4*h**7*l), 11*dt*h**8*n*(-144144*dt**2*n**2 - 1188*dt*h**2*n - h**4), 11*dt*h**2*n*(-2075673600*dt**5*n**5 - 965986560*dt**4*h**2*n**4 - 2075673600*dt**4*h*l*n**4 - 52254720*dt**3*h**4*n**3 - 274095360*dt**3*h**3*l*n**3 - 677040*dt**2*h**6*n**2 - 7015680*dt**2*h**5*l*n**2 - 2244*dt*h**8*n - 36528*dt*h**7*l*n - h**10 - 24*h**9*l), 11*dt*h**10*n*(-572*dt*n - h**2), h**2*(14529715200*dt**6*n**6 + 6865689600*dt**5*h**2*n**5 + 14529715200*dt**5*h*l*n**5 + 412352640*dt**4*h**4*n**4 + 2022451200*dt**4*h**3*l*n**4 + 6720000*dt**3*h**6*n**3 + 61084800*dt**3*h**5*l*n**3 + 33180*dt**2*h**8*n**2 + 551040*dt**2*h**7*l*n**2 + 38*dt*h**10*n + 1692*dt*h**9*l*n + h**11*l), -11*dt*h**12*n], [13*dt*h**6*n*(-665280*dt**3*n**3 - 75600*dt**2*h**2*n**2 - 840*dt*h**4*n - h**6), 13*dt*h**8*n*(-55440*dt**2*n**2 - 756*dt*h**2*n - h**4), 13*dt*h**4*n*(-7983360*dt**4*n**4 - 3568320*dt**3*h**2*n**3 - 7983360*dt**3*h*l*n**3 - 146160*dt**2*h**4*n**2 - 907200*dt**2*h**3*l*n**2 - 1104*dt*h**6*n - 10080*dt*h**5*l*n - h**8 - 12*h**7*l), 13*dt*h**10*n*(-396*dt*n - h**2), 13*dt*h**2*n*(-1117670400*dt**5*n**5 - 526176000*dt**4*h**2*n**4 - 1117670400*dt**4*h*l*n**4 - 30804480*dt**3*h**4*n**3 - 153619200*dt**3*h**3*l*n**3 - 465360*dt**2*h**6*n**2 - 4435200*dt**2*h**5*l*n**2 - 1860*dt*h**8*n - 35280*dt*h**7*l*n - h**10 - 40*h**9*l), -13*dt*h**12*n, h**2*(10059033600*dt**6*n**6 + 4775500800*dt**5*h**2*n**5 + 10059033600*dt**5*h*l*n**5 + 295747200*dt**4*h**4*n**4 + 1422489600*dt**4*h**3*l*n**4 + 5160960*dt**3*h**6*n**3 + 45118080*dt**3*h**5*l*n**3 + 28140*dt**2*h**8*n**2 + 443520*dt**2*h**7*l*n**2 + 36*dt*h**10*n + 1260*dt*h**9*l*n + h**11*l)]])
            
            

        end_values[0] = initial_value[0]
        end_values[1:] = S_inv.dot(initial_value[1:])
        return end_values