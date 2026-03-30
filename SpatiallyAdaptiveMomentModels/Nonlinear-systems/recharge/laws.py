#  Packages
import numpy as np

class HortonInfiltration(object):
    """
    Simple Horton infiltration model for the first recharge test case.

    Parameters
    ----------
    f0 : float
        Initial infiltration capacity.
    fc : float
        Final infiltration capacity.
    k : float
        Exponential decay coefficient.
    eps : float, optional
        Small positive number to avoid division by zero.
    """

    # Constructor
    def __init__(self, f0 : float, fc : float, k : float, eps : float = 1e-14):
        self.f0 = f0
        self.fc = fc
        self.k = k
        self.eps = eps

        if self.f0 < 0.0 or self.fc < 0.0 or self.k < 0.0:
            raise ValueError("Horton parameters must be non-negative.")
        if self.f0 < self.fc:
            raise ValueError("Expected f0 >= fc for classical Horton model.")
        
    def capacity(self, t):
        """
        Horton infiltration capacity (exponential decay law):
            f(t) = fc + (f0 - fc) * exp(-k * t)

        Parameters
        ----------
        t : float or np.ndarray
            Time value(s). 

        Returns
        -------
        float or np.ndarray
            Infiltration capacity at time t.
        """
        t = np.asarray(t, dtype=np.float64)
        fp = self.fc + (self.f0 - self.fc) * np.exp(-self.k * t)
        return np.float64(fp) if fp.ndim == 0 else fp
    
    def rate(self, t, rainfall, h=None, dt=None):
        """
        Actual infiltration rate used in the recharge source term.
        
        Simplest Choice: I(t) = min( f_p(t), rainfall, h / dt )

        Parameters
        ----------
        t : float or np.ndarray
            Current time.
        rainfall : float or np.ndarray
            Rainfall rate R.
        h : float or np.ndarray, optional  
            Local water depth.
        dt : float, optional
            Current timestep.
        
        Returns
        -------
        float or np.ndarray
            Actual infiltration rate at time t.
        """
        fp = np.asarray(self.capacity(t), dtype=np.float64)
        R = np.asarray(rainfall, dtype=np.float64)

        I = np.minimum(fp, R)

        if h is not None and dt is not None:
            h = np.asarray(h, dtype=np.float64)

            # Dry cell treatment via tolerance eps to avoid division by zero
            available_rate = np.maximum(h, 0.0) / max(float(dt), self.eps)
            
            I = np.minimum(I, available_rate)

        I = np.maximum(I, 0.0)
        return np.float64(I) if I.ndim == 0 else I
    
    def infiltrated_depth(self, t, rainfall, dt, h=None):
        """
        Water depth removed by infiltration over one timestep.

        Parameters
        ----------
        t : float or np.ndarray
            Current time.
        rainfall : float or np.ndarray
            Rainfall rate R.
        dt : float
            Current timestep.
        h : float or np.ndarray, optional
            Local water depth.

        Returns
        -------
        float or np.ndarray
            Infiltrated water depth over the timestep.
        """
        I = self.rate(t = t, rainfall=rainfall, h=h, dt=dt)
        depth = np.asarray(I, dtype=np.float64) * float(dt)
        return np.float64(depth) if depth.ndim == 0 else depth
    
    # String representation
    def __repr__(self):
        return (
            f"HortonInfiltration(f0={self.f0}, fc={self.fc},  "
            f"k={self.k}, eps={self.eps})"
        )
    

# Test Case & Sanity Check
if __name__ == "__main__":
    horton = HortonInfiltration(
        f0 = 1.977e-4,
        fc = 3.272e-5,
        k = 2.43e-5,
    )

    t = 0.0
    dt = 0.01
    R = 1.0e-4
    h = 1.0e-3

    fp = horton.capacity(t)
    I = horton.rate(t=t, rainfall=R, h=h, dt=dt)

    print("Horton Model:", horton)
    print(f"Capacity f_p(t) = {fp:.6e}")
    print(f"Actual Infiltration I(t) = {I:.6e}")