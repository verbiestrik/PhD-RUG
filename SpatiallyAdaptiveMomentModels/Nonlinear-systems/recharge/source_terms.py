# Packages
import numpy as np

def compute_recharge_source_n1(
        values,
        rainfall,
        infiltration_model,
        context,
        eps : float = 1e-14,
    ):
    """
    Compute the pure rainfall-infiltration source term S_{R, I}(U) for the 
    first order SWME system.

    State:
        U = [h, h * u_m, h * a_1]

    Source:
        S_{R, I}(U) = [
            R - I, 
            (R - I) * u_m - (R + I) * a_1,
            2 * (R - I) * a_1
        ]

    Parameters
    ----------
    values : array_like
        Conserved state variables [h, h * u_m, h * a_1].
    rainfall : float or callable
        Rainfall rate R. If callable, it is evaluated as rainfall(context).
    infiltration_model : object
        A customly encoded infiltration model. Must provide a method
            rate(t, rainfall, h, dt) or similar
    context : SourceContext
        Runtime source context.
    eps : float
        Small tolerance to avoid division by zero.

    Returns
    -------
    np.ndarray
        Recharge source vector of shape (3,).
    dict
        Diagnostics: R, I, h, u_m, a_1
    """
    U = np.asarray(values, dtype=np.float64)
    if U.shape[0] != 3:
        raise ValueError("N=1 SWME system expects state vector of length 3.")
    
    # Manually set the height as first state variable
    h = U[0]

    # Dry-cell treatment
    if h <= eps:
        return np.zeros(3, dtype=np.float64), {
            "R": 0.0,
            "I": 0.0,
            "h": h,
            "u_m": 0.0,
            "a_1": 0.0,
        }

    # Set the other two state variables
    hu_m = U[1]
    ha_1 = U[2]

    u_m = hu_m / h
    a1 = ha_1 / h

    # Rainfall can be a constant or a function of context
    R = rainfall(context) if callable(rainfall) else np.float64(rainfall)

    # Define the infiltration model
    I = infiltration_model.rate(
        t = context.time,
        rainfall = R,
        h = h,
        dt = context.dt,
    )

    # Build rainfall-infiltration mass production source term
    S = np.array([
        R - I,
        (R - I) * u_m - (R + I) * a1,
        2 * (R - I) * a1,
    ], dtype=np.float64)

    # Set up diagnostics for debugging
    diagnostics = {
        "R": R,
        "I": I,
        "h": h,
        "u_m": u_m,
        "a_1": a1,
    }

    # Output
    return S, diagnostics

def compute_recharge_source(
        order,
        values,
        rainfall,
        infiltration_model,
        context,
        eps : float = 1e-12
    ):
    """
    Dispatch recharge source by moment order.
    """
    if order == 1:
        return compute_recharge_source_n1(
            values = values,
            rainfall = rainfall,
            infiltration_model = infiltration_model,
            context = context,
            eps = eps,
        )
    else:
        raise NotImplementedError(
            f"Recharge source currently implemented only for N=1."
        )