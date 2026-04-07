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

def compute_friction_matrix_n1(
        values,
        f_R,
        f_I,
        viscosity,
        slip_length,
        eps : float = 1e-14,
    ):
    """
    Compute the raw friction/mixing block P(U) for the N=1 SWME model with 
    rainfall and infiltration.

    State:
        U = [h, h * u_m, h * a_1]

    Source:
        P(U) = [
        0, 
        (f_R + f_I + nu/lambda) * u_m - (f_R - f_I - nu/lambda) * a_1,
        -3 (f_R - f_I - nu/lambda) * u_m + ( 3(f_R + f_I + nu/lambda) + 12 nu / h) * a_1
        ]

    Notes
    -----
    The solver should use S_total = S_{R, I}(U) - P(U)
    """
    U = np.asarray(values, dtype=np.float64)
    if U.shape[0] != 3:
        raise ValueError("N=1 SWME system expects state vector of length 3.")
    
    h = U[0]

    # Dry-cell treatment
    if h <= eps:
        return np.zeros(3, dtype=np.float64), {
            "h" : 0.0,
            "u_m" : 0.0,
            "a_1" : 0.0,
        }
    
    if slip_length <= eps:
        raise ValueError(
            "The current friction implementation uses nu / lambda explicitly,  " \
            "so it requires slip_length > 0."
        )
    
    # Set the other two state variables
    hu_m = U[1]
    ha_1 = U[2]

    u_m = hu_m / h
    a1 = ha_1 / h

    # Define the composite friction coefficients
    # There are of the form Navier Slip Law + Mixing Induced Friction

    # One is the total sum, the other looks like a skewness condition
    coeff = f_R + f_I + viscosity / slip_length
    skew = f_R - f_I - viscosity / slip_length

    # Build the composite friction vector as a numpy array
    P = np.array([
        0.0,
        coeff * u_m - skew * a1,
        -3 * skew * u_m + (3 * coeff + 12 * viscosity / h) * a1,
    ], dtype=np.float64)

    # Diagnostics for debugging
    diagnostics = {
        "h": h,
        "u_m": u_m,
        "a_1": a1,
        "f_R": f_R,
        "f_I": f_I,
    }

    return P, diagnostics

def compute_total_source(
        values,
        rainfall,
        infiltration_model,
        context,
        f_R,
        f_I,
        viscosity,
        slip_length,
        eps : float = 1e-14,
        ):
    """
    Compute the full N=1 source:
        S_total(U) = S_{R, I}(U) - P(U)
    """
    # First build the recharge mass production vector
    S, recharge_diagnostics = compute_recharge_source_n1(
        values = values,
        rainfall = rainfall,
        infiltration_model = infiltration_model,
        context = context,
        eps = eps,
    )

    # Then build the composite friction vector
    P, friction_diagnostics = compute_friction_matrix_n1(
        values = values,
        f_R = f_R,
        f_I = f_I,
        viscosity = viscosity,
        slip_length = slip_length,
        eps = eps,
    )

    # Combine diagnostics into one dictionary
    diagnostics = {
        "recharge": recharge_diagnostics,
        "friction": friction_diagnostics,
    }

    return S - P, diagnostics

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