# Packages
import numpy as np
from recharge.context import SourceContext
from recharge.laws import HortonInfiltration
from recharge.source_terms import compute_recharge_source_n1

if __name__ == "__main__":
    # Example state U = [h, hu_m, ha_1]
    h = 1.0
    u_m = 0.75
    a1 = -0.25
    U = np.array([h, h * u_m, h * a1], dtype=np.float64)

    horton = HortonInfiltration(
        f0 = 1.977e-4,
        fc = 3.272e-5,
        k = 2.43e-3,
    )

    ctx = SourceContext(time = 0.0, dt = 0.01, cell_index = 0, x = 0.0)

    S, info = compute_recharge_source_n1(
        values = U,
        rainfall = 1.0e-3,
        infiltration_model = horton,
        context = ctx,
    )

    print("State U:", U)
    print("Diagnostics =", info)
    print("Recharge Source S =", S)