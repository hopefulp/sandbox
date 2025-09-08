import numpy as np

def generate_distances(r_min, r_eq, r_max, n_short=3, n_long=20):
    """
    Generate non-uniform distance list for diatomic molecule scan.
    - r_min: minimum distance (Å)
    - r_eq: equilibrium distance (Å)
    - r_max: maximum distance (Å)
    - n_short: number of points between r_min and r_eq
    - n_long: number of points between r_eq and r_max
    Returns: sorted numpy array of distances
    """

    # 1. Short side (r < r_eq): cluster closer to r_eq
    short_side = r_eq - (r_eq - r_min) * (np.linspace(0, 1, n_short) ** 2)
    #print(f"short side: {short_side}")

    # 2. Long side (r > r_eq): logarithmic spacing
    long_side = r_eq + (r_max - r_eq) * np.logspace(-2, 0, n_long, base=10)
    #print(f"long side: {long_side}")
    # Remove possible duplicate at r_eq
    distances = np.unique(np.concatenate([short_side, [r_eq], long_side]))

    return np.round(distances, 4)
