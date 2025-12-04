import numpy as np
from scipy.optimize import nnls

def k(a, r):
    """
    Auxiliary function for geometric calculations.
    """
    return np.maximum(a, r) + np.sqrt(np.maximum(a, r)**2 - r**2)

def phi(a_list, i, r, cutoff=1):
    """
    Basis function calculation.
    """
    if cutoff == 1:
        ans = r * np.log(k(a_list[i+1], r) / k(a_list[i], r))
        return ans
    else:
        step_function = np.where(r >= (1 - cutoff) * r[i], 1, 0)
        ans = step_function * r * np.log(k(a_list[i+1], r) / k(a_list[i], r))
        return ans

def get_psi(i, j, a_list, r, cutoff=1):
    """
    Calculate Psi matrix element.
    """
    ans = np.sum(phi(a_list, i, r, cutoff) * phi(a_list, j, r, cutoff)) * (r[1] - r[0])
    return ans

def get_b(i, a_list, r, f_r, cutoff=1):
    """
    Calculate b vector element.
    """
    ans = np.sum(phi(a_list, i, r, cutoff) * f_r) * (r[1] - r[0])
    return ans

def estimate_c(csd, r_min=0.5, cutoff=0.3, num_bins=20, linear=True):
    """
    Estimate the 3D pore radius distribution from chord size distribution.

    Parameters
    ----------
    csd : array_like
        Chord size distribution data.
    r_min : float, optional
        Minimum radius to consider.
    cutoff : float, optional
        Cutoff parameter for the phi function.
    num_bins : int, optional
        Number of bins for the histogram.
    linear : bool, optional
        If True, use linear binning. If False, use logarithmic binning.

    Returns
    -------
    c : array
        Estimated coefficients (probability density).
    bins : array
        Bin edges.
    residual : float
        Residual of the NNLS solution.
    """
    if linear:
        bins = np.linspace(max(min(csd), r_min), max(csd), num=num_bins+1)
    else:
        bins = np.logspace(np.log10(max(min(csd), r_min)), np.log10(max(csd)), num=num_bins+1)
    
    hist_log, _ = np.histogram(csd, bins=bins)
    # Avoid division by zero if sum is 0
    total_hist = sum(hist_log)
    if total_hist == 0:
        hist_log_dens = np.zeros_like(hist_log, dtype=float)
    else:
        hist_log_dens = hist_log / (bins[1:] - bins[:-1]) / total_hist

    dr = 0.001
    r = np.arange(max(min(csd), r_min), max(csd), dr)
    
    # Get step-wise P_c probability density function
    f_r = []
    current_step = 0
    for i, r_i in enumerate(r):
        if r_i < bins[-1] and r_i > bins[0]:
            if current_step < len(hist_log_dens) - 1:
                if r_i < bins[current_step+1]:
                    f_r.append(hist_log_dens[current_step])
                else:
                    current_step = current_step + 1
                    f_r.append(hist_log_dens[current_step])
            else:
                f_r.append(hist_log_dens[current_step])
        else:
            f_r.append(0)
    f_r = np.array(f_r)
    
    # Calculate Psi and b
    psi = np.zeros((num_bins, num_bins))
    b = np.zeros(num_bins)
    for i in range(num_bins):
        b[i] = get_b(i, bins, r, f_r, cutoff=cutoff)
        for j in range(num_bins):
            psi[i, j] = get_psi(i, j, bins, r, cutoff=cutoff)
            
    alpha, residual = nnls(psi, b)
    
    # Avoid division by zero
    denom = sum(alpha * (bins[1:] - bins[:-1]))
    if denom == 0:
        c = alpha # Or zeros?
    else:
        c = alpha / denom
        
    return c, bins, residual
