import numpy as np
import matplotlib.pyplot as plt

def Saltykov(diameters,
             numbins=10,
             calc_vol=None,
             text_file=None,
             return_data=False,
             log_scale= True,
             left_edge=0):
    """ Estimate the actual (3D) distribution of grain size from the population
    of apparent diameters measured in a thin section using a Saltykov-type
    algorithm (Saltykov 1967; Sahagian and Proussevitch 1998).
    """

    if isinstance(numbins, int) is False:
        raise ValueError("Numbins must be a positive integer")
    if numbins <= 0:
        raise ValueError("Numbins must be a positive integer")
    if isinstance(left_edge, (int, float)):
        if left_edge < 0:
            raise ValueError("left_edge must be a positive scalar or 'min'")

    # set histogram left edge, either automatic or set by the user
    if left_edge == "min":
        minimo = diameters.min()
    else:
        minimo = left_edge

    # compute the histogram
    if log_scale:
        bin_edges = np.logspace(np.log10(min(diameters)), np.log10(max(diameters)), num=numbins+1)
        # FIX: Use counts (density=False), not density.
        # unfold_population expects counts and will normalize to density.
        freq, _ = np.histogram(diameters, bins=bin_edges)
    else:
        # FIX: Use counts (density=False), not density.
        freq, bin_edges = np.histogram(diameters, bins=numbins, range=(minimo, diameters.max()))

    # FIX: Create array of bin widths, as they are variable in log scale
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Unfold the population of apparent diameters using the Saltykov method
    # freq3D is a probability density function (PDF)
    freq3D = unfold_population(freq, bin_edges, bin_midpoints)

    # Calculate the volume-weighted cumulative frequency distribution
    # FIX: Pass bin_widths to correctly calculate volume per bin
    cdf_norm = volume_weighted_cdf(freq3D, bin_midpoints, bin_widths)

    # Estimate the volume of a particular grain size fraction (if apply)
    if calc_vol is not None:
        calc_volume_fraction_hist(calc_vol, bin_midpoints, cdf_norm)

    # Create a text file with the midpoints, class frequencies, and
    # cumulative volumes (if apply)
    if text_file is not None:
        # FIX: Pass bin_widths to correctly calculate freqs2one
        create_tabular_file(text_file, bin_widths, bin_midpoints, freq3D, cdf_norm)

    # return data or figure
    if return_data is True:
        return bin_midpoints,bin_edges, freq3D

    elif return_data is False:
        # FIX: Update print statement
        if log_scale:
            print("Using variable log-spaced bin widths.")
        else:
            print(f"calculated bin size = {bin_widths[0]:0.2f}")
        # FIX: Pass bin_widths to plot function
        return Saltykov_plot(bin_edges[:-1], freq3D, bin_widths, bin_midpoints, cdf_norm)

    else:
        raise TypeError("return_data must be set as True or False")

def unfold_population(freq, bin_edges, mid_points, normalize=True):
    """ Applies the Saltykov algorithm to unfold the population of apparent
    (2D) diameters into the actual (3D) population of grain sizes.
    """

    d_values = np.copy(bin_edges)
    midpoints = np.copy(mid_points)
    i = len(midpoints) - 1

    while i > 0:
        j = i
        D = d_values[-1]
        Pi = wicksell_solution(D, d_values[i], d_values[i + 1])

        if freq[i] > 0:
            while j > 0:
                D = midpoints[-1]
                Pj = wicksell_solution(D, d_values[j - 1], d_values[j])
                P_norm = (Pj * freq[i]) / Pi
                np.put(freq, j - 1, freq[j - 1] - P_norm)  # replace specified elements of an array
                j -= 1

            i -= 1
            d_values = np.delete(d_values, -1)
            midpoints = np.delete(midpoints, -1)

        else:  # if the value of the current class is zero or negative move to the next class
            i -= 1
            d_values = np.delete(d_values, -1)
            midpoints = np.delete(midpoints, -1)

    if normalize is True:
        freq = np.clip(freq, a_min=0.0, a_max=None)  # replacing negative values with zero
        # check if sum is zero to avoid division by zero
        if np.sum(freq) == 0:
             return freq # return all zeros
        freq_norm = freq / np.sum(freq)              # normalize to probability mass function (sums to 1)
        bin_widths = bin_edges[1:] - bin_edges[:-1]  # get widths of each bin
        freq_norm = freq_norm / bin_widths           # normalize to probability density function (integrates to 1)
        return freq_norm

    else:
        return freq

def wicksell_solution(diameter, lower_bound, upper_bound):
    """ Estimate the cross-section size probability for a discretized population
    of spheres.
    """

    radius = diameter / 2
    r1, r2 = lower_bound / 2, upper_bound / 2
    return 1 / radius * (np.sqrt(radius**2 - r1**2) - np.sqrt(radius**2 - r2**2))

def volume_weighted_cdf(freqs, bin_midpoints, bin_widths):
    """Calculates the volume-weighted cumulative frequency
    distribution of a histogram in percentage.
    """

    # Calculate the volume for each bin assuming that
    # particles are spherical and sizes distributed homogeneously
    grain_volumes = diameter_to_volume(bin_midpoints)

    # FIX: Calculate probability mass (area) for each bin
    prob_per_bin = freqs * bin_widths

    # FIX: Weight the probability mass (not density) by the volume
    vol_weighted_freqs = prob_per_bin * grain_volumes

    # Compute the cumulative sum of the volume-weighted counts
    cumulative_volume = np.cumsum(vol_weighted_freqs)

    # Normalize the cumulative sum to get the cumulative frequency distribution
    if cumulative_volume[-1] == 0:
        return np.zeros_like(cumulative_volume) # Avoid division by zero
    vol_cfd = cumulative_volume / cumulative_volume[-1]

    return 100 * vol_cfd

def calc_volume_fraction_hist(size, bin_midpoints, cdf_norm):
    """Calculates and print the volume fraction of a
    occuppied up to a grain size specified by the user
    """
    
    # Handle edge case where size is smaller than the first midpoint
    if size < bin_midpoints[0]:
        print("=================================================")
        print(f"volume fraction (up to {size} microns) = 0.00 %")
        print("=================================================")
        return None
    
    # Handle edge case where size is larger than the last midpoint
    if size >= bin_midpoints[-1]:
        print("=================================================")
        print(f"volume fraction (up to {size} microns) = 100.00 %")
        print("=================================================")
        return None

    index = np.argmax(bin_midpoints > size)
    
    # Avoid division by zero if midpoints are identical
    if (bin_midpoints[index] - bin_midpoints[index - 1]) == 0:
        volume = cdf_norm[index - 1]
    else:
        angle = np.arctan(
            (cdf_norm[index] - cdf_norm[index - 1])
            / (bin_midpoints[index] - bin_midpoints[index - 1])
        )
        volume = cdf_norm[index - 1] + np.tan(angle) * (size - bin_midpoints[index - 1])

    if volume < 100.0:
        print("=================================================")
        print(f"volume fraction (up to {size} microns) = {volume:.2f} %")
        print("=================================================")
    else:
        print("=================================================")
        print(f"volume fraction (up to {size} microns) = 100.00 %")
        print("=================================================")

    return None

def diameter_to_volume(d):
    """Calculates the volume of a sphere with diameter d"""
    return (np.pi / 6) * d**3

def Saltykov_plot(left_edges, freq3D, bin_widths, mid_points, cdf_norm):
    """ Generate two plots once the Saltykov method is applied:
    """

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

    # frequency vs grain size plot
    ax1.bar(
        left_edges,
        freq3D,
        # FIX: Use the array of bin_widths
        width=bin_widths,
        color="xkcd:azure",
        edgecolor="#d9d9d9",
        align="edge",
    )
    ax1.set_ylabel("density", fontsize=18)
    ax1.set_xlabel(r"diameter ($\mu m$)", fontsize=18)

    # volume-weighted cumulative frequency curve
    ax2.set_ylim([-2, 105])
    ax2.plot(
        mid_points,
        cdf_norm,
        "o-",
        color="#ed4256",
        label="volume weighted CFD",
        linewidth=2,
    )
    ax2.set_ylabel("cumulative volume (%)", color="#252525")
    ax2.set_xlabel(r"diameter ($\mu m$)", color="#252525")

    fig.tight_layout()

    return fig, (ax1, ax2)

def create_tabular_file(text_file, bin_widths, bin_midpoints, freq3D, cdf_norm):
    """Generate and save a tabular file with data
    """

    if isinstance(text_file, str) is False:
        raise TypeError("text_file must be a string type")

    from pandas import DataFrame

    df = DataFrame(
        {
            "bin_midpoints": np.around(bin_midpoints, 3),
            "freqs": np.around(freq3D, 4),
            # FIX: Calculate probability (area) using the bin_widths array
            "freqs2one": np.around(freq3D * bin_widths, 3),
            "cum_vol": np.around(cdf_norm, 2),
        }
    )

    if text_file.endswith((".tsv", ".txt")):
        df.to_csv(text_file, sep="\t", index=False)
    elif text_file.endswith(".csv"):
        df.to_csv(text_file, sep=";", index=False)
    else:
        raise ValueError("text file must be specified as .csv, .tsv or .txt")

    print("=======================================")
    print(f"The file {text_file} was created")
    print("=======================================")

    return None
