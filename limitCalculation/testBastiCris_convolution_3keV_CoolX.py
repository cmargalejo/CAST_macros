import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import convolve2d
from scipy.ndimage import gaussian_filter

#def gaussian_kernel(size, sigma):
#    kernel = np.fromfunction(
#        lambda x, y: (1/(2*np.pi*sigma**2)) * np.exp(-((x - size//2)**2 + (y - size//2)**2) / (2*sigma**2)),
#        (size, size)
#    )
#    return kernel / np.sum(kernel)


def perform_interpolation(filename, csv_filename):
    # Read the text file using Pandas
    df = pd.read_csv(filename, skiprows = 2, delim_whitespace = True, names = ["x", "y", "z"])
    # Compute normalized data. Using *MEAN* of *DATA*
    df["x"] = df["x"] - df["x"].mean() + 1.3
    df["y"] = df["y"] - df["y"].mean() + 0.025
    # Sort data by *X* then *Y*
    df = df.sort_values(by = ["x", "y"], ascending = True)
    # Number of elements per axis
    size = int(np.sqrt(len(df)))
    # get unique elements
    unique_x = pd.unique(df["x"])
    unique_y = pd.unique(df["y"])
    # Construct array for interpolation
    zs = np.zeros([size, size])
    for x in range(size):
        for y in range(size):
            ## Important: First entries are `x`, second `y`. "Matrix" convention, i.e.
            ## first index is *row*, not column!
            zs[x, y] = df["z"][x * size + y] # Sorted by `x` first, hence `x` changes
                                             # *after* y.
    # Convolve the data with a Gaussian kernel
    #kernel_size = int(convolved_resolution / 0.1)
    #sigma = convolved_resolution / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to standard deviation
    #gaussian = gaussian_kernel(kernel_size, sigma)
    #convolved_zs = convolve2d(zs, gaussian, mode='same', boundary='wrap')
    sigma = convolved_resolution / 2.35482 #by definition, convolved_resolution = FWHM = 2 * sqrt(2 * ln(2)) * sigma = 2.35482 * sigma
    sigma = sigma / 100 #to convert into um
    print("Sigma = ", sigma)
    convolved_zs = gaussian_filter(zs, sigma=sigma, mode='wrap')
    # Create a RegularGridInterpolator
    interp_func = RegularGridInterpolator((unique_x, unique_y), convolved_zs, method='linear', bounds_error = False, fill_value = 0.0) #linear, nearest
    # Filter out CSV values that are within the bounds
    # x_min, y_min = unique_x.min(), unique_y.min()
    # x_max, y_max = unique_x.max(), unique_y.max()
    x_min, y_min = -10.0, -10.0
    x_max, y_max = 10.0, 10.0

    ## Number of points to use to evaluate interpolation
    NUM = 1000

    xs = np.linspace(x_min, x_max, NUM)
    ys = np.linspace(y_min, y_max, NUM)
    ## IMPORTANT: `ij` indexing means the indexing ends up being same as above, i.e.
    ## when calling `interp_func((X, Y))` it behaves as `zs[x, y]` and not `zs[y, x]`.
    (X, Y) = np.meshgrid(xs, ys, indexing = "ij")
    zzs = np.reshape(interp_func((X, Y)), [NUM, NUM])
    # Equivalent manual code:
    #zzs = np.zeros([NUM, NUM])
    #for yi, y in enumerate(ys):
    #    for xi, x in enumerate(xs):
    #        zzs[xi, yi] = interp_func([x, y])
    ## Create a plot.
    ## IMPORTANT: `lower` sets the coordinates for the *IMAGE* to start bottom
    ## left, not top as normal for *IMAGE* coordinates. This does not affect the actual *AXES*
    ## for plots put on top after!
    ## IMPORTANT 2: We *transpose* the given image, because `imshow` otherwise interprets
    ## our image the wrong way (columns and rows interchanged). The fact that we do everything
    ## correctly, we see in the interpolation and plotting of the points below.
    ## `vmax` is a heavy crop on the color scale. Otherwise goes to 7000.
    #plt.imshow(zzs.T, extent=(x_min, x_max, y_min, y_max), origin = "lower", vmax = 300)
    plt.imshow(zzs.T, extent=(x_min, x_max, y_min, y_max), origin = "lower")

    # Read positions (x, y) from the CSV file
    positions_data = np.genfromtxt(csv_filename, delimiter=',', skip_header=1, usecols=(0, 1))
    # Create specular image by negating x-coordinates
    specular_positions_data = np.copy(positions_data)
    specular_positions_data[:, 0] = -specular_positions_data[:, 0] #0 to negate x-coordinates (right-left) or 1 to negate y-coordinates (up-down)
    # Rotate the positions by 0 degrees clockwise (up to the right)
    angle_degrees = 45
    rotated_positions_data = np.dot(specular_positions_data, np.array([[np.cos(np.radians(angle_degrees)), -np.sin(np.radians(angle_degrees))],
                                                           [np.sin(np.radians(angle_degrees)), np.cos(np.radians(angle_degrees))]]))
    # Filter out CSV values that are within the bounds using numpy functions, because we are working with numpy arrays
    filtered_positions_data = rotated_positions_data[
        np.logical_and(
            np.logical_and(x_min <= rotated_positions_data[:, 0], rotated_positions_data[:, 0] <= x_max),
            np.logical_and(y_min <= rotated_positions_data[:, 1], rotated_positions_data[:, 1] <= y_max)

    )
    ]
    print("Filtered Positions:")
    print(filtered_positions_data)

    ## XXX: don't have `z` at the moment
    # Calculate contour levels for 95%, 90%, and 85% of the data
    zNonZero = df[df["z"] > 0.0]
    contour_levels = np.percentile(df["z"], [100-95,100-85,100-65]) #zs or z?
    contour_levels_no_zeroes = np.percentile(zNonZero["z"], [5,15,35]) #zs or z? [68,85,95] [100-95,100-85,100-65]
    print("contour levels = ", contour_levels_no_zeroes)
   
    #Print the percentiles used in a text box    
    percentile_values = [68, 85, 95]
    text_box_content = "\n".join([f"{percentile:.0f}%: {level:.0f}" for percentile, level in zip(percentile_values, contour_levels_no_zeroes)])

    text_box = plt.text(0.95, 0.95, text_box_content, transform=plt.gca().transAxes,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=dict(facecolor='white', alpha=0.5, edgecolor='white'))

    # Print the number of events used for interpolation
    num_events_used = len(filtered_positions_data)
    print(f"Number of events used for interpolation: {num_events_used}")

    # Plot the interpolated data

    ## XXX: we don't care about the original data
    #plt.imshow(zs, extent=(-10, 10, -10, 10), origin='lower', aspect='auto')

    plt.colorbar(label='Interpolated Value')
    plt.scatter(
        [pos[0] for pos in filtered_positions_data],
        [pos[1] for pos in filtered_positions_data],
        c='red', marker='x', label='Filtered CSV File Positions'
    )
    ## Interpolate each cluster and directly annotate
    for i, pos in enumerate(filtered_positions_data):
        interpolated_value = interp_func(pos)
        print(pos, " = ", interpolated_value)
        plt.annotate(f"{interpolated_value[0]:.2f}", pos, textcoords="offset points", xytext=(10,5), ha='center', fontsize=8, color='red')

    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Linear Interpolation using RegularGridInterpolator')
    plt.legend()

    # Plot contour lines for specified data percentages
    #plt.contour(X, Y, zzs, levels=contour_levels, colors='white')
    plt.contour(X, Y, zzs, levels=contour_levels_no_zeroes, colors='white')

    # Add the circle of calibration runs
    circle = plt.Circle((0, 0), 8.5, fill=False, edgecolor='white', linestyle='dashed', linewidth=1)
    plt.gca().add_artist(circle)

    plt.savefig("clusters_tests.pdf")
    plt.show()


# Call the function with your file names
#map_filename = '/home/cristina/GitHub/CAST_macros/limitCalculation/data/Jaime_data/2016_DEC_Final_CAST_XRT/3.00keV_2Dmap.txt'
map_filename = '/home/cristina/GitHub/CAST_macros/limitCalculation/data/Jaime_data/2016_DEC_Final_CAST_XRT/3.00keV_2Dmap_CoolX.txt'
csv_filename = 'data/cluster_candidates_tracking.csv'
convolved_resolution = 500  # Desired convolved spatial resolution (FWHM) in um (microns) is convolved_resolution FWHM = 2 * sqrt(2 * ln(2)) * sigma    
                            # 500 microns is equivalent to a physical resolution of 200 microns, Why? If FWHM=500 um, sigma=2.12*100 = 212 um. In other words, in the zs matrix, if sigma=2
                            # it means we use 2 indices in each direction to blur the image (better said, 2 indices are within the 1 sigma region, but it uses more indices).
                            # Each index is 0.1mm away so 2 indices is 0.2mm or 200 microns.
perform_interpolation(map_filename, csv_filename)

