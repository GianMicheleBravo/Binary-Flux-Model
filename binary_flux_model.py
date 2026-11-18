import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as apc

plt.rcParams.update({'font.size': 15})

Planck_law_term_1 = (8 * np.pi * apc.h * apc.c).to_value(u.Jy * u.nm**3)
Planck_law_term_2 = (apc.h * apc.c / apc.k_B).to_value(u.nm * u.K)


def Planck_law(lambd_nm, T_K):
    """Return the flux from blackbody emission.    
    Input:
        lambd_nm: the wavelength at which to calculate the flux, in nm, larger than 0
        T_K: the temperature, in Kelvin, larger than to 0

    Output:
        flux: the flux of the object at lambd_nm, in Jy/sr, larger than 0

    Example:
        >>> Planck_law(550, 3000)
        386.54307053235556"""

    return Planck_law_term_1 / (lambd_nm**3 * (np.exp(Planck_law_term_2 / (lambd_nm * T_K)) - 1))


def parameter_subset(parameter, subset):
    """If the input parameter is an array, rturn the indicated subset of that array.
    If the input parameter is scalar, return the scalar itself, as the scalar parameter is 
    assumed to be equal everywhere, incuding the subset. Use to allow for both scalar and
    array input.

    Input:
        parameter: array-like or scalar
        subset: index, slice, or other object that selects a subset of an array

    Output:
        The array subset, or the scalar itself if the input parameter is scalar

    Example:
        >>> parameter_subset([1, 2, 3], 0)
        1
        >>> parameter_subset(4, 0)
        4
        >>> arr = np.array([1, 2, 3])
        >>> subset = np.array([True, False, True])
        >>> parameter_subset(arr, subset)
        [1 3]
        >>> parameter_subset(4, subset)
        4"""

    if hasattr(parameter, '__iter__'):
        return parameter[subset]
    return parameter


def phase(t, P, t0):
    out = t - t0
    out /= P
    out = np.mod(out, 1)

    return out


def normalized_separation(phi, i):
    phi_2_pi = phi * np.pi * 2

    out = np.sin(phi_2_pi)**2
    out += (np.cos(i) * np.cos(phi_2_pi))**2
    out = np.sqrt(out)

    return out


def projected_xy_pos(phi, i, theta_max):
    phi_2_pi = phi * 2 * np.pi

    x = np.sin(phi_2_pi) * theta_max / 2
    y = np.cos(phi_2_pi) * theta_max * np.cos(i) / 2

    return x, y


def partial_circle_intersection_area(r1, r2, d):
    """Return the intersection area of two circles that partially cover each other.
    Input:
        r1: the radius of the first circle
        r2: the radius of the second circle, same unit as r1
        d: the distance between the centers of the two circles, same unit as r1

    Output:
        out: the intersect area, same unit as r1^2

    From https://mathworld.wolfram.com/Circle-CircleIntersection.html"""

    r1_squared = r1**2
    r2_squared = r2**2
    d_squared = d**2

    out = r1_squared * \
        np.arccos((d_squared + r1_squared - r2_squared) / (2 * d * r1))
    out += r2_squared * \
        np.arccos((d_squared - r1_squared + r2_squared) / (2 * d * r2))
    out -= 0.5 * np.sqrt((r1 + r2 + d) * (r1 + r2 - d)
                         * (d + r1 - r2) * (d - r1 + r2))

    return out


def input_sanity_check(lambd, P, i_deg, f1, f2):
    if np.any(lambd) <= 0:
        raise ValueError('The wavelength lambd must be larger than 0.')
    if np.any(P) <= 0:
        raise ValueError('The period P must be larger than 0.')
    if np.any(i_deg < 0):
        raise ValueError(
            'The inclination i_deg must be between 0 and 90 degrees. It is currently too small.')
    if np.any(i_deg > 90):
        raise ValueError(
            'The inclination i_deg must be between 0 and 90 degrees. It is currently too large.')
    if np.any(f1 < 0):
        raise ValueError(
            'The primary fractional radius f1 must be between 0 and 1. It is currently too small.')
    if np.any(f1 > 1):
        raise ValueError(
            'The primary fractional radius f1 must be between 0 and 1. It is currently too large.')
    if np.any(f2 < 0):
        raise ValueError(
            'The secondary fractional radius f2 must be between 0 and 1. It is currently too small.')
    if np.any(f2 > 1):
        raise ValueError(
            'The secondary fractional radius f2 must be between 0 and 1. It is currently too large.')
    if np.any(f1 + f2) > 1:
        raise ValueError(
            'The sum of the fractional radii f1+f2 must be smaller than 1. The two objects would collide otherwise.')


def eclipsing_binary_flux(t, lambd, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1=Planck_law, spectrum_function_2=Planck_law):
    """Return the light emitted by an eclipsing binary. 
        This model assumes the objects as uniform, opaque spheres in a circular orbit.
    Input:
        t: time, array-like or scalar
        lambd: wavelength, array-like with or scalar, larger than 0, default in nm
        P: period, array-like with or scalar, larger than 0, same unit as t
        t0: time of primary eclipse, array-like with or scalar, same unit as t
        theta_max: maximum projected separation angle of the binary, array-like with or scalar,
            default in radians
        i_deg: angle between the plane vector of the system and the line of sight,
            array-like or scalar, between 0 and 90, must be in degrees
        f1: the ratio between the primary radius and the binary distance, array-like or scalar,
            between 0 and 1
        f2: the ratio between the secondary radius and the binary distance,array-like or scalar,
            between 0 and 1
        T1: temperature of the primary, array-like or scalar, default in K
        T2: temperature of the secondary, array-like or scalar, default in K

        spectrum_function_1: the function used to get flux/solid angle of the primary:
            spectrum_function_1(lambd, T) -> flux
            Default: Planck_law:
                lambd_nm, T_K -> flux in Jy/sr

        spectrum_function_2: the function used to get flux/solid angle of the secondary:
            spectrum_function_2(lambd, T) -> flux
            Default: Planck_law:
                lambd_nm, T_K -> flux in Jy/sr

    Output:
        out: the flux of the system, unit specificed by
            spectrum_function_1 and spectrum_function_2, default in Jy  

    Note on units:
        The inclination angle i_deg must be in degrees.
        All time input must have the same unit.

        All other parameters (lambd, theta_max, T1 and T2) have optional units, determined
        by what spectrum_function is used. The default Planck_law needs the following input units:
            lambd: nm
            thema_max: radians
            T1: K
            T2: K
        With Planck_law and these units, the output flux is in Jy.

        If different spectrum_functions are used, different input and output units can be used.

    Example:
        >>> t = [0, 10, 20, 30, 40]
        >>> lambd_nm = 550 # green
        >>> P = 40
        >>> t0 = 10
        >>> i_deg = 85
        >>> theta_max_sr = 1/21600000 # 1 milliarcsecond
        >>> f1 = 0.3
        >>> f2 = 0.4
        >>> T1_K = 6000
        >>> T2_K = 5000
        >>> eclipsing_binary_flux(time, lambd_nm, P, t0, i_deg, theta_max_sr, f1, f2, T1_K, T2_K)
        array([40906.9806891 , 17367.02548793, 40906.9806891 , 31138.02885213,
            40906.9806891 ])"""

    # input formatting and testing
    t = np.array(t, dtype=np.float64)
    lambd = np.array(lambd, dtype=np.float64)
    P = np.array(P, dtype=np.float64)
    i_deg = np.array(i_deg, dtype=np.float64)
    theta_max = np.array(theta_max, dtype=np.float64)
    f1 = np.array(f1, dtype=np.float64)
    f2 = np.array(f2, dtype=np.float64)
    T1 = np.array(T1, dtype=np.float64)
    T2 = np.array(T2, dtype=np.float64)

    input_sanity_check(lambd, P, i_deg, f1, f2)

    # output formatting

    out_shape = np.broadcast_shapes(
        t.shape, lambd.shape, P.shape, i_deg.shape, theta_max.shape, f1.shape, f2.shape, T1.shape, T2.shape)
    out = np.zeros(out_shape)  # the output array

    # the base flux
    theta_1 = f1 * theta_max
    theta_2 = f2 * theta_max

    alpha_1 = theta_1**2 * np.pi
    alpha_2 = theta_2**2 * np.pi

    surface_luminosity_1 = spectrum_function_1(lambd, T1)
    surface_luminosity_2 = spectrum_function_2(lambd, T2)

    flux1 = surface_luminosity_1 * alpha_1
    flux2 = surface_luminosity_2 * alpha_2

    out += flux1 + flux2

    # testing whether the system eclipses

    i = i_deg * np.pi / 180  # in radians
    phi = phase(t, P, t0)
    radius_sum = f1 + f2

    sep = normalized_separation(phi, i)

    if np.any(sep < radius_sum):  # if we have any eclipse
        radius_difference = np.abs(f1 - f2)

        star_1_behind = np.bitwise_or(phi <= 0.25, phi >= 0.75)
        star_2_behind = np.bitwise_not(star_1_behind)

        full_eclipse = sep <= radius_difference
        partial_eclipse = np.bitwise_and(
            radius_difference < sep, sep < radius_sum)

        if np.any(full_eclipse):
            alpha_1_broadcasted, alpha_2_broadcasted = np.broadcast_arrays(
                alpha_1, alpha_2)

            min_alpha = np.min(
                (alpha_1_broadcasted, alpha_2_broadcasted), axis=0)

            full_eclipse_1 = np.bitwise_and(full_eclipse, star_1_behind)
            full_eclipse_2 = np.bitwise_and(full_eclipse, star_2_behind)

            if np.any(full_eclipse_1):
                out[full_eclipse_1] -= parameter_subset(
                    surface_luminosity_1, full_eclipse_1) * parameter_subset(min_alpha, full_eclipse_1)

            if np.any(full_eclipse_2):
                out[full_eclipse_2] -= parameter_subset(
                    surface_luminosity_2, full_eclipse_2) * parameter_subset(min_alpha, full_eclipse_2)

        if np.any(partial_eclipse):
            partial_eclipse_1 = np.bitwise_and(partial_eclipse, star_1_behind)
            partial_eclipse_2 = np.bitwise_and(partial_eclipse, star_2_behind)

            theta_sep = sep * theta_max

            if np.any(partial_eclipse_1):
                out[partial_eclipse_1] -= partial_circle_intersection_area(
                    parameter_subset(theta_1, partial_eclipse_1), parameter_subset(theta_2, partial_eclipse_1), parameter_subset(theta_sep, partial_eclipse_1)) * parameter_subset(surface_luminosity_1, partial_eclipse_1)

            if np.any(partial_eclipse_2):
                out[partial_eclipse_2] -= partial_circle_intersection_area(
                    parameter_subset(theta_1, partial_eclipse_2), parameter_subset(theta_2, partial_eclipse_2), parameter_subset(theta_sep, partial_eclipse_2)) * parameter_subset(surface_luminosity_2, partial_eclipse_2)

    return out


def isoscelese_triangle_angle(f):
    """Return the unique angle in an isosceles triangle.

    Input:
        f: the ratio beween the unique side and the two equal length side

    Output:
        angle: the  unique angle in an isosceles triangle, in radians"""

    angle = np.arcsin(f / (2) * np.sqrt(4 - f**2))

    return angle


def plot_red_green_blue_geometry(ax, t, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1=Planck_law, spectrum_function_2=Planck_law, lambd_red=700, lambd_green=550, lambd_blue=400, colorful_orbit=False):
    i = i_deg * np.pi / 180

    phi_t = phase(t, P, t0)

    r1 = theta_max * f1
    r2 = theta_max * f2

    rmax = np.max((r1, r2))

    margin = 0.1

    xlimit = (rmax + theta_max / 2) * (1 + margin)
    ylimit = (rmax + theta_max * np.cos(i) / 2) + (xlimit * margin)

    ax.set_aspect('equal')
    ax.set(xlim=(-xlimit, xlimit), ylim=(-ylimit, ylimit))
    ax.set_xlabel('Projected x angle in radians')
    ax.set_ylabel('Projected y angle in radians')
    ax.grid(zorder=-9)

    red1 = spectrum_function_1(lambd_red, T1)
    green1 = spectrum_function_1(lambd_green, T1)
    blue1 = spectrum_function_1(lambd_blue, T1)

    red2 = spectrum_function_2(lambd_red, T2)
    green2 = spectrum_function_2(lambd_green, T2)
    blue2 = spectrum_function_2(lambd_blue, T2)

    colors = np.array((red1, green1, blue1, red2, green2, blue2))
    colors /= np.max(colors)

    color1 = colors[:3]
    color2 = colors[3:]

    phi_t = phase(t, P, t0)

    x1, y1 = projected_xy_pos(phi_t, i, theta_max)
    x2 = -x1
    y2 = -y1

    r1 = theta_max * f1
    r2 = theta_max * f2

    star1 = plt.Circle((x1, y1), r1, color=color1, label='primary')
    star2 = plt.Circle((x2, y2), r2, color=color2, label='secondary')

    outline1 = plt.Circle((x1, y1), r1, color='black', alpha=0.5, fill=False)
    outline2 = plt.Circle((x2, y2), r2, color='black', alpha=0.5, fill=False)

    phase_angle1 = isoscelese_triangle_angle(f1 * 2) / (2 * np.pi)
    phase_angle2 = isoscelese_triangle_angle(f2 * 2) / (2 * np.pi)

    phi_span_trailing_1 = np.linspace(phi_t - 0.25, phi_t - phase_angle1)
    phi_span_leading_1 = np.linspace(phi_t + phase_angle1, phi_t + 0.25)

    phi_span_trailing_2 = np.linspace(phi_t + 0.25, phi_t + 0.5 - phase_angle2)
    phi_span_leading_2 = np.linspace(phi_t + 0.5 + phase_angle2, phi_t + 0.75)

    x_leading_1, y_leading_1 = projected_xy_pos(
        phi_span_trailing_1, i, theta_max)
    x_trailing_1, y_trailing_1 = projected_xy_pos(
        phi_span_leading_1, i, theta_max)

    x_trailing_2, y_trailing_2 = projected_xy_pos(
        phi_span_trailing_2, i, theta_max)
    x_leading_2, y_leading_2 = projected_xy_pos(
        phi_span_leading_2, i, theta_max)

    z_trailing_1 = np.cos((phi_t + phase_angle1) * 2 * np.pi)
    z1 = np.cos(phi_t * 2 * np.pi)
    z_leading_1 = np.cos((phi_t - phase_angle1) * 2 * np.pi)

    z_trailing_2 = np.cos((phi_t + 0.5 - phase_angle2) * 2 * np.pi)
    z2 = np.cos((phi_t + 0.5) * 2 * np.pi)
    z_leading_2 = np.cos((phi_t + 0.5 + phase_angle2) * 2 * np.pi)

    if colorful_orbit:
        orbit_colors = ['yellow', 'red', 'green', 'blue']
    else:
        orbit_colors = ['gray', 'gray', 'gray', 'gray']

    if z1 >= z2:  # star 1 further
        # star 1 drawn first

        if z_leading_1 >= z1:
            ax.plot(x_leading_1, y_leading_1,
                    color=orbit_colors[0], zorder=1, linewidth=4)
        if z_trailing_1 >= z1:
            ax.plot(x_trailing_1, y_trailing_1,
                    color=orbit_colors[1], zorder=1, linewidth=4)

        ax.add_artist(star1)
        ax.add_artist(outline1)

        if z_leading_1 < z1:
            ax.plot(x_leading_1, y_leading_1,
                    color=orbit_colors[0], zorder=1, linewidth=4)
        if z_trailing_1 < z1:
            ax.plot(x_trailing_1, y_trailing_1,
                    color=orbit_colors[1], zorder=1, linewidth=4)

        # star 2 drawn on top

        if z_leading_2 >= z2:
            ax.plot(x_leading_2, y_leading_2,
                    color=orbit_colors[2], zorder=1, linewidth=4)
        if z_trailing_2 >= z2:
            ax.plot(x_trailing_2, y_trailing_2,
                    color=orbit_colors[3], zorder=1, linewidth=4)

        ax.add_artist(star2)
        ax.add_artist(outline2)

        if z_leading_2 < z2:
            ax.plot(x_leading_2, y_leading_2,
                    color=orbit_colors[2], zorder=1, linewidth=4)
        if z_trailing_2 < z2:
            ax.plot(x_trailing_2, y_trailing_2,
                    color=orbit_colors[3], zorder=1, linewidth=4)

    else:
        # star 2 drawn first

        if z_leading_2 >= z2:
            ax.plot(x_leading_2, y_leading_2,
                    color=orbit_colors[2], zorder=1, linewidth=4)
        if z_trailing_2 >= z2:
            ax.plot(x_trailing_2, y_trailing_2,
                    color=orbit_colors[3], zorder=1, linewidth=4)

        ax.add_artist(star2)
        ax.add_artist(outline2)

        if z_leading_2 < z2:
            ax.plot(x_leading_2, y_leading_2,
                    color=orbit_colors[2], zorder=1, linewidth=4)
        if z_trailing_2 < z2:
            ax.plot(x_trailing_2, y_trailing_2,
                    color=orbit_colors[3], zorder=1, linewidth=4)

        # star 1 drawn on top

        if z_leading_1 >= z1:
            ax.plot(x_leading_1, y_leading_1,
                    color=orbit_colors[0], zorder=1, linewidth=4)
        if z_trailing_1 >= z1:
            ax.plot(x_trailing_1, y_trailing_1,
                    color=orbit_colors[1], zorder=1, linewidth=4)

        ax.add_artist(star1)
        ax.add_artist(outline1)

        if z_leading_1 < z1:
            ax.plot(x_leading_1, y_leading_1,
                    color=orbit_colors[0], zorder=1, linewidth=4)
        if z_trailing_1 < z1:
            ax.plot(x_trailing_1, y_trailing_1,
                    color=orbit_colors[1], zorder=1, linewidth=4)


def plot_red_green_blue_light_curve(ax, t, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1=Planck_law, spectrum_function_2=Planck_law, lambd_red=700, lambd_green=550, lambd_blue=400):
    n = np.floor((t - t0) / P + 0.25)

    t_min = (n - 0.25) * P + t0
    t_max = t_min + P

    t_span = np.linspace(t_min, t_max, 1000)

    flux_red = eclipsing_binary_flux(
        t_span, lambd_red, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1, spectrum_function_2)
    flux_green = eclipsing_binary_flux(
        t_span, lambd_green, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1, spectrum_function_2)
    flux_blue = eclipsing_binary_flux(
        t_span, lambd_blue, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1, spectrum_function_2)

    flux_red_t = np.interp(t, t_span, flux_red)
    flux_green_t = np.interp(t, t_span, flux_green)
    flux_blue_t = np.interp(t, t_span, flux_blue)

    flux_all = np.hstack((flux_red, flux_green, flux_blue))

    flux_min = np.min(flux_all)
    flux_max = np.max(flux_all)

    ax.grid(zorder=-9)
    ax.vlines([t], flux_min, flux_max, color='gray', zorder=-1)

    ax.plot(t_span, flux_red, color='red',
            label=f'flux at {round(lambd_red)} nm')
    ax.plot(t_span, flux_green, color='green',
            label=f'flux at {round(lambd_green)} nm')
    ax.plot(t_span, flux_blue, color='blue',
            label=f'flux at {round(lambd_blue)} nm')

    ax.scatter(t, flux_red_t, color='red')
    ax.scatter(t, flux_green_t, color='green')
    ax.scatter(t, flux_blue_t, color='blue')

    ax.legend(loc='upper right')
    ax.set_xlabel('Time')
    ax.set_ylabel('Flux in Jy')


def plot_eclipsing_binary_flux(t, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1=Planck_law, spectrum_function_2=Planck_law, lambd_red=700, lambd_green=550, lambd_blue=400, show=True):
    """Plots a binary system with given parameters at a given time.
    Input:
        t: the time at which to plot the system
        P: period, larger than 0, same unit as t
        t0: time of primary eclipse, same unit as t
        theta_max: maximum projected separation angle of the binary, default in radians
        i_deg: angle between the plane of the system and the line of sight, between 0 and 90 degrees
        f1: the ratio between the primary radius and the binary distance, between 0 and 1
        f2: the ratio between the secondary radius and the binary distance, between 0 and 
        T1: temperature of the primary, default in K
        T2: temperature of the secondary, default in K

        spectrum_function_1: the function used to get flux of the primary:
            spectrum_function_1(lambd, T) -> flux
            Default: Planck_law:
                lambd_nm, T_K -> flux in  W/(m^2 * nm * sr)

        spectrum_function_2: the function used to get flux of the secondary:
            spectrum_function_1(lambd, T) -> flux
            Default: Planck_law:
                lambd_nm, T_K -> flux in W/(m^2 * nm * sr)

        lambd_red: the wavelength to use for the plot's red colour channel.
            Defaut: 700 
        lambd_green: the wavelength to use for the plot's green colour channel.
            Defaut: 550
        lambd_red: the wavelength to use for the plot's red colour channel.
            Defaut: 400

        show: whether to display the plot or return fig, axs
            Default: True"""

    fig, axs = plt.subplots(2, figsize=[10, 10])
    ax_img, ax_lc = axs

    # the geometry plot

    plot_red_green_blue_geometry(ax_img, t, P, t0, i_deg, theta_max, f1, f2, T1, T2,
                                 spectrum_function_1, spectrum_function_2, lambd_red, lambd_green, lambd_blue)

    # the light curve plot
    plot_red_green_blue_light_curve(ax_lc, t, P, t0, i_deg, theta_max, f1, f2, T1, T2,
                                    spectrum_function_1, spectrum_function_2, lambd_red, lambd_green, lambd_blue)

    if show:
        plt.show()
    else:
        return fig, axs
