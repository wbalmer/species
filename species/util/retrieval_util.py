"""
Utility functions for atmospheric retrieval with ``petitRADTRANS``.
This module was put together many contributions by Paul Mollière
(MPIA).
"""

import copy
import sys
import warnings

from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import interp1d, PchipInterpolator
from scipy.ndimage import gaussian_filter
from typeguard import typechecked

from species.core import constants


@typechecked
def get_line_species() -> List[str]:
    """
    Function to get the list of the molecular and atomic line species.
    This function is not used anywhere so could be removed.

    Returns
    -------
    list(str)
        List with the line species.
    """

    return [
        "CH4",
        "CO",
        "CO_all_iso",
        "CO_all_iso_HITEMP",
        "CO_all_iso_Chubb",
        "CO2",
        "H2O",
        "H2O_HITEMP",
        "H2S",
        "HCN",
        "K",
        "K_lor_cut",
        "K_allard",
        "K_burrows",
        "NH3",
        "Na",
        "Na_lor_cut",
        "Na_allard",
        "Na_burrows",
        "OH",
        "PH3",
        "TiO",
        "TiO_all_Exomol",
        "TiO_all_Plez",
        "VO",
        "VO_Plez",
        "FeH",
        "H2O_main_iso",
        "CH4_main_iso",
    ]


@typechecked
def pt_ret_model(
    temp_3: Optional[np.ndarray],
    delta: float,
    alpha: float,
    tint: float,
    press: np.ndarray,
    metallicity: float,
    c_o_ratio: float,
    conv: bool = True,
) -> Tuple[Optional[np.ndarray], Optional[float], Optional[float]]:
    """
    Pressure-temperature profile for a self-luminous atmosphere (see
    Mollière et al. 2020).

    Parameters
    ----------
    temp_3 : np.ndarray, None
        Array with three temperature points that are added on top of
        the radiative Eddington structure (i.e. above tau = 0.1). The
        temperature nodes are connected with a spline interpolation
        and a prior is used such that t1 < t2 < t3 < t_connect. The
        three temperature points are not used if set to ``None``.
    delta : float
        Proportionality factor in tau = delta * press_cgs**alpha.
    alpha : float
        Power law index in
        :math:`\\tau = \\delta * P_\\mathrm{cgs}**\\alpha`.
        For the tau model: use the proximity to the
        :math:`\\kappa_\\mathrm{rosseland}` photosphere as prior.
    tint : float
        Internal temperature for the Eddington model.
    press : np.ndarray
        Pressure profile (bar).
    metallicity : float
        Metallicity [Fe/H]. Required for the ``nabla_ad``
        interpolation.
    c_o_ratio : float
        Carbon-to-oxygen ratio. Required for the ``nabla_ad``
        interpolation.
    conv : bool
        Enforce a convective adiabat.

    Returns
    -------
    np.ndarray
        Temperature profile (K) for ``press``.
    float
        Pressure (bar) where the optical depth is 1.
    float, None
        Pressure (bar) at the radiative-convective boundary.
    """

    # Convert pressures from bar to cgs units
    press_cgs = press * 1e6

    # Calculate the optical depth
    tau = delta * press_cgs**alpha

    # Calculate the Eddington temperature
    tedd = (3.0 / 4.0 * tint**4.0 * (2.0 / 3.0 + tau)) ** 0.25

    # Import interpol_abundances here because it slows down importing
    # species otherwise. Importing interpol_abundances is only slow
    # the first time, which occurs at the start of the run_multinest
    # method of AtmosphericRetrieval

    if "poor_mans_nonequ_chem" in sys.modules:
        from poor_mans_nonequ_chem.poor_mans_nonequ_chem import interpol_abundances
    else:
        from petitRADTRANS.poor_mans_nonequ_chem.poor_mans_nonequ_chem import (
            interpol_abundances,
        )

    ab = interpol_abundances(
        np.full(tedd.shape[0], c_o_ratio),
        np.full(tedd.shape[0], metallicity),
        tedd,
        press,
    )

    nabla_ad = ab["nabla_ad"]

    # Enforce convective adiabat
    if conv:
        # Calculate the current, radiative temperature gradient
        nab_rad = np.diff(np.log(tedd)) / np.diff(np.log(press_cgs))

        # Extend to array of same length as pressure structure
        nabla_rad = np.ones_like(tedd)
        nabla_rad[0] = nab_rad[0]
        nabla_rad[-1] = nab_rad[-1]
        nabla_rad[1:-1] = (nab_rad[1:] + nab_rad[:-1]) / 2.0

        # Where is the atmosphere convectively unstable?
        conv_index = nabla_rad > nabla_ad

        if np.argwhere(conv_index).size == 0:
            conv_press = None

        else:
            conv_bound = np.amin(np.argwhere(conv_index))
            conv_press = press[conv_bound]

        tfinal = None

        for i in range(10):
            if i == 0:
                t_take = copy.copy(tedd)
            else:
                t_take = copy.copy(tfinal)

            ab = interpol_abundances(
                np.full(t_take.shape[0], c_o_ratio),
                np.full(t_take.shape[0], metallicity),
                t_take,
                press,
            )

            nabla_ad = ab["nabla_ad"]

            # Calculate the average nabla_ad between the layers
            nabla_ad_mean = nabla_ad
            nabla_ad_mean[1:] = (nabla_ad[1:] + nabla_ad[:-1]) / 2.0

            # What are the increments in temperature due to convection
            tnew = nabla_ad_mean[conv_index] * np.mean(np.diff(np.log(press_cgs)))

            # What is the last radiative temperature?
            tstart = np.log(t_take[~conv_index][-1])

            # Integrate and translate to temperature
            # from log(temperature)
            tnew = np.exp(np.cumsum(tnew) + tstart)

            # Add upper radiative and lower covective
            # part into one single array
            tfinal = copy.copy(t_take)
            tfinal[conv_index] = tnew

            if np.max(np.abs(t_take - tfinal) / t_take) < 0.01:
                # print('n_ad', 1./(1.-nabla_ad[conv_index]))
                break

    else:
        tfinal = tedd
        conv_press = None

    # Add the three temperature-point P-T description above tau = 0.1

    @typechecked
    def press_tau(tau: float) -> float:
        """
        Function to return the pressure in cgs units at a given
        optical depth.

        Parameters
        ----------
        tau : float
            Optical depth.

        Returns
        -------
        float
            Pressure (cgs) at optical depth ``tau``.
        """

        return (tau / delta) ** (1.0 / alpha)

    # Where is the uppermost pressure of the
    # Eddington radiative structure?
    p_bot_spline = press_tau(0.1)

    if temp_3 is None:
        tret = tfinal

    else:
        for i_intp in range(2):
            if i_intp == 0:
                # Create the pressure coordinates for the spline
                # support nodes at low pressure
                support_points_low = np.logspace(
                    np.log10(press_cgs[0]), np.log10(p_bot_spline), 4
                )

                # Create the pressure coordinates for the spline
                # support nodes at high pressure, the corresponding
                # temperatures for these nodes will be taken from
                # the radiative-convective solution
                support_points_high = 10.0 ** np.arange(
                    np.log10(p_bot_spline),
                    np.log10(press_cgs[-1]),
                    np.diff(np.log10(support_points_low))[0],
                )

                # Combine into one support node array, don't add
                # the p_bot_spline point twice.
                support_points = np.zeros(
                    len(support_points_low) + len(support_points_high) - 1
                )

                support_points[:4] = support_points_low
                support_points[4:] = support_points_high[1:]

            else:
                # Create the pressure coordinates for the spline
                # support nodes at low pressure
                support_points_low = np.logspace(
                    np.log10(press_cgs[0]), np.log10(p_bot_spline), 7
                )

                # Create the pressure coordinates for the spline
                # support nodes at high pressure, the corresponding
                # temperatures for these nodes will be taken from
                # the radiative-convective solution
                support_points_high = np.logspace(
                    np.log10(p_bot_spline), np.log10(press_cgs[-1]), 7
                )

                # Combine into one support node array, don't add
                # the p_bot_spline point twice.
                support_points = np.zeros(
                    len(support_points_low) + len(support_points_high) - 1
                )
                support_points[:7] = support_points_low
                support_points[7:] = support_points_high[1:]

            # Define the temperature values at the node points
            t_support = np.zeros_like(support_points)

            if i_intp == 0:
                tfintp = interp1d(press_cgs, tfinal)

                # The temperature at p_bot_spline (from the
                # radiative-convective solution)
                t_support[len(support_points_low) - 1] = tfintp(p_bot_spline)

                # if temp_3 is not None:
                # The temperature at pressures below
                # p_bot_spline (free parameters)
                t_support[: len(support_points_low) - 1] = temp_3

                # else:
                #     t_support[:3] = tfintp(support_points_low[:3])

                # The temperature at pressures above p_bot_spline
                # (from the radiative-convective solution)
                t_support[len(support_points_low) :] = tfintp(
                    support_points[len(support_points_low) :]
                )

            else:
                tfintp1 = interp1d(press_cgs, tret)

                t_support[: len(support_points_low) - 1] = tfintp1(
                    support_points[: len(support_points_low) - 1]
                )

                tfintp = interp1d(press_cgs, tfinal)

                # The temperature at p_bot_spline (from
                # the radiative-convective solution)
                t_support[len(support_points_low) - 1] = tfintp(p_bot_spline)

                try:
                    t_support[len(support_points_low) :] = tfintp(
                        support_points[len(support_points_low) :]
                    )

                except ValueError:
                    return None, None, None

            # Make the temperature spline interpolation to be returned
            # to the user tret = spline(np.log10(support_points),
            # t_support, np.log10(press_cgs), order = 3)
            cs = PchipInterpolator(np.log10(support_points), t_support)
            tret = cs(np.log10(press_cgs))

    # Return the temperature, the pressure at tau = 1
    # The temperature at the connection point: tfintp(p_bot_spline)
    # The last two are needed for the priors on the P-T profile.
    return tret, press_tau(1.0) / 1e6, conv_press


@typechecked
def pt_spline_interp(
    knot_press: np.ndarray,
    knot_temp: np.ndarray,
    pressure: np.ndarray,
    pt_smooth: Optional[float] = 0.3,
) -> np.ndarray:
    """
    Function for interpolating the P-T nodes with a PCHIP 1-D monotonic
    cubic interpolation. The interpolated temperature is smoothed with
    a Gaussian kernel of width 0.3 dex in pressure (see Piette &
    Madhusudhan 2020).

    Parameters
    ----------
    knot_press : np.ndarray
        Pressure knots (bar).
    knot_temp : np.ndarray
        Temperature knots (K).
    pressure : np.ndarray
        Pressure points (bar) at which the temperatures is
        interpolated.
    pt_smooth : float, dict, None
        Standard deviation of the Gaussian kernel that is used for
        smoothing the P-T profile, after the temperature nodes
        have been interpolated to a higher pressure resolution.
        The argument should be given as
        :math:`\\log10{P/\\mathrm{bar}}`, with the default value
        set to 0.3 dex. No smoothing is applied if the argument
        is set to ``None``.

    Returns
    -------
    np.ndarray
        Interpolated, smoothed temperature points (K).
    """

    pt_interp = PchipInterpolator(np.log10(knot_press), knot_temp)

    temp_interp = pt_interp(np.log10(pressure))

    log_press = np.log10(pressure)
    log_diff = np.mean(np.diff(log_press))

    if np.std(np.diff(log_press)) / log_diff > 1e-6:
        raise ValueError("Expecting equally spaced pressures in log space.")

    if pt_smooth is not None:
        temp_interp = gaussian_filter(
            temp_interp, sigma=pt_smooth / log_diff, mode="nearest"
        )

    return temp_interp


@typechecked
def create_pt_profile(
    cube,
    cube_index: Dict[str, float],
    pt_profile: str,
    pressure: np.ndarray,
    knot_press: Optional[np.ndarray],
    metallicity: float,
    c_o_ratio: float,
    pt_smooth: Optional[Union[float, Dict[str, float]]] = 0.3,
) -> Tuple[
    Optional[np.ndarray], Optional[np.ndarray], Optional[float], Optional[float]
]:
    """
    Function for creating a pressure-temperature profile.

    Parameters
    ----------
    cube : LP_c_double
        Unit cube.
    cube_index : dict
        Dictionary with the index of each parameter in the ``cube``.
    pt_profile : str
        The parametrization for the pressure-temperature profile
        ('molliere', 'free', 'monotonic', 'eddington', 'gradient').
    pressure : np.ndarray
        Pressure points (bar) at which the temperatures is
        interpolated.
    knot_press : np.ndarray, None
        Pressure knots (bar), which are required when the argument of
        ``pt_profile`` is either 'free' or 'monotonic'.
    metallicity : float
        Metallicity [Fe/H].
    c_o_ratio : float, None
        Carbon-to-oxgen ratio.
    pt_smooth : float, dict, None
        Standard deviation of the Gaussian kernel that is used for
        smoothing the P-T profile, after the temperature nodes
        have been interpolated to a higher pressure resolution.
        The argument should be given as
        :math:`\\log10{P/\\mathrm{bar}}`, with the default value
        set to 0.3 dex. No smoothing is applied if the argument
        is set to ``None``.

    Returns
    -------
    np.ndarray
        Temperatures (K).
    np.ndarray, None
        Temperature at the knots (K). A ``None`` is returned if
        ``pt_profile`` is set to 'molliere' or 'eddington'.
    float, None
        Pressure (bar) where the optical depth is 1.
    float, None
        Pressure (bar) at the radiative-convective boundary.
    """

    knot_temp = None

    if pt_profile == "molliere":
        temp, phot_press, conv_press = pt_ret_model(
            np.array(
                [cube[cube_index["t1"]], cube[cube_index["t2"]], cube[cube_index["t3"]]]
            ),
            10.0 ** cube[cube_index["log_delta"]],
            cube[cube_index["alpha"]],
            cube[cube_index["tint"]],
            pressure,
            metallicity,
            c_o_ratio,
        )

    elif pt_profile == "mod-molliere":
        temp, phot_press, conv_press = pt_ret_model(
            None,
            10.0 ** cube[cube_index["log_delta"]],
            cube[cube_index["alpha"]],
            cube[cube_index["tint"]],
            pressure,
            metallicity,
            c_o_ratio,
        )

    elif pt_profile == "gradient":
        num_layer = 6  # could make a variable in the future
        layer_pt_slopes = np.ones(num_layer) * np.nan

        for index in range(num_layer):
            layer_pt_slopes[index] = cube[cube_index[f"PTslope_{num_layer - index}"]]

        try:
            from petitRADTRANS.physics import dTdP_temperature_profile

        except ImportError as import_error:
            raise ImportError(
                "Can't import the dTdP profile function from "
                "petitRADTRANS, check that the version of pRT "
                "includes this function in petitRADTRANS.physics",
                import_error,
            )

        temp = dTdP_temperature_profile(
            pressure,
            num_layer,  # could change in the future
            layer_pt_slopes,
            cube[cube_index["T_bottom"]],
        )

        phot_press = None
        conv_press = None

    elif pt_profile in ["free", "monotonic"]:
        knot_temp = []
        for i in range(knot_press.shape[0]):
            knot_temp.append(cube[cube_index[f"t{i}"]])

        knot_temp = np.asarray(knot_temp)

        nan_count = np.sum(np.isnan(knot_temp))

        if nan_count > 0:
            # TODO Not clear why this can happen
            warnings.warn(f"Found {nan_count} NaN values in sampled temperature nodes.")

            return None, None, None, None

        temp = pt_spline_interp(knot_press, knot_temp, pressure, pt_smooth)

        phot_press = None
        conv_press = None

    elif pt_profile == "eddington":
        # Eddington approximation
        # delta = kappa_ir/gravity
        tau = pressure * 1e6 * 10.0 ** cube[cube_index["log_delta"]]
        temp = (0.75 * cube[cube_index["tint"]] ** 4.0 * (2.0 / 3.0 + tau)) ** 0.25

        phot_press = None
        conv_press = None

    return temp, knot_temp, phot_press, conv_press


@typechecked
def make_half_pressure_better(
    p_base: Dict[str, float], pressure: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for reducing the number of pressure layers from 1440 to
    ~100 (depending on the number of cloud species) with a refinement
    around the cloud decks.

    Parameters
    ----------
    p_base : dict
        Dictionary with the base of the cloud deck for all cloud
        species. The keys in the dictionary are included for example
        as MgSiO3(c).
    pressure : np.ndarray
        Pressures (bar) at high resolution (1440 points).

    Returns
    -------
    np.ndarray
        Pressures (bar) at lower resolution (60 points) but with a
        refinement around the position of the cloud decks.
    np.ndarray, None
        The indices of the pressures that have been selected from
        the input array ``pressure``.
    """

    press_plus_index = np.zeros(len(pressure) * 2).reshape(len(pressure), 2)
    press_plus_index[:, 0] = pressure
    press_plus_index[:, 1] = range(len(pressure))

    press_small = press_plus_index[::24, :]
    press_plus_index = press_plus_index[::2, :]

    indexes_small = press_small[:, 0] > 0.0
    indexes = press_plus_index[:, 0] > 0.0

    for key, P_cloud in p_base.items():
        indexes_small = indexes_small & (
            (np.log10(press_small[:, 0] / P_cloud) > 0.05)
            | (np.log10(press_small[:, 0] / P_cloud) < -0.3)
        )

        indexes = indexes & (
            (np.log10(press_plus_index[:, 0] / P_cloud) > 0.05)
            | (np.log10(press_plus_index[:, 0] / P_cloud) < -0.3)
        )

    press_cut = press_plus_index[~indexes, :]
    press_small_cut = press_small[indexes_small, :]

    press_out = np.zeros((len(press_cut) + len(press_small_cut)) * 2).reshape(
        (len(press_cut) + len(press_small_cut)), 2
    )

    press_out[: len(press_small_cut), :] = press_small_cut
    press_out[len(press_small_cut) :, :] = press_cut

    press_out = np.sort(press_out, axis=0)

    return press_out[:, 0], press_out[:, 1].astype("int")


@typechecked
def create_abund_dict(
    abund_in: Dict[str, np.ndarray],
    line_species: List[str],
    chemistry: str,
    pressure_grid: str = "smaller",
    indices: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """
    Function to update the names in the abundance dictionary.

    Parameters
    ----------
    abund_in : dict
        Dictionary with the mass fractions.
    line_species : list
        List with the line species.
    chemistry : str
        Chemistry type ('equilibrium' or 'free').
    pressure_grid : str
        The type of pressure grid that is used for the radiative
        transfer. Either 'standard', to use 180 layers both for the
        atmospheric structure (e.g. when interpolating the abundances)
        and 180 layers with the radiative transfer, or 'smaller' to
        use 60 (instead of 180) with the radiative transfer, or 'clouds'
        to start with 1440 layers but resample to ~100 layers (depending
        on the number of cloud species) with a refinement around the
        cloud decks. For cloudless atmospheres it is recommended to use
        'smaller', which runs faster than 'standard' and provides
        sufficient accuracy. For cloudy atmosphere, one can test with
        'smaller' but it is recommended to use 'clouds' for improved
        accuracy fluxes.
    indices : np.ndarray, None
        Pressure indices from the adaptive refinement in a cloudy
        atmosphere. Only required with ``pressure_grid='clouds'``.
        Otherwise, the argument can be set to ``None``.

    Returns
    -------
    dict
        Dictionary with the updated names of the abundances.
    """

    # Create a dictionary with the updated abundance names

    abund_out = {}

    if indices is not None:
        for item in line_species:
            if chemistry == "equilibrium":
                item_replace = item.replace("_R_10", "")
                item_replace = item_replace.replace("_R_30", "")
                item_replace = item_replace.replace("_all_iso_HITEMP", "")
                item_replace = item_replace.replace("_all_iso_Chubb", "")
                item_replace = item_replace.replace("_all_iso", "")
                item_replace = item_replace.replace("_HITEMP", "")
                item_replace = item_replace.replace("_main_iso", "")
                item_replace = item_replace.replace("_lor_cut", "")
                item_replace = item_replace.replace("_allard", "")
                item_replace = item_replace.replace("_burrows", "")
                item_replace = item_replace.replace("_all_Plez", "")
                item_replace = item_replace.replace("_all_Exomol", "")
                item_replace = item_replace.replace("_Plez", "")

                abund_out[item] = abund_in[item_replace][indices]

            elif chemistry == "free":
                abund_out[item] = abund_in[item][indices]

        if "Fe(c)" in abund_in:
            abund_out["Fe(c)"] = abund_in["Fe(c)"][indices]

        if "MgSiO3(c)" in abund_in:
            abund_out["MgSiO3(c)"] = abund_in["MgSiO3(c)"][indices]

        if "Mg2SiO4(c)" in abund_in:
            abund_out["Mg2SiO4(c)"] = abund_in["Mg2SiO4(c)"][indices]

        if "Al2O3(c)" in abund_in:
            abund_out["Al2O3(c)"] = abund_in["Al2O3(c)"][indices]

        if "Na2S(c)" in abund_in:
            abund_out["Na2S(c)"] = abund_in["Na2S(c)"][indices]

        if "KCL(c)" in abund_in:
            abund_out["KCL(c)"] = abund_in["KCL(c)"][indices]

        abund_out["H2"] = abund_in["H2"][indices]
        abund_out["He"] = abund_in["He"][indices]

    elif pressure_grid == "smaller":
        for item in line_species:
            if chemistry == "equilibrium":
                item_replace = item.replace("_R_10", "")
                item_replace = item_replace.replace("_R_30", "")
                item_replace = item_replace.replace("_all_iso_HITEMP", "")
                item_replace = item_replace.replace("_all_iso_Chubb", "")
                item_replace = item_replace.replace("_all_iso", "")
                item_replace = item_replace.replace("_HITEMP", "")
                item_replace = item_replace.replace("_main_iso", "")
                item_replace = item_replace.replace("_lor_cut", "")
                item_replace = item_replace.replace("_allard", "")
                item_replace = item_replace.replace("_burrows", "")
                item_replace = item_replace.replace("_all_Plez", "")
                item_replace = item_replace.replace("_all_Exomol", "")
                item_replace = item_replace.replace("_Plez", "")

                abund_out[item] = abund_in[item_replace][::3]

            elif chemistry == "free":
                abund_out[item] = abund_in[item][::3]

        if "Fe(c)" in abund_in:
            abund_out["Fe(c)"] = abund_in["Fe(c)"][::3]

        if "MgSiO3(c)" in abund_in:
            abund_out["MgSiO3(c)"] = abund_in["MgSiO3(c)"][::3]

        if "Mg2SiO4(c)" in abund_in:
            abund_out["Mg2SiO4(c)"] = abund_in["Mg2SiO4(c)"][::3]

        if "Al2O3(c)" in abund_in:
            abund_out["Al2O3(c)"] = abund_in["Al2O3(c)"][::3]

        if "Na2S(c)" in abund_in:
            abund_out["Na2S(c)"] = abund_in["Na2S(c)"][::3]

        if "KCL(c)" in abund_in:
            abund_out["KCL(c)"] = abund_in["KCL(c)"][::3]

        abund_out["H2"] = abund_in["H2"][::3]
        abund_out["He"] = abund_in["He"][::3]

    else:
        for item in line_species:
            if chemistry == "equilibrium":
                item_replace = item.replace("_R_10", "")
                item_replace = item_replace.replace("_R_30", "")
                item_replace = item_replace.replace("_all_iso_HITEMP", "")
                item_replace = item_replace.replace("_all_iso_Chubb", "")
                item_replace = item_replace.replace("_all_iso", "")
                item_replace = item_replace.replace("_HITEMP", "")
                item_replace = item_replace.replace("_main_iso", "")
                item_replace = item_replace.replace("_lor_cut", "")
                item_replace = item_replace.replace("_allard", "")
                item_replace = item_replace.replace("_burrows", "")
                item_replace = item_replace.replace("_all_Plez", "")
                item_replace = item_replace.replace("_all_Exomol", "")
                item_replace = item_replace.replace("_Plez", "")

                abund_out[item] = abund_in[item_replace]

            elif chemistry == "free":
                abund_out[item] = abund_in[item]

        if "Fe(c)" in abund_in:
            abund_out["Fe(c)"] = abund_in["Fe(c)"]

        if "MgSiO3(c)" in abund_in:
            abund_out["MgSiO3(c)"] = abund_in["MgSiO3(c)"]

        if "Mg2SiO4(c)" in abund_in:
            abund_out["Mg2SiO4(c)"] = abund_in["Mg2SiO4(c)"]

        if "Al2O3(c)" in abund_in:
            abund_out["Al2O3(c)"] = abund_in["Al2O3(c)"]

        if "Na2S(c)" in abund_in:
            abund_out["Na2S(c)"] = abund_in["Na2S(c)"]

        if "KCL(c)" in abund_in:
            abund_out["KCL(c)"] = abund_in["KCL(c)"]

        abund_out["H2"] = abund_in["H2"]
        abund_out["He"] = abund_in["He"]

    # Correction for the nuclear spin degeneracy that was not included
    # in the partition function. See Charnay et al. (2018)

    if "FeH" in abund_out:
        abund_out["FeH"] = abund_out["FeH"] / 2.0

    return abund_out


@typechecked
def calc_spectrum_clear(
    rt_object,
    pressure: np.ndarray,
    temperature: np.ndarray,
    log_g: float,
    c_o_ratio: Optional[float],
    metallicity: Optional[float],
    p_quench: Optional[float],
    log_x_abund: Optional[dict],
    chemistry: str,
    knot_press_abund: Optional[np.ndarray],
    abund_smooth: Optional[float],
    pressure_grid: str = "smaller",
    contribution: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Function to simulate an emission spectrum of a clear atmosphere.
    The function supports both equilibrium chemistry
    (``chemistry='equilibrium'``) and free abundances
    (``chemistry='free'``).

    rt_object : petitRADTRANS.radtrans.Radtrans
        Instance of ``Radtrans``.
    pressure : np.ndarray
        Array with the pressure points (bar).
    temperature : np.ndarray
        Array with the temperature points (K) corresponding to
        ``pressure``.
    log_g : float
        Log10 of the surface gravity (cm s-2).
    c_o_ratio : float, None
        Carbon-to-oxygen ratio.
    metallicity : float, None
        Metallicity.
    p_quench : float, None
        Quenching pressure (bar).
    log_x_abund : dict, None
        Dictionary with the log10 of the abundances. Only required when
        ``chemistry='free'``.
    chemistry : str
        Chemistry type ('equilibrium' or 'free').
    knot_press_abund : np.ndarray, None
        Pressure nodes at which the abundances are sampled. Only
        required when ``chemistry='free'``.
    abund_smooth : float, None
        Standard deviation of the Gaussian kernel that is used for
        smoothing the abundance profiles, after the abundance nodes
        have been interpolated to a higher pressure resolution.
        Only required with ```chemistry='free'``` and
        ``knot_press_abund`` is not set to ``None``. The argument
        should be given as :math:`\\log10{P/\\mathrm{bar}}`. No
        smoothing is applied if the argument if set to 0 or ``None``.
    pressure_grid : str
        The type of pressure grid that is used for the radiative
        transfer. Either 'standard', to use 180 layers both for the
        atmospheric structure (e.g. when interpolating the abundances)
        and 180 layers with the radiative transfer, or 'smaller' to use
        60 (instead of 180) with the radiative transfer, or 'clouds' to
        start with 1440 layers but resample to ~100 layers (depending
        on the number of cloud species) with a refinement around the
        cloud decks. For cloudless atmospheres it is recommended to use
        'smaller', which runs faster than 'standard' and provides
        sufficient accuracy. For cloudy atmosphere, one can test with
        'smaller' but it is recommended to use 'clouds' for improved
        accuracy fluxes.
    contribution : bool
        Calculate the emission contribution.

    Returns
    -------
    np.ndarray
        Wavelength (um).
    np.ndarray
        Flux (W m-2 um-1).
    np.ndarray, None
        Emission contribution.
    """

    if knot_press_abund is None:
        abund_nodes = None
    else:
        abund_nodes = knot_press_abund.size

    if chemistry == "equilibrium":
        # Import interpol_abundances here because it slows down
        # importing species otherwise. Importing interpol_abundances
        # is only slow the first time, which occurs at the start of
        # the run_multinest method of AtmosphericRetrieval
        if "poor_mans_nonequ_chem" in sys.modules:
            from poor_mans_nonequ_chem.poor_mans_nonequ_chem import interpol_abundances
        else:
            from petitRADTRANS.poor_mans_nonequ_chem.poor_mans_nonequ_chem import (
                interpol_abundances,
            )

        # Chemical equilibrium
        abund_in = interpol_abundances(
            np.full(pressure.shape, c_o_ratio),
            np.full(pressure.shape, metallicity),
            temperature,
            pressure,
            Pquench_carbon=p_quench,
        )
        # Extract the mean molecular weight
        mmw = abund_in["MMW"]

    elif chemistry == "free":
        # Free abundances

        # Create a dictionary with all mass fractions
        abund_in = mass_fractions(log_x_abund, rt_object.line_species, abund_nodes)

        # Create list of all species
        abund_species = rt_object.line_species.copy()
        abund_species.append("H2")
        abund_species.append("He")

        if abund_nodes is None:
            # Mean molecular weight

            mmw_float = mean_molecular_weight(abund_in)
            mmw = np.full(pressure.size, mmw_float)

            # Create arrays of constant abundances with pressure

            for abund_item in abund_species:
                abund_in[abund_item] = np.full(pressure.size, abund_in[abund_item])

        else:
            # Create arrays of pressure-dependent abundances

            for abund_item in abund_species:
                knot_abund = np.zeros(knot_press_abund.shape)

                for node_idx in range(abund_nodes):
                    knot_abund[node_idx] = abund_in[f"{abund_item}_{node_idx}"]
                    del abund_in[f"{abund_item}_{node_idx}"]

                nan_count = np.sum(np.isnan(knot_abund))

                if nan_count > 0:
                    warnings.warn(
                        f"Found {nan_count} NaN values in sampled abundance nodes."
                    )

                    return None, None, None

                abund_in[abund_item] = pt_spline_interp(
                    knot_press_abund, knot_abund, pressure, pt_smooth=abund_smooth
                )

            # Mean molecular weight

            mmw = np.zeros(pressure.size)

            for press_idx in range(pressure.size):
                abund_dict = {}
                for abund_item, abund_value in abund_in.items():
                    abund_dict[abund_item] = abund_value[press_idx]

                mmw[press_idx] = mean_molecular_weight(abund_dict)

    # Extract every three levels when pressure_grid is set to 'smaller'

    if pressure_grid == "smaller":
        temperature = temperature[::3]
        pressure = pressure[::3]
        mmw = mmw[::3]

    abundances = create_abund_dict(
        abund_in,
        rt_object.line_species,
        chemistry,
        pressure_grid=pressure_grid,
        indices=None,
    )

    # calculate the emission spectrum
    rt_object.calc_flux(
        temperature, abundances, 10.0**log_g, mmw, contribution=contribution
    )

    # convert frequency (Hz) to wavelength (cm)
    wavel = constants.LIGHT * 1e2 / rt_object.freq

    # optionally return the emission contribution
    if contribution:
        contr_em = rt_object.contr_em
    else:
        contr_em = None

    # return wavelength (micron), flux (W m-2 um-1),
    # and emission contribution
    return (
        1e4 * wavel,
        1e-7 * rt_object.flux * constants.LIGHT * 1e2 / wavel**2.0,
        contr_em,
    )


@typechecked
def calc_spectrum_clouds(
    rt_object,
    pressure: np.ndarray,
    temperature: np.ndarray,
    c_o_ratio: float,
    metallicity: float,
    p_quench: Optional[float],
    log_x_abund: Optional[dict],
    log_x_base: Optional[dict],
    cloud_dict: Dict[str, Optional[float]],
    log_g: float,
    chemistry: str,
    knot_press_abund: Optional[np.ndarray],
    abund_smooth: Optional[float],
    pressure_grid: str = "smaller",
    plotting: bool = False,
    contribution: bool = False,
    tau_cloud: Optional[float] = None,
    cloud_wavel: Optional[Tuple[float, float]] = None,
) -> Tuple[
    Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], np.ndarray
]:
    """
    Function to simulate an emission spectrum of a cloudy atmosphere.

    Parameters
    ----------
    rt_object : petitRADTRANS.radtrans.Radtrans
        Instance of ``Radtrans``.
    pressure : np.ndarray
        Array with the pressure points (bar).
    temperature : np.ndarray
        Array with the temperature points (K) corresponding to
        ``pressure``.
    c_o_ratio : float
        Carbon-to-oxygen ratio.
    metallicity : float
        Metallicity.
    p_quench : float, None
        Quenching pressure (bar).
    log_x_abund : dict, None
        Dictionary with the log10 of the abundances. Only required
        when ``chemistry='free'``.
    log_x_base : dict, None
        Dictionary with the log10 of the mass fractions at the cloud
        base. Only required when the ``cloud_dict`` contains ``fsed``,
        ``log_kzz``, and ``sigma_lnorm``.
    cloud_dict : dict
        Dictionary with the cloud parameters.
    log_g : float
        Log10 of the surface gravity (cm s-2).
    chemistry : str
        Chemistry type ('equilibrium' or 'free').
    knot_press_abund : np.ndarray, None
        Pressure nodes at which the abundances are sampled. Only
        required when ``chemistry='free'``.
    abund_smooth : float, None
        Standard deviation of the Gaussian kernel that is used for
        smoothing the abundance profiles, after the abundance nodes
        have been interpolated to a higher pressure resolution.
        Only required with ```chemistry='free'``` and
        ``knot_press_abund`` is not set to ``None``. The argument
        should be given as :math:`\\log10{P/\\mathrm{bar}}`. No
        smoothing is applied if the argument if set to 0 or ``None``.
    pressure_grid : str
        The type of pressure grid that is used for the radiative
        transfer. Either 'standard', to use 180 layers both for the
        atmospheric structure (e.g. when interpolating the abundances)
        and 180 layers with the radiative transfer, or 'smaller' to
        use 60 (instead of 180) with the radiative transfer, or
        'clouds' to start with 1440 layers but resample to ~100 layers
        (depending on the number of cloud species) with a refinement
        around the cloud decks. For cloudless atmospheres it is
        recommended to use 'smaller', which runs faster than 'standard'
        and provides sufficient accuracy. For cloudy atmosphere, one
        can test with 'smaller' but it is recommended to use 'clouds'
        for improved accuracy fluxes.
    plotting : bool
        Create plots.
    contribution : bool
        Calculate the emission contribution.
    tau_cloud : float, None
        Total cloud optical that will be used for scaling the cloud
        mass fractions. The mass fractions will not be scaled if the
        parameter is set to ``None``.
    cloud_wavel : tuple(float, float), None
        Tuple with the wavelength range (um) that is used for
        calculating the median optical depth of the clouds at the
        gas-only photosphere and then scaling the cloud optical
        depth to the value of ``log_tau_cloud``. The range of
        ``cloud_wavel`` should be encompassed by the range of
        ``wavel_range``.  The full wavelength range (i.e.
        ``wavel_range``) is used if the argument is set to ``None``.

    Returns
    -------
    np.ndarray, None
        Wavelength (um).
    np.ndarray, None
        Flux (W m-2 um-1).
    np.ndarray, None
        Emission contribution.
    np.ndarray
        Array with mean molecular weight.
    """

    if knot_press_abund is None:
        abund_nodes = None
    else:
        abund_nodes = knot_press_abund.size

    if chemistry == "equilibrium":
        # Import interpol_abundances here because it slows down
        # importing species otherwise. Importing interpol_abundances
        # is only slow the first time, which occurs at the start
        # of the run_multinest method of AtmosphericRetrieval

        if "poor_mans_nonequ_chem" in sys.modules:
            from poor_mans_nonequ_chem.poor_mans_nonequ_chem import interpol_abundances
        else:
            from petitRADTRANS.poor_mans_nonequ_chem.poor_mans_nonequ_chem import (
                interpol_abundances,
            )

        # Interpolate the abundances, following chemical equilibrium
        abund_in = interpol_abundances(
            np.full(pressure.shape, c_o_ratio),
            np.full(pressure.shape, metallicity),
            temperature,
            pressure,
            Pquench_carbon=p_quench,
        )

        # Extract the mean molecular weight
        mmw = abund_in["MMW"]

    elif chemistry == "free":
        # Free abundances

        # Create a dictionary with all mass fractions
        abund_in = mass_fractions(log_x_abund, rt_object.line_species, abund_nodes)

        # Create list of all species
        abund_species = rt_object.line_species.copy()
        abund_species.append("H2")
        abund_species.append("He")

        if abund_nodes is None:
            # Mean molecular weight

            mmw_float = mean_molecular_weight(abund_in)
            mmw = np.full(pressure.size, mmw_float)

            # Create arrays of constant abundances with pressure

            for abund_item in abund_species:
                abund_in[abund_item] = np.full(pressure.size, abund_in[abund_item])

        else:
            # Create arrays of pressure-dependent abundances

            for abund_item in abund_species:
                knot_abund = np.zeros(knot_press_abund.shape)

                for node_idx in range(abund_nodes):
                    knot_abund[node_idx] = abund_in[f"{abund_item}_{node_idx}"]
                    del abund_in[f"{abund_item}_{node_idx}"]

                nan_count = np.sum(np.isnan(knot_abund))

                if nan_count > 0:
                    warnings.warn(
                        f"Found {nan_count} NaN values in sampled abundance nodes."
                    )

                    return None, None, None, np.zeros(1)

                abund_in[abund_item] = pt_spline_interp(
                    knot_press_abund, knot_abund, pressure, pt_smooth=abund_smooth
                )

            # Mean molecular weight

            mmw = np.zeros(pressure.size)

            for press_idx in range(pressure.size):
                abund_dict = {}
                for abund_item, abund_value in abund_in.items():
                    abund_dict[abund_item] = abund_value[press_idx]

                mmw[press_idx] = mean_molecular_weight(abund_dict)

    if log_x_base is not None:
        p_base = {}

        for cloud_item in log_x_base:
            if f"log_p_base_{cloud_item}(c)" in cloud_dict:
                p_base_item = 10.0 ** cloud_dict[f"log_p_base_{cloud_item}(c)"]

            else:
                p_base_item = find_cloud_deck(
                    cloud_item,
                    pressure,
                    temperature,
                    metallicity,
                    c_o_ratio,
                    mmw=np.mean(mmw),
                    plotting=plotting,
                )

            p_base[f"{cloud_item}(c)"] = p_base_item

            abund_in[f"{cloud_item}(c)"] = np.zeros_like(temperature)

            if "fsed" in cloud_dict:
                f_sed = cloud_dict["fsed"]
            else:
                f_sed = cloud_dict[f"fsed_{cloud_item}(c)"]

            abund_in[f"{cloud_item}(c)"][pressure < p_base_item] = (
                10.0 ** log_x_base[cloud_item]
                * (pressure[pressure <= p_base_item] / p_base_item) ** f_sed
            )

    # Adaptive pressure refinement around the cloud base

    if pressure_grid == "clouds":
        _, indices = make_half_pressure_better(p_base, pressure)
    else:
        indices = None

    # Update the abundance dictionary

    abundances = create_abund_dict(
        abund_in,
        rt_object.line_species,
        chemistry,
        pressure_grid=pressure_grid,
        indices=indices,
    )

    # Create dictionary with sedimentation parameters
    # Use the same value for all cloud species

    fseds = {}
    for cloud_item in rt_object.cloud_species:
        # The cloud_item has the form of e.g. MgSiO3(c). For
        # parametrized cloud opacities, then the number of
        # cloud_species is zero so the fseds dictionary is empty
        if "fsed" in cloud_dict:
            fseds[cloud_item] = cloud_dict["fsed"]
        else:
            fseds[cloud_item] = cloud_dict[f"fsed_{cloud_item}"]

    # Create an array with a constant eddy diffusion coefficient (cm2 s-1)

    if "log_kzz" in cloud_dict:
        Kzz_use = np.full(pressure.shape, 10.0 ** cloud_dict["log_kzz"])
    else:
        Kzz_use = None

    # Adjust number of atmospheric levels

    if pressure_grid == "smaller":
        temperature = temperature[::3]
        pressure = pressure[::3]
        mmw = mmw[::3]

        if "log_kzz" in cloud_dict:
            Kzz_use = Kzz_use[::3]

    elif pressure_grid == "clouds":
        temperature = temperature[indices]
        pressure = pressure[indices]
        mmw = mmw[indices]

        if "log_kzz" in cloud_dict:
            Kzz_use = Kzz_use[indices]

    # Optionally plot the cloud properties

    if (
        plotting
        and Kzz_use is not None
        and (
            rt_object.wlen_bords_micron[0] != 0.5
            and rt_object.wlen_bords_micron[1] != 30.0
        )
    ):
        if "CO_all_iso" in abundances:
            plt.plot(abundances["CO_all_iso"], pressure, label="CO")
        if "CO_all_iso_HITEMP" in abundances:
            plt.plot(abundances["CO_all_iso_HITEMP"], pressure, label="CO")
        if "CO_all_iso_Chubb" in abundances:
            plt.plot(abundances["CO_all_iso_Chubb"], pressure, label="CO")
        if "CH4" in abundances:
            plt.plot(abundances["CH4"], pressure, label="CH4")
        if "H2O" in abundances:
            plt.plot(abundances["H2O"], pressure, label="H2O")
        if "H2O_HITEMP" in abundances:
            plt.plot(abundances["H2O_HITEMP"], pressure, label="H2O")
        if "Na" in abundances:
            plt.plot(abundances["Na"], pressure, label="Na")
        if "K" in abundances:
            plt.plot(abundances["K"], pressure, label="K")
        if "Na_allard" in abundances:
            plt.plot(abundances["Na_allard"], pressure, label="Na")
        if "K_allard" in abundances:
            plt.plot(abundances["K_allard"], pressure, label="K")
        if "Na_burrows" in abundances:
            plt.plot(abundances["Na_burrows"], pressure, label="Na")
        if "K_burrows" in abundances:
            plt.plot(abundances["K_burrows"], pressure, label="K")
        plt.xlim(1e-10, 1.0)
        plt.ylim(pressure[-1], pressure[0])
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("Mass fraction")
        plt.ylabel("Pressure (bar)")
        if p_quench is not None:
            plt.axhline(p_quench, ls="--", color="black")
        plt.legend(loc="best")
        plt.savefig("abundances.png", bbox_inches="tight")
        plt.clf()

        plt.plot(temperature, pressure, "o", ls="none", ms=2.0)

        for item in log_x_base:
            plt.axhline(
                p_base[f"{item}(c)"], label=f"Cloud deck {item}", ls="--", color="black"
            )

        plt.yscale("log")
        plt.ylim(1e3, 1e-6)
        plt.xlim(0.0, 4000.0)
        plt.savefig("pt_cloud_deck.png", bbox_inches="tight")
        plt.clf()

        for item in log_x_base:
            plt.plot(abundances[f"{item}(c)"], pressure)
            plt.axhline(p_base[f"{item}(c)"])
            plt.yscale("log")
            if np.count_nonzero(abundances[f"{item}(c)"]) > 0:
                plt.xscale("log")
            plt.ylim(1e3, 1e-6)
            plt.xlim(1e-10, 1.0)
            log_x_base_item = log_x_base[item]
            if "fsed" in cloud_dict:
                fsed = cloud_dict["fsed"]
            else:
                fsed = cloud_dict[f"fsed_{item}(c)"]
            log_kzz = cloud_dict["log_kzz"]
            plt.title(
                f"fsed = {fsed:.2f}, log(Kzz) = {log_kzz:.2f}, "
                + f"X_b = {log_x_base_item:.2f}"
            )
            plt.savefig(f"{item.lower()}_clouds.png", bbox_inches="tight")
            plt.clf()

    # Turn clouds off
    # abundances['MgSiO3(c)'] = np.zeros_like(pressure)
    # abundances['Fe(c)'] = np.zeros_like(pressure)

    # Reinitiate the pressure layers after make_half_pressure_better
    if pressure_grid == "clouds":
        rt_object.setup_opa_structure(pressure)

    # Width of cloud particle distribution
    if "sigma_lnorm" in cloud_dict:
        sigma_lnorm = cloud_dict["sigma_lnorm"]
    else:
        sigma_lnorm = None

    if "log_kappa_0" in cloud_dict:
        # Cloud model 2

        @typechecked
        def kappa_abs(wavel_micron: np.ndarray, press_bar: np.ndarray) -> np.ndarray:
            p_base = 10.0 ** cloud_dict["log_p_base"]  # (bar)
            kappa_0 = 10.0 ** cloud_dict["log_kappa_0"]  # (cm2 g-1)

            # Opacity at 1 um (cm2 g-1) as function of pressure (bar)
            # See Eq. 5 in Mollière et al. 2020
            kappa_p = kappa_0 * (press_bar / p_base) ** cloud_dict["fsed"]

            # Opacity (cm2 g-1) as function of wavelength (um)
            # See Eq. 4 in Mollière et al. 2020
            kappa_grid, wavel_grid = np.meshgrid(kappa_p, wavel_micron, sparse=True)
            kappa_tot = kappa_grid * wavel_grid ** cloud_dict["opa_index"]
            kappa_tot[:, press_bar > p_base] = 0.0

            # if (
            #     cloud_dict["opa_knee"] > wavel_micron[0]
            #     and cloud_dict["opa_knee"] < wavel_micron[-1]
            # ):
            #     indices = np.where(wavel_micron > cloud_dict["opa_knee"])[0]
            #     for i in range(press_bar.size):
            #         kappa_tot[indices, i] = (
            #             kappa_tot[indices[0], i]
            #             * (wavel_micron[indices] / wavel_micron[indices[0]]) ** -4.0
            #         )

            return (1.0 - cloud_dict["albedo"]) * kappa_tot

        @typechecked
        def kappa_scat(wavel_micron: np.ndarray, press_bar: np.ndarray):
            p_base = 10.0 ** cloud_dict["log_p_base"]  # (bar)
            kappa_0 = 10.0 ** cloud_dict["log_kappa_0"]  # (cm2 g-1)

            # Opacity at 1 um (cm2 g-1) as function of pressure (bar)
            # See Eq. 5 in Mollière et al. 2020
            kappa_p = kappa_0 * (press_bar / p_base) ** cloud_dict["fsed"]

            # Opacity (cm2 g-1) as function of wavelength (um)
            # See Eq. 4 in Mollière et al. 2020
            kappa_grid, wavel_grid = np.meshgrid(kappa_p, wavel_micron, sparse=True)
            kappa_tot = kappa_grid * wavel_grid ** cloud_dict["opa_index"]
            kappa_tot[:, press_bar > p_base] = 0.0

            # if (
            #     cloud_dict["opa_knee"] > wavel_micron[0]
            #     and cloud_dict["opa_knee"] < wavel_micron[-1]
            # ):
            #     indices = np.where(wavel_micron > cloud_dict["opa_knee"])[0]
            #     for i in range(press_bar.size):
            #         kappa_tot[indices, i] = (
            #             kappa_tot[indices[0], i]
            #             * (wavel_micron[indices] / wavel_micron[indices[0]]) ** -4.0
            #         )

            return cloud_dict["albedo"] * kappa_tot

    elif "log_kappa_abs" in cloud_dict:
        # Powerlaw absorption and scattering opacities

        @typechecked
        def kappa_abs(wavel_micron: np.ndarray, press_bar: np.ndarray) -> np.ndarray:
            p_base = 10.0 ** cloud_dict["log_p_base"]  # (bar)
            kappa_0 = 10.0 ** cloud_dict["log_kappa_abs"]  # (cm2 g-1)

            # Opacity at 1 um (cm2 g-1) as function of pressure (bar)
            kappa_p = kappa_0 * (press_bar / p_base) ** cloud_dict["fsed"]

            # Opacity (cm2 g-1) as function of wavelength (um)
            kappa_grid, wavel_grid = np.meshgrid(kappa_p, wavel_micron, sparse=True)
            kappa_abs = kappa_grid * wavel_grid ** cloud_dict["opa_abs_index"]
            kappa_abs[:, press_bar > p_base] = 0.0

            return kappa_abs

        if "log_kappa_sca" in cloud_dict:

            @typechecked
            def kappa_scat(wavel_micron: np.ndarray, press_bar: np.ndarray):
                p_base = 10.0 ** cloud_dict["log_p_base"]  # (bar)
                kappa_0 = 10.0 ** cloud_dict["log_kappa_sca"]  # (cm2 g-1)

                # Opacity at 1 um (cm2 g-1) as function of pressure (bar)
                kappa_p = kappa_0 * (press_bar / p_base) ** cloud_dict["fsed"]

                # Opacity (cm2 g-1) as function of wavelength (um)
                kappa_grid, wavel_grid = np.meshgrid(kappa_p, wavel_micron, sparse=True)
                kappa_sca = kappa_grid * wavel_grid ** cloud_dict["opa_sca_index"]
                kappa_sca[:, press_bar > p_base] = 0.0

                if (
                    cloud_dict["lambda_ray"] > wavel_micron[0]
                    and cloud_dict["lambda_ray"] < wavel_micron[-1]
                ):
                    indices = np.where(wavel_micron > cloud_dict["lambda_ray"])[0]
                    for i in range(press_bar.size):
                        kappa_sca[indices, i] = (
                            kappa_sca[indices[0], i]
                            * (wavel_micron[indices] / wavel_micron[indices[0]]) ** -4.0
                        )

                return kappa_sca

        else:
            kappa_scat = None

    elif "log_kappa_gray" in cloud_dict:
        # Gray clouds with cloud top

        @typechecked
        def kappa_abs(wavel_micron: np.ndarray, press_bar: np.ndarray) -> np.ndarray:
            p_top = 10.0 ** cloud_dict["log_cloud_top"]  # (bar)
            kappa_gray = 10.0 ** cloud_dict["log_kappa_gray"]  # (cm2 g-1)

            opa_abs = np.full((wavel_micron.size, press_bar.size), kappa_gray)
            opa_abs[:, press_bar < p_top] = 0.0

            return opa_abs

        # Add optional scattering opacity

        if "albedo" in cloud_dict:

            @typechecked
            def kappa_scat(
                wavel_micron: np.ndarray, press_bar: np.ndarray
            ) -> np.ndarray:
                # Absorption opacity (cm2 g-1)
                opa_abs = kappa_abs(wavel_micron, press_bar)

                # Scattering opacity (cm2 g-1)
                opa_scat = cloud_dict["albedo"] * opa_abs / (1.0 - cloud_dict["albedo"])

                return opa_scat

        else:
            kappa_scat = None

    else:
        kappa_abs = None
        kappa_scat = None

    # Calculate the emission spectrum

    rt_object.calc_flux(
        temperature,
        abundances,
        10.0**log_g,
        mmw,
        sigma_lnorm=sigma_lnorm,
        Kzz=Kzz_use,
        fsed=fseds,
        radius=None,
        contribution=contribution,
        gray_opacity=None,
        Pcloud=None,
        kappa_zero=None,
        gamma_scat=None,
        add_cloud_scat_as_abs=False,
        hack_cloud_photospheric_tau=tau_cloud,
        give_absorption_opacity=kappa_abs,
        give_scattering_opacity=kappa_scat,
        cloud_wlen=cloud_wavel,
    )

    # if (
    #     hasattr(rt_object, "scaling_physicality")
    #     and rt_object.scaling_physicality > 1.0
    # ):
    #     # cloud_scaling_factor > 2 * (fsed + 1)
    #     # Set to None such that -inf will be returned as ln_like
    #     wavel = None
    #     f_lambda = None
    #     contr_em = None
    #
    # else:
    #     wavel = 1e6 * constants.LIGHT / rt_object.freq  # (um)
    #
    #     # (erg s-1 cm-2 Hz-1) -> (erg s-1 m-2 Hz-1)
    #     f_lambda = 1e4 * rt_object.flux
    #
    #     # (erg s-1 m-2 Hz-1) -> (erg s-1 m-2 m-1)
    #     f_lambda *= constants.LIGHT / (1e-6 * wavel) ** 2.0
    #
    #     # (erg s-1 m-2 m-1) -> (erg s-1 m-2 um-1)
    #     f_lambda *= 1e-6
    #
    #     # (erg s-1 m-2 um-1) -> (W m-2 um-1)
    #     f_lambda *= 1e-7
    #
    #     # Optionally return the emission contribution
    #     if contribution:
    #         contr_em = rt_object.contr_em
    #     else:
    #         contr_em = None

    if rt_object.flux is None:
        wavel = None
        f_lambda = None
        contr_em = None

    else:
        wavel = 1e6 * constants.LIGHT / rt_object.freq  # (um)

        # (erg s-1 cm-2 Hz-1) -> (erg s-1 m-2 Hz-1)
        f_lambda = 1e4 * rt_object.flux

        # (erg s-1 m-2 Hz-1) -> (erg s-1 m-2 m-1)
        f_lambda *= constants.LIGHT / (1e-6 * wavel) ** 2.0

        # (erg s-1 m-2 m-1) -> (erg s-1 m-2 um-1)
        f_lambda *= 1e-6

        # (erg s-1 m-2 um-1) -> (W m-2 um-1)
        f_lambda *= 1e-7

        # Optionally return the emission contribution
        if contribution:
            contr_em = rt_object.contr_em
        else:
            contr_em = None

    # if plotting and Kzz_use is None and hasattr(rt_object, "continuum_opa"):
    #     plt.plot(wavel, rt_object.continuum_opa[:, 0], label="Total continuum opacity")
    #     plt.plot(
    #         wavel,
    #         rt_object.continuum_opa[:, 0] - rt_object.continuum_opa_scat[:, 0],
    #         label="Absorption continuum opacity",
    #     )
    #     plt.plot(
    #         wavel,
    #         rt_object.continuum_opa_scat[:, 0],
    #         label="Scattering continuum opacity",
    #     )
    #     plt.xlabel(r"Wavelength ($\mu$m)")
    #     plt.ylabel("Opacity at smallest pressure")
    #     plt.yscale("log")
    #     plt.legend(loc="best")
    #     plt.savefig("continuum_opacity.png", bbox_inches="tight")
    #     plt.clf()

    return wavel, f_lambda, contr_em, mmw


@typechecked
def mass_fractions(
    log_x_abund: Dict[str, float],
    line_species: List[str],
    abund_nodes: Optional[int] = None,
) -> Dict[str, float]:
    """
    Function to return a dictionary with the mass fractions of
    all species, including molecular hydrogen and helium.

    Parameters
    ----------
    log_x_abund : dict
        Dictionary with the log10 of the mass fractions of metals.
    line_species : list(str)
        List with the line species.
    abund_nodes : int, None
        Number of abundance nodes. The argument can be set to
        ``None`` if abundances are constant with pressure.

    Returns
    -------
    dict
        Dictionary with the mass fractions of all species.
    """

    if abund_nodes is None:
        # Initiate abundance dictionary
        abund = {}

        # Initiate the total mass fraction of the metals
        metal_sum = 0.0

        for line_item in line_species:
            if line_item in log_x_abund:
                # Add the mass fraction to the dictionary
                abund[line_item] = 10.0 ** log_x_abund[line_item]

                # Update the total mass fraction of the metals
                metal_sum += abund[line_item]

        # Mass fraction of H2 and He
        ab_h2_he = 1.0 - metal_sum

        # Add H2 and He mass fraction to the dictionary
        abund["H2"] = ab_h2_he * 0.75
        abund["He"] = ab_h2_he * 0.25

    else:
        # Initiate abundance dictionary
        abund = {}

        for node_idx in range(abund_nodes):
            # Initiate the total mass fraction of the metals
            metal_sum = 0.0

            for line_item in line_species:
                if f"{line_item}_{node_idx}" in log_x_abund:
                    # Add the mass fraction to the dictionary
                    abund[f"{line_item}_{node_idx}"] = (
                        10.0 ** log_x_abund[f"{line_item}_{node_idx}"]
                    )

                    # Update the total mass fraction of the metals
                    metal_sum += abund[f"{line_item}_{node_idx}"]

            # Mass fraction of H2 and He
            ab_h2_he = 1.0 - metal_sum

            # Add H2 and He mass fraction to the dictionary
            abund[f"H2_{node_idx}"] = ab_h2_he * 0.75
            abund[f"He_{node_idx}"] = ab_h2_he * 0.25

    return abund


@typechecked
def calc_metal_ratio(
    log_x_abund: Dict[str, float],
    line_species: List[str],
) -> Tuple[float, float, float]:
    """
    Function for calculating [C/H], [O/H], and C/O for a given set
    of abundances.

    Parameters
    ----------
    log_x_abund : dict
        Dictionary with the log10 mass fractions.
    line_species : list(str)
        List with the line species.

    Returns
    -------
    float
        Carbon-to-hydrogen ratio, relative to solar.
    float
        Oxygen-to-hydrogen ratio, relative to solar.
    float
        Carbon-to-oxygen ratio.
    """

    # Solar C/H from Asplund et al. (2009)

    c_h_solar = 10.0 ** (8.43 - 12.0)

    # Solar O/H from Asplund et al. (2009)

    o_h_solar = 10.0 ** (8.69 - 12.0)

    # Get the atomic masses

    masses = atomic_masses()

    # Create a dictionary with all mass fractions

    abund = mass_fractions(log_x_abund, line_species, abund_nodes=None)

    # Calculate the mean molecular weight from the input mass fractions

    mmw = mean_molecular_weight(abund)

    # Initiate the C, H, and O abundance

    c_abund = 0.0
    o_abund = 0.0
    h_abund = 0.0

    for abund_item in abund:
        abund_split = abund_item.split("_")[0]

        # Calculate the total C abundance

        if abund_split == "CO":
            c_abund += abund[abund_item] * mmw / masses["CO"]

        elif abund_split == "CO2":
            c_abund += abund[abund_item] * mmw / masses["CO2"]

        elif abund_split == "CH4":
            c_abund += abund[abund_item] * mmw / masses["CH4"]

        # Calculate the total O abundance

        if abund_split == "CO":
            o_abund += abund[abund_item] * mmw / masses["CO"]

        elif abund_split == "CO2":
            o_abund += 2.0 * abund[abund_item] * mmw / masses["CO2"]

        elif abund_split == "H2O":
            o_abund += abund[abund_item] * mmw / masses["H2O"]

        # Calculate the total H abundance

        if abund_split == "H2":
            h_abund += 2.0 * abund[abund_item] * mmw / masses["H2"]

        elif abund_split == "CH4":
            h_abund += 4.0 * abund[abund_item] * mmw / masses["CH4"]

        elif abund_split == "H2O":
            h_abund += 2.0 * abund[abund_item] * mmw / masses["H2O"]

        elif abund_split == "NH3":
            h_abund += 3.0 * abund[abund_item] * mmw / masses["NH3"]

        elif abund_split == "H2S":
            h_abund += 2.0 * abund[abund_item] * mmw / masses["H2S"]

    return (
        np.log10(c_abund / h_abund / c_h_solar),
        np.log10(o_abund / h_abund / o_h_solar),
        c_abund / o_abund,
    )


@typechecked
def mean_molecular_weight(mass_frac: Dict[str, float]) -> float:
    """
    Function to calculate the mean molecular weight from a
    dictionary with mass fractions.

    Parameters
    ----------
    mass_frac : dict
        Dictionary with the mass fraction of each species.

    Returns
    -------
    float
        Mean molecular weight in atomic mass units.
    """

    masses = atomic_masses()
    mmw_sum = 0.0

    for abund_item in mass_frac:
        if "_" in abund_item:
            mmw_sum += mass_frac[abund_item] / masses[abund_item.split("_")[0]]
        else:
            mmw_sum += mass_frac[abund_item] / masses[abund_item]

    return 1.0 / mmw_sum


@typechecked
def potassium_abundance(
    log_x_abund: Dict[str, float],
    line_species: List[str],
    abund_nodes: Optional[int] = None,
) -> Union[float, List[float]]:
    """
    Function to calculate the mass fraction of potassium at a solar
    ratio of the sodium and potassium abundances.

    Parameters
    ----------
    log_x_abund : dict
        Dictionary with the log10 of the mass fractions.
    line_species : list(str)
        List with the line species.
    abund_nodes : int, None
        Number of abundance nodes. The argument can be set to
        ``None`` if abundances are constant with pressure.

    Returns
    -------
    float, list(float)
        Log10 of the mass fraction of potassium or a list with the
        log10 mass fraction of potassium in case ``abund_nodes`` is
        not ``None``. The length of the list is equal to the argument
        of ``abund_nodes`` in that case.
    """

    # Solar volume mixing ratios of Na and K (Asplund et al. 2009)
    n_na_solar = 1.60008694353205e-06
    n_k_solar = 9.86605611925677e-08

    # Get the atomic masses
    masses = atomic_masses()

    # Create a dictionary with all mass fractions
    x_abund = mass_fractions(log_x_abund, line_species, abund_nodes)

    if abund_nodes is None:
        # Calculate mean molecular weight from the mass fractions

        mmw = mean_molecular_weight(x_abund)

        # Volume mixing ratio of sodium

        if "Na" in log_x_abund:
            n_na_abund = x_abund["Na"] * mmw / masses["Na"]

        elif "Na_lor_cut" in log_x_abund:
            n_na_abund = x_abund["Na_lor_cut"] * mmw / masses["Na"]

        elif "Na_allard" in log_x_abund:
            n_na_abund = x_abund["Na_allard"] * mmw / masses["Na"]

        elif "Na_burrows" in log_x_abund:
            n_na_abund = x_abund["Na_burrows"] * mmw / masses["Na"]

        # Volume mixing ratio of potassium

        n_k_abund = n_na_abund * n_k_solar / n_na_solar

        # Mass fraction of potassium
        log_x_k = np.log10(n_k_abund * masses["K"] / mmw)

    else:
        log_x_k = []

        for node_idx in range(abund_nodes):
            # Calculate mean molecular weight from the mass fractions

            abund_dict = {}
            for abund_item in line_species:
                abund_dict[abund_item] = x_abund[f"{abund_item}_{node_idx}"]

            mmw = mean_molecular_weight(abund_dict)

            # Volume mixing ratio of sodium

            if f"Na_{node_idx}" in log_x_abund:
                n_na_abund = x_abund[f"Na_{node_idx}"] * mmw / masses["Na"]

            elif f"Na_lor_cut_{node_idx}" in log_x_abund:
                n_na_abund = x_abund[f"Na_lor_cut_{node_idx}"] * mmw / masses["Na"]

            elif f"Na_allard_{node_idx}" in log_x_abund:
                n_na_abund = x_abund[f"Na_allard_{node_idx}"] * mmw / masses["Na"]

            elif f"Na_burrows_{node_idx}" in log_x_abund:
                n_na_abund = x_abund[f"Na_burrows_{node_idx}"] * mmw / masses["Na"]

            # Volume mixing ratio of potassium

            n_k_abund = n_na_abund * n_k_solar / n_na_solar

            # Mass fraction of potassium
            log_x_k.append(np.log10(n_k_abund * masses["K"] / mmw))

    return log_x_k


@typechecked
def log_x_cloud_base(
    c_o_ratio: float, metallicity: float, cloud_fractions: Dict[str, float]
) -> Dict[str, float]:
    """
    Function for returning a dictionary with the log10 mass fractions
    at the cloud base.

    Parameters
    ----------
    c_o_ratio : float
        C/O ratio.
    metallicity : float
        Metallicity, [Fe/H].
    cloud_fractions : dict
        Dictionary with the log10 mass fractions at the cloud base,
        relative to the maximum values allowed from elemental
        abundances. The dictionary keys are the cloud species without
        the structure and shape index (e.g. Na2S(c) instead of
        Na2S(c)_cd).

    Returns
    -------
    dict
        Dictionary with the log10 mass fractions at the cloud base.
        Compared to the keys of ``cloud_fractions``, the keys in the
        returned dictionary are provided without ``(c)`` (e.g. Na2S
        instead of Na2S(c)).
    """

    log_x_base = {}

    for item in cloud_fractions:
        # Mass fraction
        x_cloud = cloud_mass_fraction(f"{item[:-3]}", metallicity, c_o_ratio)

        # Log10 of the mass fraction at the cloud base
        log_x_base[f"{item[:-3]}"] = np.log10(10.0 ** cloud_fractions[item] * x_cloud)

    return log_x_base


@typechecked
def solar_mixing_ratios() -> dict:
    """
    Function which returns the volume mixing ratios for solar elemental
    abundances (i.e. [Fe/H] = 0); adopted from Asplund et al. (2009).

    Returns
    -------
    dict
        Dictionary with the solar number fractions (i.e. volume
        mixing ratios).
    """

    n_fracs = {}
    n_fracs["H"] = 0.9207539305
    n_fracs["He"] = 0.0783688694
    n_fracs["C"] = 0.0002478241
    n_fracs["N"] = 6.22506056949881e-05
    n_fracs["O"] = 0.0004509658
    n_fracs["Na"] = 1.60008694353205e-06
    n_fracs["Mg"] = 3.66558742055362e-05
    n_fracs["Al"] = 2.595e-06
    n_fracs["Si"] = 2.9795e-05
    n_fracs["P"] = 2.36670201997668e-07
    n_fracs["S"] = 1.2137900734604e-05
    n_fracs["Cl"] = 2.91167958499589e-07
    n_fracs["K"] = 9.86605611925677e-08
    n_fracs["Ca"] = 2.01439011429255e-06
    n_fracs["Ti"] = 8.20622804366359e-08
    n_fracs["V"] = 7.83688694089992e-09
    n_fracs["Fe"] = 2.91167958499589e-05
    n_fracs["Ni"] = 1.52807116806281e-06

    return n_fracs


@typechecked
def atomic_masses() -> dict:
    """
    Function which returns the atomic and molecular masses.

    Returns
    -------
    dict
        Dictionary with the atomic and molecular masses.
    """

    masses = {}

    # Atoms
    masses["H"] = 1.0
    masses["He"] = 4.0
    masses["C"] = 12.0
    masses["N"] = 14.0
    masses["O"] = 16.0
    masses["Na"] = 23.0
    masses["Na_lor_cur"] = 23.0
    masses["Na_allard"] = 23.0
    masses["Na_burrows"] = 23.0
    masses["Mg"] = 24.3
    masses["Al"] = 27.0
    masses["Si"] = 28.0
    masses["P"] = 31.0
    masses["S"] = 32.0
    masses["Cl"] = 35.45
    masses["K"] = 39.1
    masses["K_lor_cut"] = 39.1
    masses["K_allard"] = 39.1
    masses["K_burrows"] = 39.1
    masses["Ca"] = 40.0
    masses["Ti"] = 47.9
    masses["V"] = 51.0
    masses["Fe"] = 55.8
    masses["Ni"] = 58.7

    # Molecules
    masses["H2"] = 2.0
    masses["H2O"] = 18.0
    masses["H2O_HITEMP"] = 18.0
    masses["H2O_main_iso"] = 18.0
    masses["CH4"] = 16.0
    masses["CH4_main_iso"] = 16.0
    masses["CO2"] = 44.0
    masses["CO2_main_iso"] = 44.0
    masses["CO"] = 28.0
    masses["CO_all_iso"] = 28.0
    masses["CO_all_iso_Chubb"] = 28.0
    masses["CO_all_iso_HITEMP"] = 28.0
    masses["NH3"] = 17.0
    masses["NH3_main_iso"] = 17.0
    masses["HCN"] = 27.0
    masses["C2H2,acetylene"] = 26.0
    masses["PH3"] = 34.0
    masses["PH3_main_iso"] = 34.0
    masses["H2S"] = 34.0
    masses["H2S_main_iso"] = 34.0
    masses["VO"] = 67.0
    masses["VO_Plez"] = 67.0
    masses["TiO"] = 64.0
    masses["TiO_all_Exomol"] = 64.0
    masses["TiO_all_Plez"] = 64.0
    masses["FeH"] = 57.0
    masses["FeH_main_iso"] = 57.0
    masses["OH"] = 17.0

    return masses


@typechecked
def cloud_mass_fraction(
    composition: str, metallicity: float, c_o_ratio: float
) -> float:
    """
    Function to calculate the mass fraction for a cloud species.

    Parameters
    ----------
    composition : str
        Cloud composition ('Fe', 'MgSiO3', 'Mg2SiO4', 'Al2O3',
        'Na2S', or 'KCL').
    metallicity : float
        Metallicity [Fe/H].
    c_o_ratio : float
        Carbon-to-oxygen ratio.

    Returns
    -------
    float
        Mass fraction.
    """

    # Solar fractional number densities (i.e. volume mixing ratios; VMR)
    nfracs = solar_mixing_ratios()

    # Atomic masses
    masses = atomic_masses()

    # Make a copy of the dictionary with the solar number densities
    nfracs_use = copy.copy(nfracs)

    # Scale the solar number densities by the [Fe/H], except H and He
    for item in nfracs:
        if item != "H" and item != "He":
            nfracs_use[item] = nfracs[item] * 10.0**metallicity

    # Adjust the VMR of O with the C/O ratio
    nfracs_use["O"] = nfracs_use["C"] / c_o_ratio

    if composition == "Fe":
        nfrac_cloud = nfracs_use["Fe"]
        mass_cloud = masses["Fe"]

    elif composition == "MgSiO3":
        nfrac_cloud = np.min(
            [nfracs_use["Mg"], nfracs_use["Si"], nfracs_use["O"] / 3.0]
        )
        mass_cloud = masses["Mg"] + masses["Si"] + 3.0 * masses["O"]

    elif composition == "Mg2SiO4":
        nfrac_cloud = np.min(
            [nfracs_use["Mg"] / 2.0, nfracs_use["Si"], nfracs_use["O"] / 4.0]
        )
        mass_cloud = 2.0 * masses["Mg"] + masses["Si"] + 4.0 * masses["O"]

    elif composition == "Al2O3":
        nfrac_cloud = np.min([nfracs_use["Al"] / 2.0, nfracs_use["O"] / 3.0])
        mass_cloud = 2.0 * masses["Al"] + 3.0 * masses["O"]

    elif composition == "Na2S":
        nfrac_cloud = np.min([nfracs_use["Na"] / 2.0, nfracs_use["S"]])
        mass_cloud = 2.0 * masses["Na"] + masses["S"]

    elif composition == "KCL":
        nfrac_cloud = np.min([nfracs_use["K"], nfracs_use["Cl"]])
        mass_cloud = masses["K"] + masses["Cl"]

    # Cloud mass fraction
    x_cloud = mass_cloud * nfrac_cloud

    mass_norm = 0.0
    for item in nfracs_use:
        # Sum up the mass fractions of all species
        mass_norm += masses[item] * nfracs_use[item]

    # Normalize the cloud mass fraction by the total mass fraction
    return x_cloud / mass_norm


@typechecked
def get_condensation_curve(
    composition: str,
    press: np.ndarray,
    metallicity: float,
    c_o_ratio: float,
    mmw: float = 2.33,
) -> np.ndarray:
    """
    Function to find the base of the cloud deck by intersecting the
    P-T profile with the saturation vapor pressure.

    Parameters
    ----------
    composition : str
        Cloud composition ('Fe', 'MgSiO3', 'Mg2SiO4', 'Al2O3',
        'Na2S', or 'KCL').
    press : np.ndarray
        Pressures (bar).
    metallicity : float
        Metallicity [Fe/H].
    c_o_ratio : float
        Carbon-to-oxygen ratio.
    mmw : float
        Mean molecular weight.

    Returns
    -------
    np.array
        Condensation temperatures (K) for the provided input pressures.
    """

    if composition == "Fe":
        Pc, Tc = return_T_cond_Fe_comb(metallicity, c_o_ratio, mmw)

    elif composition == "MgSiO3":
        Pc, Tc = return_T_cond_MgSiO3(metallicity, c_o_ratio, mmw)

    elif composition == "Mg2SiO4":
        Pc, Tc = return_T_cond_Mg2SiO4(metallicity)

    elif composition == "Al2O3":
        Pc, Tc = return_T_cond_Al2O3(metallicity)

    elif composition == "Na2S":
        Pc, Tc = return_T_cond_Na2S(metallicity, c_o_ratio, mmw)

    elif composition == "KCL":
        Pc, Tc = return_T_cond_KCl(metallicity, c_o_ratio, mmw)

    else:
        raise ValueError(
            f"The '{composition}' composition is not "
            "supported by get_condensation_curve."
        )

    index = (Pc > 1e-8) & (Pc < 1e5)
    Pc, Tc = Pc[index], Tc[index]

    tcond_p = interp1d(Pc, Tc)

    return tcond_p(press)


@typechecked
def find_cloud_deck(
    composition: str,
    press: np.ndarray,
    temp: np.ndarray,
    metallicity: float,
    c_o_ratio: float,
    mmw: float = 2.33,
    plotting: bool = False,
) -> float:
    """
    Function to find the base of the cloud deck by intersecting the
    P-T profile with the saturation vapor pressure.

    Parameters
    ----------
    composition : str
        Cloud composition ('Fe', 'MgSiO3', 'Mg2SiO4', 'Al2O3',
        'Na2S', or 'KCL').
    press : np.ndarray
        Pressures (bar).
    temp : np.ndarray
        Temperatures (K).
    metallicity : float
        Metallicity [Fe/H].
    c_o_ratio : float
        Carbon-to-oxygen ratio.
    mmw : float
        Mean molecular weight.
    plotting : bool
        Create a plot.

    Returns
    -------
    float
        Pressure (bar) at the base of the cloud deck.
    """

    Tcond_on_input_grid = get_condensation_curve(
        composition=composition,
        press=press,
        metallicity=metallicity,
        c_o_ratio=c_o_ratio,
        mmw=mmw,
    )

    Tdiff = Tcond_on_input_grid - temp
    diff_vec = Tdiff[1:] * Tdiff[:-1]
    ind_cdf = diff_vec < 0.0

    if len(diff_vec[ind_cdf]) > 0:
        P_clouds = (press[1:] + press[:-1])[ind_cdf] / 2.0
        P_cloud = float(P_clouds[-1])

    else:
        P_cloud = 1e-8

    if plotting:
        plt.plot(temp, press)
        plt.plot(Tcond_on_input_grid, press)
        plt.axhline(P_cloud, color="red", linestyle="--")
        plt.yscale("log")
        plt.xlim(0.0, 3000.0)
        plt.ylim(1e2, 1e-6)
        plt.savefig(f"{composition.lower()}_clouds_cdf.png", bbox_inches="tight")
        plt.clf()

    return P_cloud


@typechecked
def scale_cloud_abund(
    params: Dict[str, float],
    rt_object,
    pressure: np.ndarray,
    temperature: np.ndarray,
    mmw: np.ndarray,
    chemistry: str,
    abund_in: Dict[str, np.ndarray],
    composition: str,
    tau_cloud: float,
    pressure_grid: str,
) -> float:
    """
    Function to scale the mass fraction of a cloud species to the
    requested optical depth.

    Parameters
    ----------
    params : dict
        Dictionary with the model parameters.
    rt_object : petitRADTRANS.radtrans.Radtrans
        Instance of ``Radtrans``.
    pressure : np.ndarray
        Array with the pressure points (bar).
    temperature : np.ndarray
        Array with the temperature points (K) corresponding
        to ``pressure``.
    mmw : np.ndarray
        Array with the mean molecular weights corresponding
        to ``pressure``.
    chemistry : str
        Chemistry type (only ``'equilibrium'`` is supported).
    abund_in : dict
        Dictionary with arrays that contain the pressure-dependent,
        equilibrium mass fractions of the line species.
    composition : sr
        Cloud composition ('Fe(c)', 'MgSiO3(c)', 'Mg2SiO4(c)',
        'Al2O3(c)', 'Na2S(c)', 'KCl(c)').
    tau_cloud : float
        Optical depth of the clouds. The returned mass fraction is
        scaled such that the optical depth at the shortest wavelength
        is equal to ``tau_cloud``.
    pressure_grid : str
        The type of pressure grid that is used for the radiative
        transfer. Either 'standard', to use 180 layers both for the
        atmospheric structure (e.g. when interpolating the abundances)
        and 180 layers with the radiative transfer, or 'smaller' to
        use 60 (instead of 180) with the radiative transfer, or
        'clouds' to start with 1440 layers but resample to ~100 layers
        (depending on the number of cloud species) with a refinement
        around the cloud decks. For cloudless atmospheres it is
        recommended to use 'smaller', which runs faster than 'standard'
        and provides sufficient accuracy. For cloudy atmosphere, one
        can test with 'smaller' but it is recommended to use 'clouds' for
        improved accuracy fluxes.

    Returns
    -------
    float
        Mass fraction relative to the maximum value allowed from
        elemental abundances. The value has been scaled to the
        requested optical depth ``tau_cloud`` (at the shortest
        wavelength).
    """

    # Dictionary with the requested cloud composition and setting the
    # log10 of the mass fraction (relative to the maximum value
    # allowed from elemental abundances) equal to zero
    cloud_fractions = {composition: 0.0}

    # Create a dictionary with the log10 of
    # the mass fraction at the cloud base
    log_x_base = log_x_cloud_base(
        params["c_o_ratio"], params["metallicity"], cloud_fractions
    )

    # Get the pressure (bar) of the cloud base
    p_base = find_cloud_deck(
        composition[:-3],
        pressure,
        temperature,
        params["metallicity"],
        params["c_o_ratio"],
        mmw=np.mean(mmw),
        plotting=False,
    )

    # Initialize the cloud abundance in
    # the dictionary with mass fractions
    abund_in[composition] = np.zeros_like(temperature)

    # Set the cloud abundances by scaling
    # from the base with the f_sed parameter
    abund_in[composition][pressure < p_base] = (
        10.0 ** log_x_base[composition[:-3]]
        * (pressure[pressure <= p_base] / p_base) ** params["fsed"]
    )

    # Adaptive pressure refinement around the cloud base
    if pressure_grid == "clouds":
        _, indices = make_half_pressure_better({composition: p_base}, pressure)
    else:
        indices = None

    # Update the abundance dictionary

    abundances = create_abund_dict(
        abund_in,
        rt_object.line_species,
        chemistry,
        pressure_grid=pressure_grid,
        indices=indices,
    )

    # Interpolate the line opacities to the temperature structure

    if pressure_grid == "standard":
        rt_object.interpolate_species_opa(temperature)

        mmw_select = mmw.copy()

        if "log_kzz" in params:
            kzz_select = np.full(pressure.size, 10.0 ** params["log_kzz"])
        else:
            # Backward compatibility
            kzz_select = np.full(pressure.size, 10.0 ** params["kzz"])

    elif pressure_grid == "smaller":
        rt_object.interpolate_species_opa(temperature[::3])

        mmw_select = mmw[::3]

        if "log_kzz" in params:
            kzz_select = np.full(pressure[::3].size, 10.0 ** params["log_kzz"])
        else:
            # Backward compatibility
            kzz_select = np.full(pressure[::3].size, 10.0 ** params["kzz"])

    elif pressure_grid == "clouds":
        # Reinitiate the pressure structure
        # after make_half_pressure_better
        rt_object.setup_opa_structure(pressure[indices])
        rt_object.interpolate_species_opa(temperature[indices])

        mmw_select = mmw[indices]

        if "log_kzz" in params:
            kzz_select = np.full(pressure[indices].size, 10.0 ** params["log_kzz"])
        else:
            # Backward compatibility
            kzz_select = np.full(pressure[indices].size, 10.0 ** params["kzz"])

    # Set the continuum opacities to zero because
    # calc_cloud_opacity adds to existing opacities
    rt_object.continuum_opa = np.zeros_like(rt_object.continuum_opa)
    rt_object.continuum_opa_scat = np.zeros_like(rt_object.continuum_opa_scat)
    rt_object.continuum_opa_scat_emis = np.zeros_like(rt_object.continuum_opa_scat_emis)

    # Calculate the cloud opacities for
    # the defined atmospheric structure
    rt_object.calc_cloud_opacity(
        abundances,
        mmw_select,
        10.0 ** params["logg"],
        params["sigma_lnorm"],
        fsed=params["fsed"],
        Kzz=kzz_select,
        radius=None,
        add_cloud_scat_as_abs=False,
    )

    # Calculate the cloud optical depth and set the tau_cloud attribute
    rt_object.calc_tau_cloud(10.0 ** params["logg"])

    # Extract the wavelength-averaged optical
    # depth at the largest pressure
    tau_current = np.mean(rt_object.tau_cloud[0, :, 0, -1])

    # Set the continuum opacities again to zero
    rt_object.continuum_opa = np.zeros_like(rt_object.continuum_opa)
    rt_object.continuum_opa_scat = np.zeros_like(rt_object.continuum_opa_scat)
    rt_object.continuum_opa_scat_emis = np.zeros_like(rt_object.continuum_opa_scat_emis)

    if tau_current > 0.0:
        # Scale the mass fraction
        log_x_scaled = np.log10(tau_cloud / tau_current)
    else:
        log_x_scaled = 100.0

    return log_x_scaled


@typechecked
def cube_to_dict(cube, cube_index: Dict[str, float]) -> Dict[str, float]:
    """
    Function to convert the parameter cube into a dictionary.

    Parameters
    ----------
    cube : LP_c_double
        Cube with the parameters.
    cube_index : dict
        Dictionary with the index of each parameter in the ``cube``.

    Returns
    -------
    dict
        Dictionary with the parameters.
    """

    params = {}

    for key, value in cube_index.items():
        params[key] = cube[value]

    return params


@typechecked
def list_to_dict(param_list: List[str], sample_val: np.ndarray) -> Dict[str, float]:
    """
    Function to convert the parameter cube into a dictionary.

    Parameters
    ----------
    param_list : list(str)
        List with the parameter labels.
    sample_val : np.ndarray
        Array with the parameter values, in the same order as
        ``param_list``.

    Returns
    -------
    dict
        Dictionary with the parameters.
    """

    sample_dict = {}

    for item in param_list:
        sample_dict[item] = sample_val[param_list.index(item)]

    return sample_dict


@typechecked
def return_T_cond_Fe(
    FeH: float, CO: float, MMW: float = 2.33
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for calculating the saturation pressure for solid Fe.

    Parameters
    ----------
    FeH : float
        Metallicity.
    CO : float
        Carbon-to-oxygen ratio.
    MMW : float
        Mean molecular weight.

    Returns
    -------
    np.ndarray
        Saturation pressure (bar).
    np.ndarray
        Temperature (K).
    """

    masses = atomic_masses()

    T = np.linspace(100.0, 10000.0, 1000)

    # Taken from Ackerman & Marley (2001)
    # including erratum (P_vap is in bar, not cgs!)
    P_vap = lambda x: np.exp(15.71 - 47664.0 / x)

    x_cloud = cloud_mass_fraction("Fe", FeH, CO)

    return P_vap(T) / (x_cloud * MMW / masses["Fe"]), T


@typechecked
def return_T_cond_Fe_l(
    FeH: float, CO: float, MMW: float = 2.33
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for calculating the saturation pressure for liquid Fe.

    Parameters
    ----------
    FeH : float
        Metallicity.
    CO : float
        Carbon-to-oxygen ratio.
    MMW : float
        Mean molecular weight.

    Returns
    -------
    np.ndarray
        Saturation pressure (bar).
    np.ndarray
        Temperature (K).
    """

    masses = atomic_masses()

    T = np.linspace(100.0, 10000.0, 1000)

    # Taken from Ackerman & Marley (2001)
    # including erratum (P_vap is in bar, not cgs!)
    P_vap = lambda x: np.exp(9.86 - 37120.0 / x)

    x_cloud = cloud_mass_fraction("Fe", FeH, CO)

    return P_vap(T) / (x_cloud * MMW / masses["Fe"]), T


@typechecked
def return_T_cond_Fe_comb(
    FeH: float, CO: float, MMW: float = 2.33
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for calculating the saturation pressure for Fe.

    Parameters
    ----------
    FeH : float
        Metallicity.
    CO : float
        Carbon-to-oxygen ratio.
    MMW : float
        Mean molecular weight.

    Returns
    -------
    np.ndarray
        Saturation pressure (bar).
    np.ndarray
        Temperature (K).
    """

    P1, T1 = return_T_cond_Fe(FeH, CO, MMW)
    P2, T2 = return_T_cond_Fe_l(FeH, CO, MMW)

    retP = np.zeros_like(P1)
    index = P1 < P2
    retP[index] = P1[index]
    retP[~index] = P2[~index]

    return retP, T2


@typechecked
def return_T_cond_MgSiO3(
    FeH: float, CO: float, MMW: float = 2.33
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for calculating the saturation pressure for MgSiO3.

    Parameters
    ----------
    FeH : float
        Metallicity.
    CO : float
        Carbon-to-oxygen ratio.
    MMW : float
        Mean molecular weight.

    Returns
    -------
    np.ndarray
        Saturation pressure (bar).
    np.ndarray
        Temperature (K).
    """

    masses = atomic_masses()

    T = np.linspace(100.0, 10000.0, 1000)

    # Taken from Ackerman & Marley (2001)
    # including erratum (P_vap is in bar, not cgs!)
    P_vap = lambda x: np.exp(25.37 - 58663.0 / x)

    x_cloud = cloud_mass_fraction("MgSiO3", FeH, CO)

    m_mgsio3 = masses["Mg"] + masses["Si"] + 3.0 * masses["O"]

    return P_vap(T) / (x_cloud * MMW / m_mgsio3), T


@typechecked
def return_T_cond_Mg2SiO4(FeH: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for calculating the saturation pressure for Mg2SiO4.

    Parameters
    ----------
    FeH : float
        Metallicity.

    Returns
    -------
    np.ndarray
        Saturation pressure (bar).
    np.ndarray
        Temperature (K).
    """

    # Condensation temperature of Mg2SiO4
    # See Eq. 18 in Visscher et al. 2010

    temp = np.linspace(100.0, 10000.0, 1000)
    P_vap = lambda x: 10.0 ** ((5.89 - 1e4 / temp - 0.73 * FeH) / 0.37)

    return P_vap(temp), temp


@typechecked
def return_T_cond_Al2O3(FeH: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for calculating the condensation temperature for Al2O3.

    Parameters
    ----------
    FeH : float
        Metallicity.
    CO : float
        Carbon-to-oxygen ratio.
    MMW : float
        Mean molecular weight.

    Returns
    -------
    np.ndarray
        Saturation pressure (bar).
    np.ndarray
        Temperature (K).
    """

    # Return dictionary with atomic masses
    # masses = atomic_masses()

    # Create pressures (bar)
    pressure = np.logspace(-6, 3, 1000)

    # Equilibrium mass fraction of Al2O3
    # x_cloud = cloud_mass_fraction('Al2O3', FeH, CO)

    # Molecular mass of Al2O3
    # m_al2o3 = 3. * masses['Al'] + 2. * masses['O']

    # Partial pressure of Al2O3
    # part_press = pressure/(x_cloud*MMW/m_al2o3)

    # Condensation temperature of Al2O3
    # (see Eq. 4 in Wakeford et al. 2017)
    t_cond = 1e4 / (
        5.014
        - 0.2179 * np.log10(pressure)
        + 2.264e-3 * np.log10(pressure) ** 2
        - 0.580 * FeH
    )

    return pressure, t_cond


@typechecked
def return_T_cond_Na2S(
    FeH: float, CO: float, MMW: float = 2.33
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for calculating the saturation pressure for Na2S.

    Parameters
    ----------
    FeH : float
        Metallicity.
    CO : float
        Carbon-to-oxygen ratio.
    MMW : float
        Mean molecular weight.

    Returns
    -------
    np.ndarray
        Saturation pressure (bar).
    np.ndarray
        Temperature (K).
    """

    masses = atomic_masses()

    # Taken from Charnay+2018
    T = np.linspace(100.0, 10000.0, 1000)

    # This is the partial pressure of Na, so divide by factor 2 to get
    # the partial pressure of the hypothetical Na2S gas particles, this
    # is OK: there are more S than Na atoms at solar abundance ratios.
    P_vap = lambda x: 1e1 ** (8.55 - 13889.0 / x - 0.5 * FeH) / 2.0

    x_cloud = cloud_mass_fraction("Na2S", FeH, CO)

    m_na2s = 2.0 * masses["Na"] + masses["S"]

    return P_vap(T) / (x_cloud * MMW / m_na2s), T


@typechecked
def return_T_cond_KCl(
    FeH: float, CO: float, MMW: float = 2.33
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function for calculating the saturation pressure for KCl.

    Parameters
    ----------
    FeH : float
        Metallicity.
    CO : float
        Carbon-to-oxygen ratio.
    MMW : float
        Mean molecular weight.

    Returns
    -------
    np.ndarray
        Saturation pressure (bar).
    np.ndarray
        Temperature (K).
    """

    masses = atomic_masses()

    T = np.linspace(100.0, 10000.0, 1000)

    # Taken from Charnay+2018
    P_vap = lambda x: 1e1 ** (7.611 - 11382.0 / T)

    x_cloud = cloud_mass_fraction("KCL", FeH, CO)

    m_kcl = masses["K"] + masses["Cl"]

    return P_vap(T) / (x_cloud * MMW / m_kcl), T


@typechecked
def convolve_spectrum(
    input_wavel: np.ndarray, input_flux: np.ndarray, spec_res: float
) -> np.ndarray:
    """
    Function to convolve a spectrum with a Gaussian filter.

    Parameters
    ----------
    input_wavel : np.ndarray
        Input wavelengths.
    input_flux : np.ndarray
        Input flux
    spec_res : float
        Spectral resolution of the Gaussian filter.

    Returns
    -------
    np.ndarray
        Convolved spectrum.
    """

    # delta_lambda of resolution element is
    # FWHM of the LSF's standard deviation
    sigma_lsf = 1.0 / spec_res / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    # The input spacing of petitRADTRANS is 1e3, but just compute
    # it to be sure, or more versatile in the future.
    # Also, we have a log-spaced grid, so the spacing is constant
    # as a function of wavelength
    spacing = np.mean(2.0 * np.diff(input_wavel) / (input_wavel[1:] + input_wavel[:-1]))

    # Calculate the sigma to be used with the Gaussian filter
    # in units of the input wavelength bins
    sigma_filter = sigma_lsf / spacing

    return gaussian_filter(input_flux, sigma=sigma_filter, mode="nearest")


@typechecked
def quench_pressure(
    pressure: np.ndarray,
    temperature: np.ndarray,
    metallicity: float,
    c_o_ratio: float,
    log_g: float,
    log_kzz: float,
) -> Optional[float]:
    """
    Function to determine the CO/CH$_4$ quenching pressure by
    intersecting the pressure-dependent timescales of the
    vertical mixing and the CO/CH$_4$ reaction rates.

    Parameters
    ----------
    pressure : np.ndarray
        Array with the pressures (bar).
    temperature : np.ndarray
        Array with the temperatures (K) corresponding to ``pressure``.
    metallicity : float
        Metallicity [Fe/H].
    c_o_ratio : float
        Carbon-to-oxygen ratio.
    log_g : float
        Log10 of the surface gravity (cm s-2).
    log_kzz : float
        Log10 of the eddy diffusion coefficient (cm2 s-1).

    Returns
    -------
    float, None
        Quenching pressure (bar).
    """

    # Interpolate the equilibbrium abundances

    co_array = np.full(pressure.shape[0], c_o_ratio)
    feh_array = np.full(pressure.shape[0], metallicity)

    if "poor_mans_nonequ_chem" in sys.modules:
        from poor_mans_nonequ_chem.poor_mans_nonequ_chem import interpol_abundances
    else:
        from petitRADTRANS.poor_mans_nonequ_chem.poor_mans_nonequ_chem import (
            interpol_abundances,
        )

    abund_eq = interpol_abundances(
        co_array, feh_array, temperature, pressure, Pquench_carbon=None
    )

    # Surface gravity (m s-2)
    gravity = 1e-2 * 10.0**log_g

    # Mean molecular weight (kg)
    mmw = abund_eq["MMW"] * constants.ATOMIC_MASS

    # Pressure scale height (m)
    h_scale = constants.BOLTZMANN * temperature / (mmw * gravity)

    # Diffusion coefficient (m2 s-1)
    kzz = 1e-4 * 10.0**log_kzz

    # Mixing timescale (s)
    t_mix = h_scale**2 / kzz

    # Chemical timescale (see Eq. 12 from Zahnle & Marley 2014)
    metal = 10.0**metallicity
    t_chem = 1.5e-6 * pressure**-1.0 * metal**-0.7 * np.exp(42000.0 / temperature)

    # Determine pressure at which t_mix = t_chem

    t_diff = t_mix - t_chem
    diff_product = t_diff[1:] * t_diff[:-1]

    # If t_mix and t_chem intersect then there
    # is 1 negative value in diff_product
    indices = diff_product < 0.0

    if np.sum(indices) == 1:
        p_quench = (pressure[1:] + pressure[:-1])[indices] / 2.0
        p_quench = p_quench[0]

    elif np.sum(indices) == 0:
        p_quench = None

    else:
        raise ValueError(
            f"Encountered unexpected number of indices "
            f"({np.sum(indices)}) when determining the "
            f"intersection of t_mix and t_chem."
        )

    return p_quench


def convective_flux(
    press: np.ndarray,
    temp: np.ndarray,
    mmw: np.ndarray,
    nabla_ad: np.ndarray,
    kappa_r: np.ndarray,
    density: np.ndarray,
    c_p: np.ndarray,
    gravity: float,
    f_bol: float,
    mix_length: float = 1.0,
) -> np.ndarray:
    """
    Function for calculating the convective flux with mixing-length
    theory. This function has been adopted from petitCODE (Paul
    Mollière, MPIA) and was converted from Fortran to Python.

    Parameters
    ----------
    press : np.ndarray
        Array with the pressures (Pa).
    temp : np.ndarray
        Array with the temperatures (K) at ``pressure``.
    mmw : np.ndarray
        Array with the mean molecular weights at ``pressure``.
    nabla_ad : np.ndarray
        Array with the adiabatic temperature gradient at ``pressure``.
    kappa_r : np.ndarray
        Array with the Rosseland mean opacity (m2 kg-1) at
        ``pressure``.
    density : np.ndarray
        Array with the density (kg m-3) at ``pressure``.
    c_p : np.ndarray
        Array with the specific heat capacity (J kg-1 K-1) at
        constant pressure, ``pressure``.
    gravity : float
        Surface gravity (m s-2).
    f_bol : float
        Bolometric flux (W m-2) at the top of the atmosphere,
        calculated from the low-resolution spectrum.
    mix_length : float
        Mixing length for the convection in units of the pressure
        scale height (default: 1.0).

    Returns
    -------
    np.ndarray
        Convective flux (W m-2) at each pressure.
    """

    t_transp = (f_bol / constants.SIGMA_SB) ** 0.25  # (K)
    nabla_rad = (
        3.0 * kappa_r * press * t_transp**4.0 / 16.0 / gravity / temp**4.0
    )  # (dimensionless)
    h_press = (
        constants.BOLTZMANN * temp / (mmw * constants.ATOMIC_MASS * gravity)
    )  # (m)
    l_mix = mix_length * h_press  # (m)

    U = (
        (12.0 * constants.SIGMA_SB * temp**3.0)
        / (c_p * density**2.0 * kappa_r * l_mix**2.0)
        * np.sqrt(8.0 * h_press / gravity)
    )

    W = nabla_rad - nabla_ad

    # FIXME thesis: 2336U^4W
    A = (
        1168.0 * U**3.0
        + 2187 * U * W
        + 27.0
        * np.sqrt(
            3.0 * (2048.0 * U**6.0 + 2236.0 * U**4.0 * W + 2187.0 * U**2.0 * W**2.0)
        )
    ) ** (1.0 / 3.0)

    xi = (
        19.0 / 27.0 * U
        - 184.0 / 27.0 * 2.0 ** (1.0 / 3.0) * U**2.0 / A
        + 2.0 ** (2.0 / 3.0) / 27.0 * A
    )

    nabla = xi**2.0 + nabla_ad - U**2.0
    nabla_e = nabla_ad + 2.0 * U * xi - 2.0 * U**2.0

    f_conv = (
        density
        * c_p
        * temp
        * np.sqrt(gravity)
        * (mix_length * h_press) ** 2.0
        / (4.0 * np.sqrt(2.0))
        * h_press**-1.5
        * (nabla - nabla_e) ** 1.5
    )

    f_conv[np.isnan(f_conv)] = 0.0

    return f_conv  # (W m-2)
