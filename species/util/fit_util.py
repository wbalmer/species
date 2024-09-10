"""
Utility functions for fit results.
"""

import os
import warnings

from configparser import ConfigParser
from typing import Dict, List, Optional, Union

import h5py
import numpy as np
import spectres

from typeguard import typechecked

from species.core.box import ObjectBox, ResidualsBox, SynphotBox, create_box
from species.phot.syn_phot import SyntheticPhotometry
from species.read.read_calibration import ReadCalibration
from species.read.read_filter import ReadFilter
from species.read.read_model import ReadModel
from species.read.read_planck import ReadPlanck
from species.read.read_radtrans import ReadRadtrans
from species.util.core_util import print_section
from species.util.model_util import binary_to_single, powerlaw_spectrum
from species.util.retrieval_util import convolve_spectrum


@typechecked
def multi_photometry(
    datatype: str,
    spectrum: str,
    filters: List[str],
    parameters: Dict[str, float],
    radtrans: Optional[ReadRadtrans] = None,
    verbose: bool = True,
) -> SynphotBox:
    """
    Function for calculating synthetic photometry for a list of
    filters and a specified atmosphere model and related parameters.
    This function can for example be used for calculating the
    synthetic photometry from a best-fit model spectrum. It returns
    a :class:`~species.core.box_types.SynphotBox` that can be provided
    as input to :func:`~species.plot.plot_spectrum.plot_spectrum`.

    Parameters
    ----------
    datatype : str
        Data type ('model' or 'calibration').
    spectrum : str
        Spectrum name (e.g., 'drift-phoenix', 'bt-settl-cifist',
        planck', 'powerlaw', 'petitradtrans').
    filters : list(str)
        List with the filter names.
    parameters : dict
        Dictionary with the model parameters.
    radtrans : ReadRadtrans, None
        Instance of :class:`~species.read.read_radtrans.ReadRadtrans`.
        Only required with ``spectrum='petitradtrans'`. Make sure that
        the ``wavel_range`` of the ``ReadRadtrans`` instance is
        sufficiently broad to cover all the ``filters``. The argument
        can be set to ``None`` for any other model than petitRADTRANS.
    verbose : bool
        Print output.

    Returns
    -------
    species.core.box.SynphotBox
        Box with synthetic photometry.
    """

    if verbose:
        print_section("Calculate multi-photometry")

        print(f"Data type: {datatype}")
        print(f"Spectrum name: {spectrum}")

        print("\nParameters:")
        for param_key, param_value in parameters.items():
            if -0.1 < param_value < 0.1:
                print(f"   - {param_key} = {param_value:.2e}")
            else:
                print(f"   - {param_key} = {param_value:.2f}")

    mean_wavel = {}

    for filter_item in filters:
        read_filt = ReadFilter(filter_item)
        mean_wavel[filter_item] = read_filt.mean_wavelength()

    flux = {}

    if datatype == "model":
        if spectrum == "petitradtrans":
            # Calculate the petitRADTRANS spectrum only once
            radtrans_box = radtrans.get_model(parameters)

        for item in filters:
            if spectrum == "petitradtrans":
                # Use an instance of SyntheticPhotometry instead
                # of get_flux from ReadRadtrans in order to not
                # recalculate the spectrum
                syn_phot = SyntheticPhotometry(item)

                flux[item], _ = syn_phot.spectrum_to_flux(
                    radtrans_box.wavelength, radtrans_box.flux
                )

            elif spectrum == "powerlaw":
                synphot = SyntheticPhotometry(item)

                # Set the wavel_range attribute
                synphot.zero_point()

                powerl_box = powerlaw_spectrum(synphot.wavel_range, parameters)
                flux[item] = synphot.spectrum_to_flux(
                    powerl_box.wavelength, powerl_box.flux
                )[0]

            else:
                if spectrum == "planck":
                    readmodel = ReadPlanck(filter_name=item)

                else:
                    readmodel = ReadModel(spectrum, filter_name=item)

                try:
                    if (
                        spectrum != "planck"
                        and "teff_0" in parameters
                        and "teff_1" in parameters
                    ):
                        # Binary system

                        param_0 = binary_to_single(parameters, 0)
                        model_flux_0 = readmodel.get_flux(param_0)[0]

                        param_1 = binary_to_single(parameters, 1)
                        model_flux_1 = readmodel.get_flux(param_1)[0]

                        # Weighted flux of two spectra for atmospheric asymmetries
                        # Or simply the same in case of an actual binary system

                        if "spec_weight" in parameters:
                            flux[item] = (
                                parameters["spec_weight"] * model_flux_0
                                + (1.0 - parameters["spec_weight"]) * model_flux_1
                            )

                        else:
                            flux[item] = model_flux_0 + model_flux_1

                    else:
                        # Single object

                        flux[item] = readmodel.get_flux(parameters)[0]

                except IndexError:
                    flux[item] = np.nan

                    warnings.warn(
                        f"The wavelength range of the {item} filter does not "
                        f"match with the wavelength range of {spectrum}. The "
                        f"flux is set to NaN."
                    )

    elif datatype == "calibration":
        for item in filters:
            readcalib = ReadCalibration(spectrum, filter_name=item)
            flux[item] = readcalib.get_flux(parameters)[0]

    app_mag = {}
    abs_mag = {}

    for key, value in flux.items():
        syn_phot = SyntheticPhotometry(key)
        if "parallax" in parameters:
            app_mag[key], abs_mag[key] = syn_phot.flux_to_magnitude(
                flux=value, error=None, parallax=(parameters["parallax"], None)
            )

        elif "distance" in parameters:
            app_mag[key], abs_mag[key] = syn_phot.flux_to_magnitude(
                flux=value, error=None, distance=(parameters["distance"], None)
            )

        else:
            app_mag[key], abs_mag[key] = syn_phot.flux_to_magnitude(
                flux=value, error=None
            )

    if verbose:
        print("\nMagnitudes:")
        for key, value in app_mag.items():
            if value[1] is None:
                print(f"   - {key} = {value[0]:.2f}")
            else:
                print(f"   - {key} = {value[0]:.2f} +/- {value[1]:.2f}")

        print("\nFluxes (W m-2 um-1):")
        for key, value in flux.items():
            print(f"   - {key} = {value:.2e}")

    return create_box(
        "synphot",
        name="synphot",
        flux=flux,
        wavelength=mean_wavel,
        app_mag=app_mag,
        abs_mag=abs_mag,
    )


@typechecked
def get_residuals(
    tag: str,
    parameters: Dict[str, float],
    objectbox: ObjectBox,
    inc_phot: Union[bool, List[str]] = True,
    inc_spec: Union[bool, List[str]] = True,
    radtrans: Optional[ReadRadtrans] = None,
    datatype: Optional[str] = None,
    spectrum: Optional[str] = None,
) -> ResidualsBox:
    """
    Function for calculating the residuals from fitting model or
    calibration spectra to a set of spectra and/or photometry.

    Parameters
    ----------
    tag: str
        Database tag with the sampling results.
    parameters : dict
        Parameters and values for the spectrum
    objectbox : ObjectBox
        Box with the photometry and/or spectra of an object. A scaling
        and/or error inflation of the spectra should be applied with
        :func:`~species.util.read_util.update_objectbox` beforehand.
    inc_phot : bool, list(str)
        Include photometric data in the fit. If a boolean, either all
        (``True``) or none (``False``) of the data are selected. If a
        list, a subset of filter names (as stored in the database) can
        be provided.
    inc_spec : bool, list(str)
        Include spectroscopic data in the fit. If a boolean, either all
        (``True``) or none (``False``) of the data are selected. If a
        list, a subset of spectrum names (as stored in the database
        with :func:`~species.data.database.Database.add_object`) can be
        provided.
    radtrans : ReadRadtrans, None
        Instance of :class:`~species.read.read_radtrans.ReadRadtrans`.
        Only required with ``spectrum='petitradtrans'``. Make sure that
        the ``wavel_range`` of the ``ReadRadtrans`` instance is
        sufficiently broad to cover all the photometric and
        spectroscopic data of ``inc_phot`` and ``inc_spec``. Not used
        if the argument is set to ``None``.

    Returns
    -------
    species.core.box.ResidualsBox
        Box with the residuals.
    """

    # Check deprecated parameters

    if datatype is not None:
        warnings.warn(
            "The 'datatype' parameter is no longer "
            "used by the 'get_residuals' function. "
            "Instead, the 'tag' parameter should be set, "
            "which points to the sampling results as "
            "stored in the database.",
            DeprecationWarning,
        )

    if spectrum is not None:
        warnings.warn(
            "The 'spectrum' parameter is no longer "
            "used by the 'get_residuals' function. "
            "Instead, the 'tag' parameter should be set, "
            "which points to the sampling results as "
            "stored in the database.",
            DeprecationWarning,
        )

    # Read sampling results

    print_section("Calculate residuals")

    config_file = os.path.join(os.getcwd(), "species_config.ini")

    config = ConfigParser()
    config.read(config_file)

    database_path = config["species"]["database"]

    with h5py.File(database_path, "r") as hdf5_file:
        dset = hdf5_file[f"results/fit/{tag}/samples"]
        spectrum = dset.attrs["spectrum"]
        binary = dset.attrs["binary"]
        print(f"Model: {spectrum}")
        print(f"Binary: {binary}")

        n_param = dset.attrs["n_param"]

        if "n_fixed" in dset.attrs:
            n_fixed = dset.attrs["n_fixed"]

        else:
            n_fixed = 0

            warnings.warn(
                "The 'fixed_param' group is not found in "
                f"the results of {tag}. Probably the "
                "results were obtained with an older "
                "version of the package. Please rerun "
                "FitModel to update the results. Setting "
                "the number of fixed parameters to zero."
            )

        print("\nModel parameters:")
        for param_idx in range(n_param):
            param_item = dset.attrs[f"parameter{param_idx}"]
            print(f"   - {param_item}")

        if n_fixed == 0:
            print("\nFixed parameters: none")

        else:
            dset = hdf5_file[f"results/fit/{tag}/fixed_param"]

            print("\nFixed parameters:")
            for param_item in hdf5_file[f"results/fit/{tag}/fixed_param"]:
                print(f"   - {param_item}")

    print(f"\nInclude photometry: {inc_phot}")
    print(f"Include spectra: {inc_spec}")

    res_phot = None
    res_spec = None

    # Photometry residuals

    if inc_phot and objectbox.filters is not None:
        if isinstance(inc_phot, bool) and inc_phot:
            inc_phot = objectbox.filters

        model_phot = multi_photometry(
            datatype="model",
            spectrum=spectrum,
            filters=inc_phot,
            parameters=parameters,
            radtrans=radtrans,
            verbose=False,
        )

        res_phot = {}

        for item in inc_phot:
            transmission = ReadFilter(item)
            res_phot[item] = np.zeros(objectbox.flux[item].shape)

            if objectbox.flux[item].ndim == 1:
                res_phot[item][0] = transmission.mean_wavelength()
                res_phot[item][1] = (
                    objectbox.flux[item][0] - model_phot.flux[item]
                ) / objectbox.flux[item][1]

            elif objectbox.flux[item].ndim == 2:
                for j in range(objectbox.flux[item].shape[1]):
                    res_phot[item][0, j] = transmission.mean_wavelength()
                    res_phot[item][1, j] = (
                        objectbox.flux[item][0, j] - model_phot.flux[item]
                    ) / objectbox.flux[item][1, j]

    # Spectra residuals

    if inc_spec and objectbox.spectrum is not None:
        res_spec = {}

        if spectrum == "petitradtrans":
            # Calculate the petitRADTRANS spectrum only once
            # Smoothing and resampling not with get_model
            model = radtrans.get_model(parameters)

        for key in objectbox.spectrum:
            if isinstance(inc_spec, bool) or key in inc_spec:
                wavel_range = (
                    0.9 * objectbox.spectrum[key][0][0, 0],
                    1.1 * objectbox.spectrum[key][0][-1, 0],
                )

                wl_new = objectbox.spectrum[key][0][:, 0]
                spec_res = objectbox.spectrum[key][3]

                if spectrum == "planck":
                    readmodel = ReadPlanck(wavel_range=wavel_range)

                    model = readmodel.get_spectrum(
                        model_param=parameters, spec_res=1000.0
                    )

                    # Separate resampling to the new wavelength points

                    flux_new = spectres.spectres(
                        wl_new,
                        model.wavelength,
                        model.flux,
                        spec_errs=None,
                        fill=0.0,
                        verbose=True,
                    )

                elif spectrum == "petitradtrans":
                    # Smoothing to the instrument resolution
                    flux_smooth = convolve_spectrum(
                        model.wavelength, model.flux, spec_res
                    )

                    # Resampling to the new wavelength points
                    flux_new = spectres.spectres(
                        wl_new,
                        model.wavelength,
                        flux_smooth,
                        spec_errs=None,
                        fill=0.0,
                        verbose=True,
                    )

                else:
                    # Resampling to the new wavelength points
                    # is done by the get_model method

                    readmodel = ReadModel(spectrum, wavel_range=wavel_range)

                    if "teff_0" in parameters and "teff_1" in parameters:
                        # Binary system

                        param_0 = binary_to_single(parameters, 0)

                        model_spec_0 = readmodel.get_model(
                            param_0,
                            spec_res=spec_res,
                            wavel_resample=wl_new,
                        )

                        param_1 = binary_to_single(parameters, 1)

                        model_spec_1 = readmodel.get_model(
                            param_1,
                            spec_res=spec_res,
                            wavel_resample=wl_new,
                        )

                        # Weighted flux of two spectra for atmospheric asymmetries
                        # Or simply the same in case of an actual binary system

                        if "spec_weight" in parameters:
                            flux_comb = (
                                parameters["spec_weight"] * model_spec_0.flux
                                + (1.0 - parameters["spec_weight"]) * model_spec_1.flux
                            )

                        else:
                            flux_comb = model_spec_0.flux + model_spec_1.flux

                        model_spec = create_box(
                            boxtype="model",
                            model=spectrum,
                            wavelength=wl_new,
                            flux=flux_comb,
                            parameters=parameters,
                            quantity="flux",
                        )

                    else:
                        # Single object

                        model_spec = readmodel.get_model(
                            parameters,
                            spec_res=spec_res,
                            wavel_resample=wl_new,
                        )

                    flux_new = model_spec.flux

                data_spec = objectbox.spectrum[key][0]
                res_tmp = (data_spec[:, 1] - flux_new) / data_spec[:, 2]

                res_spec[key] = np.column_stack([wl_new, res_tmp])

    print("\nResiduals (sigma):")

    if res_phot is not None:
        for item in inc_phot:
            if res_phot[item].ndim == 1:
                print(f"   - {item} = {res_phot[item][1]:.2f}")

            elif res_phot[item].ndim == 2:
                for j in range(res_phot[item].shape[1]):
                    print(f"   - {item} = {res_phot[item][1, j]:.2f}")

    if res_spec is not None:
        for key in objectbox.spectrum:
            if isinstance(inc_spec, bool) or key in inc_spec:
                print(
                    f"   - {key}: min = {np.nanmin(res_spec[key]):.2f}, "
                    f"max = {np.nanmax(res_spec[key]):.2f}"
                )

    chi2_stat = 0
    n_dof = 0

    if res_phot is not None:
        for key, value in res_phot.items():
            if value.ndim == 1:
                chi2_stat += value[1] ** 2
                n_dof += 1

            elif value.ndim == 2:
                for i in range(value.shape[1]):
                    chi2_stat += value[1][i] ** 2
                    n_dof += 1

    if res_spec is not None:
        for key, value in res_spec.items():
            count_nan = np.sum(np.isnan(value[:, 1]))
            chi2_stat += np.nansum(value[:, 1] ** 2)

            n_dof += value.shape[0]
            n_dof -= count_nan

    print(f"\nNumber of data points = {n_dof}")
    print(f"Number of model parameters = {n_param}")
    print(f"Number of fixed parameters = {n_fixed}")

    n_dof -= n_param - n_fixed

    print(f"Number of degrees of freedom = {n_dof}")

    chi2_red = chi2_stat / n_dof

    print(f"\nchi2 = {chi2_stat:.2f}")
    print(f"reduced chi2 = {chi2_red:.2f}")

    return create_box(
        boxtype="residuals",
        name=objectbox.name,
        photometry=res_phot,
        spectrum=res_spec,
        chi2=chi2_stat,
        n_dof=n_dof,
        chi2_red=chi2_red,
    )
