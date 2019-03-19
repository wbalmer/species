"""
Utility functions for photometry.
"""

import spectres
import numpy as np

from species.core import box
from species.read import read_model, read_calibration, read_filter


def multi_photometry(datatype,
                     spectrum,
                     filters,
                     model_par):
    """
    Parameters
    ----------
    datatype : str
        Data type ('model' or 'calibration').
    spectrum : str
        Spectrum name (e.g., 'drift-phoenix').
    filters : tuple(str, )
        Filter IDs.
    model_par : dict
        Model parameter values.

    Returns
    -------
    species.core.box.SynphotBox
        Box with synthetic photometry.
    """

    flux = {}

    if datatype == 'model':
        for item in filters:
            readmodel = read_model.ReadModel(spectrum, item)
            flux[item] = readmodel.get_photometry(model_par)

    elif datatype == 'calibration':
        for item in filters:
            readcalib = read_calibration.ReadCalibration(spectrum, item)
            flux[item] = readcalib.get_photometry(model_par)

    return box.create_box('synphot', name='synphot', flux=flux)


def apparent_to_absolute(app_mag,
                         distance):
    """
    Parameters
    ----------
    app_mag : float or numpy.ndarray
        Apparent magnitude (mag).
    distance : float or numpy.ndarray
        Distance (pc).

    Returns
    -------
    float or numpy.ndarray
        Absolute magnitude (mag).
    """

    return app_mag - 5.*np.log10(distance) + 5.

def get_residuals(model,
                  model_par,
                  filters,
                  objectbox):
    """
    Parameters
    ----------
    model : str
        Atmospheric model.
    model_par : dict
        Model parameters and values.
    filters : tuple(str, )
        Filter IDs. All available photometry of the object is used if set to None.
    objectbox : species.core.box.ObjectBox
        Box with the photometry and/or spectrum of an object.

    Returns
    -------
    species.core.box.ResidualsBox
        Box with the photometry and/or spectrum residuals.
    """

    if filters is None:
        filters = objectbox.filter

    model_phot = multi_photometry(datatype='model',
                                  spectrum=model,
                                  filters=filters,
                                  model_par=model_par)

    if objectbox.flux is not None:
        res_phot = np.zeros((2, len(objectbox.flux)))

        for i, item in enumerate(objectbox.flux):
            transmission = read_filter.ReadFilter(item)

            res_phot[0, i] = transmission.mean_wavelength()
            res_phot[1, i] = (objectbox.flux[item][0]-model_phot.flux[item])/objectbox.flux[item][1]

    else:
        res_phot = None

    if objectbox.spectrum is not None:
        wl_range = (0.9*objectbox.spectrum[0, 0], 1.1*objectbox.spectrum[-1, 0])

        readmodel = read_model.ReadModel(model, wl_range)
        model = readmodel.get_model(model_par)

        wl_new = objectbox.spectrum[:, 0]

        flux_new = spectres.spectres(new_spec_wavs=wl_new,
                                     old_spec_wavs=model.wavelength,
                                     spec_fluxes=model.flux,
                                     spec_errs=None)

        res_spec = np.zeros((2, objectbox.spectrum.shape[0]))

        res_spec[0, :] = wl_new
        res_spec[1, :] = (objectbox.spectrum[:, 1]-flux_new)/objectbox.spectrum[:, 2]

    else:
        res_spec = None

    return box.create_box(boxtype='residuals',
                          name=objectbox.name,
                          photometry=res_phot,
                          spectrum=res_spec)