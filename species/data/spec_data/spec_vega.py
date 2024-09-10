"""
Module for adding a flux-calibrated spectrum of Vega to the database.
"""

from pathlib import Path

import h5py
import numpy as np
import pooch
import requests

from astropy.io import fits
from typeguard import typechecked


@typechecked
def add_vega(input_path: str, database: h5py._hl.files.File) -> None:
    """
    Function for adding a flux-calibrated spectrum of Vega to the
    database. The latest spectrum (alpha_lyr_stis_011.fits) is
    downloaded from the STScI archive (see `CALSPEC page <https://
    www.stsci.edu/hst/instrumentation/reference-data-for-calibration
    -and-tools/astronomical-catalogs/calspec>`_ for details).

    Parameters
    ----------
    input_path : str
        Path of the data folder.
    database : h5py._hl.files.File
        Database.

    Returns
    -------
    NoneType
        None
    """

    data_file = Path(input_path) / "alpha_lyr_stis_011.fits"

    if not data_file.exists():
        print()

        try:
            url = (
                "https://archive.stsci.edu/hlsps/reference-atlases"
                "/cdbs/current_calspec/alpha_lyr_stis_011.fits"
            )

            pooch.retrieve(
                url=url,
                known_hash="60aebf5c193223f69061cd176d6309730c3210051fffad0dd6ad44475199ceaa",
                fname="alpha_lyr_stis_011.fits",
                path=input_path,
                progressbar=True,
            )

        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError):
            url = (
                "https://home.strw.leidenuniv.nl/~stolker/"
                "species/alpha_lyr_stis_011.fits"
            )

            pooch.retrieve(
                url=url,
                known_hash="60aebf5c193223f69061cd176d6309730c3210051fffad0dd6ad44475199ceaa",
                fname="alpha_lyr_stis_011.fits",
                path=input_path,
                progressbar=True,
            )

    if "spectra/calibration/vega" in database:
        del database["spectra/calibration/vega"]

    vega_data = fits.getdata(data_file, ext=1)
    wavelength = vega_data["WAVELENGTH"]  # (Angstrom)
    flux = vega_data["FLUX"]  # (erg s-1 cm-2 A-1)
    error_stat = vega_data["STATERROR"]  # (erg s-1 cm-2 A-1)
    error_sys = vega_data["SYSERROR"]  # (erg s-1 cm-2 A-1)

    wavelength *= 1e-4  # (Angstrom) -> (um)
    flux *= 1e-3 * 1e4  # (erg s-1 cm-2 A-1) -> (W m-2 um-1)
    error_stat *= 1e-3 * 1e4  # (erg s-1 cm-2 A-1) -> (W m-2 um-1)
    error_sys *= 1e-3 * 1e4  # (erg s-1 cm-2 A-1) -> (W m-2 um-1)

    print("\nAdding spectrum: Vega", end="", flush=True)

    database.create_dataset(
        "spectra/calibration/vega", data=np.vstack((wavelength, flux, error_stat))
    )

    print("\nReference: Bohlin et al. 2014, PASP, 126")
    print("URL: https://ui.adsabs.harvard.edu/abs/2014PASP..126..711B/abstract")
