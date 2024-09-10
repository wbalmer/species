"""
Module for the photometric data and parallaxes from the
Database of Ultracool Parallaxes.
"""

from pathlib import Path

import h5py
import numpy as np
import pooch

from astropy.io import fits
from typeguard import typechecked

from species.util.data_util import update_sptype


@typechecked
def add_vlm_plx(input_path: str, database: h5py._hl.files.File) -> None:
    """
    Function for adding the Database of Ultracool Parallaxes to the
    database. The FITS file with data was originally downloaded from
    http://www.as.utexas.edu/~tdupuy/plx/Database_of_Ultracool_Parallaxes_files/vlm-plx-all.fits
    but that website has been taken offline, probably because of the
    new table at http://bit.ly/UltracoolSheet.

    Parameters
    ----------
    input_path : str
        Data folder.
    database : h5py._hl.files.File
        The HDF5 database that has been opened.

    Returns
    -------
    NoneType
        None
    """

    input_file = "vlm-plx-all.fits"
    data_file = Path(input_path) / "vlm-plx-all.fits"
    url = "https://home.strw.leidenuniv.nl/~stolker/species/vlm-plx-all.fits"

    if not data_file.exists():
        print()

        pooch.retrieve(
            url=url,
            known_hash="d31bb3162d7de890c09ebf9f0497d51159889b5f5e7c4da1ddf01f24d0c2b36f",
            fname=input_file,
            path=input_path,
            progressbar=True,
        )

    database.create_group("photometry/vlm-plx")

    with fits.open(data_file) as hdu_list:
        phot_data = hdu_list[1].data

    parallax = phot_data["PLX"]  # (mas)
    parallax_error = phot_data["EPLX"]  # (mas)

    name = phot_data["NAME"]
    name = np.core.defchararray.strip(name)

    sptype = phot_data["OSPTSTR"]
    sptype = np.core.defchararray.strip(sptype)

    sptype_nir = phot_data["ISPTSTR"]
    sptype_nir = np.core.defchararray.strip(sptype_nir)

    for i, item in enumerate(sptype):
        if item == "null":
            sptype[i] = sptype_nir[i]

    flag = phot_data["FLAG"]
    flag = np.core.defchararray.strip(flag)

    sptype = update_sptype(sptype)

    dtype = h5py.special_dtype(vlen=str)

    dset = database.create_dataset(
        "photometry/vlm-plx/name", (np.size(name),), dtype=dtype
    )
    dset[...] = name

    dset = database.create_dataset(
        "photometry/vlm-plx/sptype", (np.size(sptype),), dtype=dtype
    )
    dset[...] = sptype

    dset = database.create_dataset(
        "photometry/vlm-plx/flag", (np.size(flag),), dtype=dtype
    )
    dset[...] = flag

    database.create_dataset("photometry/vlm-plx/ra", data=phot_data["RA"])  # (deg)
    database.create_dataset("photometry/vlm-plx/dec", data=phot_data["DEC"])  # (deg)
    database.create_dataset("photometry/vlm-plx/parallax", data=parallax)
    database.create_dataset("photometry/vlm-plx/parallax_error", data=parallax_error)
    database.create_dataset("photometry/vlm-plx/Keck/NIRC.Y", data=phot_data["YMAG"])
    database.create_dataset("photometry/vlm-plx/MKO/NSFCam.J", data=phot_data["JMAG"])
    database.create_dataset("photometry/vlm-plx/MKO/NSFCam.H", data=phot_data["HMAG"])
    database.create_dataset("photometry/vlm-plx/MKO/NSFCam.K", data=phot_data["KMAG"])
    database.create_dataset("photometry/vlm-plx/MKO/NSFCam.Lp", data=phot_data["LMAG"])
    database.create_dataset("photometry/vlm-plx/MKO/NSFCam.Mp", data=phot_data["MMAG"])
    database.create_dataset("photometry/vlm-plx/2MASS/2MASS.J", data=phot_data["J2MAG"])
    database.create_dataset("photometry/vlm-plx/2MASS/2MASS.H", data=phot_data["H2MAG"])
    database.create_dataset(
        "photometry/vlm-plx/2MASS/2MASS.Ks", data=phot_data["K2MAG"]
    )

    database.close()
