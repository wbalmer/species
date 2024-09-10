"""
Module for adding AKARI L- and T-type dwarf spectra from
`Sorahana et al. (2013) to the database.
"""

import gzip
import os
import shutil
import urllib.request

import h5py
import numpy as np
import pandas as pd

from astropy.io import fits
from typeguard import typechecked

from species.util.data_util import extract_tarfile
from species.util.data_util import convert_units


@typechecked
def add_sorahana2013(input_path: str, database: h5py._hl.files.File) -> None:
    """
    Function for adding the AKARI spectra of L- and T-type
    dwarfs from `Sorahana et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJ...767...77S/abstract>`_ to the database.

    Parameters
    ----------
    input_path : str
        Path of the data folder.
    database : h5py._hl.files.File
        The HDF5 database.

    Returns
    -------
    NoneType
        None
    """

    print_text = "spectra of AKARI L-T type objects from Sorahana et al. 2013"

    data_url = None
    data_file = os.path.join(input_path, "sorahana_2013.tar")
    data_folder = os.path.join(input_path, "Sorahana_2013/")

    # if not os.path.isfile(data_file):
    #     print(f"File not uploaded yet, sorry", end="", flush=True)
    #     # urllib.request.urlretrieve(data_url, data_file)
    #     print(" [DONE]")

    # if os.path.exists(data_folder):
    #     shutil.rmtree(data_folder)

    # print(f"Unpacking {print_text} (300 KB)...", end="", flush=True)
    # extract_tarfile(data_file, data_folder)
    # print(" [DONE]")

    spec_dict = {}

    table1 = pd.read_csv(os.path.join(data_folder, 'sorahana13_table1.csv'))

    for idx in table1.index:
        row = table1.iloc[idx]
        print(row)
        name = row.id
        surv, coord = name.split(" ")
        if "+" in coord:
            dec = '+'
            coord1, coord2 = coord.split(dec)
            files = surv+'-'+coord1[0:5]+dec+coord2[0:4]+'.spc'
        elif "-" in coord:
            dec = '-'
            coord1, coord2 = coord.split(dec)
            files = surv+'-'+coord1[0:5]+dec+coord2[0:4]+'.spc'
        else:
            files = surv+'-'+coord+'.spc'
        sptype = row.spt
        plx,plx_err = row.plx.split("(")
        plx_err = plx_err.replace(")","")

        if "." in sptype:
            sptype = sptype[:4]
        else:
            sptype = sptype[:2]

        spec_dict[name] = {"name": name, "sptype": sptype, "files": files, "parallax":float(plx), "parallax_error":float(plx_err)}

    database.create_group("spectra/sorahana+2013")

    fits_folder =data_folder

    print_message = ""

    data_dict = {}

    for _, _, files in os.walk(fits_folder):
        for _, filename in enumerate(files):
            fname_split = filename.split("_")

            if '.spc' in filename:
                data = pd.read_csv(data_folder+filename, names=['w','f', 'fe1', 'fe2'], sep='\s+')
                data = data.to_numpy()
                print(data.shape)

                for name, value in spec_dict.items():
                    if filename in value["files"]:
                        print(filename)
                        wl = np.flip(data[:,0])
                        print(wl.shape)
                        fl_mjy = np.flip(data[:,1])
                        print(fl_mjy.shape)
                        fl_err_mjy = np.abs(np.flip(data[:,1]-np.median([data[:,2],data[:,3]])))
                        print(fl_err_mjy.shape)

                        flux_in = np.array([[wl, fl_mjy, fl_err_mjy]])

                        flux_out = convert_units(
                            flux_in, ("µm", "mJy"), convert_from=True
                        )[0].T

                        print(flux_out.shape)

                        data_dict[name] = flux_out

    for name, value in spec_dict.items():
        empty_message = len(print_message) * " "
        print(f"\r{empty_message}", end="")

        print_message = f"Adding spectra... {name}"
        print(f"\r{print_message}", end="")

        if name in data_dict:

            sp_data = data_dict[name]

            dset = database.create_dataset(f"spectra/sorahana+2013/{name}", data=sp_data)

            dset.attrs["name"] = str(name).encode()
            dset.attrs["sptype"] = str(value["sptype"]).encode()
            dset.attrs["parallax"] =value["parallax"]
            dset.attrs["parallax_error"] = value["parallax_error"]
        else:
            print('No file for this one, sorry')

    empty_message = len(print_message) * " "
    print(f"\r{empty_message}", end="")

    print_message = "Adding spectra... [DONE]"
    print(f"\r{print_message}")
