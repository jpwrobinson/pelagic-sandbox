

import copernicusmarine
copernicusmarine.login(username='jrobinson12', password='ft3@etJ8X.Lxg4c')


copernicusmarine.subset(
    dataset_id="cmems_obs-oc_glo_bgc-plankton_my_l4-olci-300m_P1M",
    variables=["CHL"],
    minimum_longitude=50.0,   # west edge (Amirantes)
    maximum_longitude=57.0,   # east edge (inner islands)
    minimum_latitude=-9.0,    # south
    maximum_latitude=-2.0,    # north
    start_datetime="2016-04-01",
    end_datetime="2025-10-01",
    file_format="netcdf",
    output_directory="./seychelles_chl_300m"
)