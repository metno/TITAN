# TITAN - auTomatIc daTa quAlity coNtrol

Automatic quality control of in-situ observations with an emphasis on spatial controls.

TITAN is designed to test all the observations referring to the same observation time simultaneously.
The statistics of the individual station time series is not considered.

Installation Instructions
-------------------------
Ensure the following R-libraries (and their dependencies) are installed:

   * argparser (tested with version 0.4. It does not work with version 0.6)
   * sp
   * raster
   * rgdal
   * ncdf4

Optional, only used for some functions:
   * igraph (used to post-process radar-derived precipitation)
   * RANN (used with background files)

The tests are based on [titanlib](https://github.com/metno/titanlib), which is required. 
The path to the files `titanlib.R` and `titanlib.so` must be specified in the configuration files as `titanlib_path` (see the examples in the test directory).


Running the program
-------------------

The environmental variable `TITANR_PATH` must be specified. This is the path to the TITAN `functions` directory (see the examples in the test directory).

A detailed description of how to run the program, with specific examples for most of the tests, is provided in the Wiki [here](https://github.com/metno/TITAN/wiki). 

Copyright and license
---------------------
Copyright (C) 2021 MET Norway. TITAN is licensed under [GPL version 3](https://github.com/metno/TITAN/blob/master/LICENSE) or (at your option) any later version.

Contact
-------
E-mail: cristianl@met.no

