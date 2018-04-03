# TITAN - data quality control

Automatic quality control of in-situ observations with an emphasis on spatial controls.

TITAN is designed to test all the observations referring to the same observation time simulatenously.
Currently, the statitics of the individual station time series is not considered.

Available checks are (applied sequentially as in this list):

* Plausibility check

* Climatological check, predefined range for each month (optional)

* Buddy-check

* Isolated event test (STEVE) both over and under multiple thresholds (optional)

* Check against a deterministic first-guess field (optional)

* Check against an ensemble of first-guess fields (optional)

* Spatial Consistency Test (SCT)

* Check elevations against digital elevation model (optional)

* Isolated event test (STEVE) both over and under multiple thresholds, 2nd application (optional)

* Detect isolated observations

Possibility to have observation black-list and keep(-it-no-matter-what)-list.

Installation Instructions
-------------------------
Ensure the following R-libraries (and their dependencies) are installed:

   * argparser
   * sp
   * raster
   * rgdal
   * ncdf4 (optional, used only if additional geographical information are required)
   * tripack (optional, used by STEVE)
   * igraph (optional, used to post-process radar-derived precipitation)


Running the program
-------------------
To see program options, run:

```
   titan.R --help
```

Copyright and license
---------------------
Copyright (C) 2017 MET Norway. TITAN is licensed under [GPL
version 3](https://github.com/metno/TITAN/blob/master/LICENSE) or (at
your option) any later version.

Contact
-------
E-mail: cristianl@met.no

