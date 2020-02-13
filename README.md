# TITAN - auTomatIc daTa quAlity coNtrol

Automatic quality control of in-situ observations with an emphasis on spatial controls.

TITAN is designed to test all the observations referring to the same observation time simultaneously.
Currently, the statistics of the individual station time series is not considered.

Available checks are (applied sequentially as in this list):

* Precipitation (in-situ) and temperature (field) cross-check (optional)

* Check elevations against digital elevation model (optional)

* Plausibility check

* Climatological check, predefined range for each month (optional)

* Buddy-check (event-based)

* Buddy-check

* Check against a deterministic first-guess field (optional)

* Check against an ensemble of first-guess fields (optional)

* Spatial Consistency Test (SCT)

* Check fOr hOLes in the field (COOL test) (optional)

* Detect isolated observations

Possibility to have observation black-list and keep(-it-no-matter-what)-list.

In case of precipitation, the program can adjust the values for the wind-induced loss.


Installation Instructions
-------------------------
Ensure the following R-libraries (and their dependencies) are installed:

   * argparser (tested with version 0.4. It does not work with version 0.6)
   * sp
   * raster
   * rgdal
   * ncdf4 (optional, used only if additional geographical information are required)
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

