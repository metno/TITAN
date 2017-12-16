# TITAN - Temperature spatIal daTa quAlity coNtrol

Automatic quality control of in-situ temperature observations over a geographical region.

Available checks are (applied sequentially as in this list):

* Plausibility check

* Buddy-check

* Spatial Consistency Test (SCT)

* check elevations against digital elevation model

* detect isolated observations


Installation Instructions
-------------------------
Ensure the following R-libraries are installed:

   * argparser
   * sp
   * raster
   * rgdal
   * ncdf4 (optional, used only if additional geographical information are required)


Running the program
-------------------
To see program options, run:

```
   titan.R --help
```

run a test case with:

```
./titan.R test/TA_2017072112.txt test/dqc_2017072112.txt -c -iS 3 -iB 3 -v
```

(titan.R input output convert_lat_lon_2_km iterate_SCT_3_times iterate_BuddyCheck_3_times verbose)

Copyright and license
---------------------
Copyright (C) 2017 MET Norway. TITAN is licensed under [GPL
version 3](https://github.com/metno/TITAN/blob/master/LICENSE) or (at
your option) any later version.

Contact
-------
E-mail: cristianl@met.no

