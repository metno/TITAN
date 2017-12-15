# TITAN - Temperature spatIal daTa quAlity coNtrol

Automatic quality control of in-situ temperature observations over a geographical region.

Available checks are:

* Plausibility check

* Buddy-check

* Spatial Consistency Test (SCT)

* detect isolated observations


Installation Instructions
-------------------------
Ensure the following R-libraries are installed:

   * argparser
   * sp
   * raster
   * rgdal


Running the program
-------------------
To see program options, run:

```
   titan.R --help
```

Input / Output
-------------------
Arguments:

* input_file. text file, colums separated by ";"
  * column names = lat,lon,elev,value (any order is accepted)
* output_file. text file, colums separated by ";" (same number of rows as input_file)
  * header: lat;lon;x;y;z;val;dqc;

Data quality control (dqc) codes:

- 0 = ok
- 1 = missing metadata
- 2 = plausibility test failed
- 3 = buddy check failed
- 4 = SCT failed
- 5 = isolated station 

Copyright and license
---------------------
Copyright (C) 2017 MET Norway. TITAN is licensed under [GPL
version 3](https://github.com/metno/TITAN/blob/master/LICENSE) or (at
your option) any later version.

Contact
-------
E-mail: cristianl@met.no

