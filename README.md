# TITAN - Temperature spatIal daTa quAlity coNtrol

command line
>R --vanilla file\_input file\_output < titan.R

Arguments:

- file\_input, text file with colums separated by ";" \
column names = lat,lon,elev,value (any order is accepted)
- file\_output, text file with colums separated by ";" \
same number of rows as file\_in \
header: lat;lon;x;y;z;val;dqc;\

data quality control (dqc) codes:

- 0 = ok
- 1 = missing metadata
- 2 = plausibility test failed
- 3 = buddy check failed
- 4 = SCT failed
- 5 = isolated station 

