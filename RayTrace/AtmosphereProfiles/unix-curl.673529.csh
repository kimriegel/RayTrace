#! /bin/csh -f
#
# c-shell script to download selected files from rda.ucar.edu using curl
# NOTE: if you want to run under a different shell, make sure you change
#       the 'set' commands according to your shell's syntax
# after you save the file, don't forget to make it executable
#   i.e. - "chmod 755 <name_of_script>"
#
# you can add cURL options here (progress bars, etc.)
set opts = ""
#
# download the file(s)
# NOTE:  if you get 403 Forbidden errors when downloading the data files, check
#        the contents of the file 'auth_status.rda.ucar.edu'
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202202.00Z.grb2.nc -o pgbh.gdas.202202.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202201.00Z.grb2.nc -o pgbh.gdas.202201.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202201.18Z.grb2.nc -o pgbh.gdas.202201.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202201.06Z.grb2.nc -o pgbh.gdas.202201.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202201.12Z.grb2.nc -o pgbh.gdas.202201.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202202.06Z.grb2.nc -o pgbh.gdas.202202.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202202.18Z.grb2.nc -o pgbh.gdas.202202.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202202.12Z.grb2.nc -o pgbh.gdas.202202.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202203.06Z.grb2.nc -o pgbh.gdas.202203.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202203.00Z.grb2.nc -o pgbh.gdas.202203.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202203.12Z.grb2.nc -o pgbh.gdas.202203.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202204.12Z.grb2.nc -o pgbh.gdas.202204.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202203.18Z.grb2.nc -o pgbh.gdas.202203.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202204.00Z.grb2.nc -o pgbh.gdas.202204.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202204.06Z.grb2.nc -o pgbh.gdas.202204.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202204.18Z.grb2.nc -o pgbh.gdas.202204.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202205.00Z.grb2.nc -o pgbh.gdas.202205.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202205.06Z.grb2.nc -o pgbh.gdas.202205.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202205.12Z.grb2.nc -o pgbh.gdas.202205.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202205.18Z.grb2.nc -o pgbh.gdas.202205.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202206.00Z.grb2.nc -o pgbh.gdas.202206.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202206.12Z.grb2.nc -o pgbh.gdas.202206.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202206.06Z.grb2.nc -o pgbh.gdas.202206.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202207.06Z.grb2.nc -o pgbh.gdas.202207.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202207.00Z.grb2.nc -o pgbh.gdas.202207.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202206.18Z.grb2.nc -o pgbh.gdas.202206.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202207.18Z.grb2.nc -o pgbh.gdas.202207.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202207.12Z.grb2.nc -o pgbh.gdas.202207.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202208.00Z.grb2.nc -o pgbh.gdas.202208.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202208.06Z.grb2.nc -o pgbh.gdas.202208.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202208.12Z.grb2.nc -o pgbh.gdas.202208.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202208.18Z.grb2.nc -o pgbh.gdas.202208.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202209.00Z.grb2.nc -o pgbh.gdas.202209.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202209.06Z.grb2.nc -o pgbh.gdas.202209.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202209.18Z.grb2.nc -o pgbh.gdas.202209.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202209.12Z.grb2.nc -o pgbh.gdas.202209.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202210.00Z.grb2.nc -o pgbh.gdas.202210.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202210.06Z.grb2.nc -o pgbh.gdas.202210.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202210.12Z.grb2.nc -o pgbh.gdas.202210.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202210.18Z.grb2.nc -o pgbh.gdas.202210.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202211.00Z.grb2.nc -o pgbh.gdas.202211.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202211.06Z.grb2.nc -o pgbh.gdas.202211.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202212.00Z.grb2.nc -o pgbh.gdas.202212.00Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202212.18Z.grb2.nc -o pgbh.gdas.202212.18Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202212.12Z.grb2.nc -o pgbh.gdas.202212.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202212.06Z.grb2.nc -o pgbh.gdas.202212.06Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202211.12Z.grb2.nc -o pgbh.gdas.202211.12Z.grb2.nc
curl  $opts https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202211.18Z.grb2.nc -o pgbh.gdas.202211.18Z.grb2.nc
