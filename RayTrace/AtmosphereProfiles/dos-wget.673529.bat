:: DOS batch script to download selected files from rda.ucar.edu using Wget
::
:: Experienced Wget Users: add additional command-line flags here
::   Use the -r (--recursive) option with care
set opts=-N
::
set cert_opt=
:: If you get a certificate verification error (version 1.10 or higher),
:: uncomment the following line:
::set cert_opt=--no-check-certificate
::
:: download the file(s)
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202202.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202201.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202201.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202201.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202201.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202202.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202202.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202202.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202203.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202203.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202203.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202204.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202203.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202204.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202204.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202204.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202205.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202205.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202205.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202205.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202206.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202206.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202206.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202207.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202207.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202206.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202207.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202207.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202208.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202208.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202208.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202208.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202209.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202209.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202209.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202209.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202210.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202210.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202210.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202210.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202211.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202211.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202212.00Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202212.18Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202212.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202212.06Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202211.12Z.grb2.nc
wget %cert_opt% %opts% https://rda.ucar.edu/dsrqst/BOUZA673529/pgbh.gdas.202211.18Z.grb2.nc
