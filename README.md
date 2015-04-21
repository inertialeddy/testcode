# testcode

To link to gotm, insert correct path for 1) gotm libraries and 2) gotm modules in 'optfile' by changing "path-to-gotm-libraries" and "path-to-gotm-modules" to the directories where you place the libraries and modules.

The libraries (in gotm_ifort_libraries) and modules (in gotm_ifort_modules) were created on a Mac with ifort. If you want to compile
and create your own GOTM libraries and modules, follow these steps:

cd gotm_allfiles

cd src

make distclean

make 

This should create the libraries and modules. Check whether these really exist using:

ls lib/IFORT/

ls modules/IFORT/ 

You should see a bunch of files in both these directories. If you want to be doubly sure these were not lying around from before the compilation, 
delete these libraries (and modules) and go through the compilation steps once again (starting from "cd gotm_allfiles").

-------------------------------------PRECIPITATION AND EVAPORATION FLUXES--------------------------------

These are set in fluxes.f90 using the variables 'swr', 'qloss', 'qlatent' and 'rainrate', defined as follows:

swr=  shortwave (W/m^2)

qloss = qnet - swr (W/m^2)

qlatent = Latent heat flux (W/m^2)

rainrate = Precipitation (mm/hr) 

The latent heat flux and rainrate are converted into evaporative and precipitative fluxes, respectively, in the file mixing_vertical.f90. Inside mixing_vertical.f90 the variable 'evap_flux' holds the evaporative flux and 'precip_flux' holds
the precipitation flux. 
