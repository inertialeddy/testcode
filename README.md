# testcode

To link to gotm, insert correct path for 1) gotm libraries and 2) gotm modules in 'optfile'

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
