# Notes
- PLASMA needs lapacke, note the 'e'. This is difficult to obtain on Gentoo Systems(unsure about the others)
  The debian/ubuntu package is called liblapacke-dev
  Hence I compiled lapack-3.7.0 myself and copied lpapcke and its headers to the respective locations in /usr/local

- tests fail on my system. make lib creates something.

- various pdfs for older Versions (I use Plasma-3.0.0 which is a port to OpenMP) can be found here: http://icl.cs.utk.edu/projectsfiles/plasma/pdf/

- Plasma has no pivoted QR decomposition... This means that we will use one without pivoting and hope that the improved stability they claim due to the tiling procedure affects us...

- Plasma uses lua scripts for the tuning parameters. To get rid of the Error message during calling plasma_init do this:

  set the environment variable PLASMA_TUNING_FILENAME to point to the default tuning file  (found  in tuning/default.lua).

  export PLASMA_TUNING_FILENAME=~/path_plasma/tuning/default.lua
  
  In my current setup its this: export PLASMA_TUNING_FILENAME=../plasma/tuning/default.lua