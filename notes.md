# Notes
- PLASMA needs lapacke, note the 'e'. This is difficult to obtain on Gentoo Systems(unsure about the others)
  Hence I compiled lapack-3.7.0 myself and copied lpapcke and its headers to the respective locations in /usr/local

- tests fail on my system. make lib creates something.

- various pdfs for older Versions (I use Plasma-3.0.0 which is a port to OpenMP) can be found here: http://icl.cs.utk.edu/projectsfiles/plasma/pdf/

- Plasma has no pivoted QR decomposition... This means that we will use one without pivoting and hope that the improved stability they claim due to the tiling procedure affects us...