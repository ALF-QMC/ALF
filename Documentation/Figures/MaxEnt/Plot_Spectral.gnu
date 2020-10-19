#set terminal postscript color eps 25 enhanced
#set output "Spectral_F.eps"
set pm3d 
set title ""
set size 1.0,0.8
set logscale cb
set cbrange [0.01:10]
set xrange [0.0:23.0] 
set yrange [-8:8] 
#set ylabel " {/Symbol w}/t "
#set xlabel "q"
set tics out
#s= 2* 3.14159 *3/4
set pm3d corners2color c4
#set xtics ( '{/Symbol G}' 0, \
#            'K' 6, \
#            'M' 9, \
#            '{/Symbol G}'  18  )
plot "Green_Spectral.dat" u 1:($2):($4)  w image t ""

#!epstopdf Spectral_F.eps
#!open Spectral_F.pdf
