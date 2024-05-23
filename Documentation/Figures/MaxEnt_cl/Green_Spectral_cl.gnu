set terminal postscript color eps 25 enhanced
set output "Green_Spectral_cl.eps"
set pm3d 
set title "Green Classic"
set size 1.0,0.8
set logscale cb
set cbrange [0.01:10]
set xrange [0.0:23] 
set yrange [-8.0:8.0] 
#set ylabel " {/Symbol w}/t "
#set xlabel "q"
set tics out
#s= 2* 3.14159 *3/4
set pm3d corners2color c4
set xtics ( '0' 0, \
            '{/Symbol p}/2' 11.5, \
            '{/Symbol p}' 23 )
plot "Green_Spectral_cl.dat" u 1:($2):($4)  w image t ""

!epstopdf Green_Spectral_cl.eps
!open Green_Spectral_cl.pdf
