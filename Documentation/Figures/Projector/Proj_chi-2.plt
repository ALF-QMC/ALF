set size 0.8,0.8
set terminal epslatex color standalone
set output "Proj_chi-2.tex"
set key center bottom
#set key at graph 0.7,0.99
set multiplot
#set origin 0.2,0.5
#set key at graph 0.5,-0.5
set pointsize 1.5
#set ytics 0.0 0.5  0.1
set yrange [0.0:0.5]
set xrange [0.0:0.21] 

set ylabel "$  \\chi_{c}  $"
#set title "$L=6,  U/t=2$ "
set xlabel " $ \\frac{1}{\\beta t} $ "
#set x2label " $ \\frac{1}{2 \\theta t + \\beta t} $ "
plot "Finite_T/Chi.dat" u (1/$1):($1*$2):($1*$3) w e lc rgb "red" lt 1 pt 8 t "Finite T"
#plot "Langevin_chi.dat" u 1:2:3 w e lc rgb "black" lt 1 pt 7  lw 2 t "", f(x)  w l lt 1 lc rgb "blue"  t ""
unset multiplot
set output 
!pdflatex Proj_chi-2.tex
#!open Proj_chi-2.pdf
