set terminal epslatex color standalone
set output "Langevin.tex"
set key center top
#set key at graph 0.7,0.99
set multiplot
#set origin 0.2,0.5
set size 0.9,0.9
#set key at graph 0.5,-0.5
set pointsize 1.5
set ytics -2.95,0.01, -2.85
set yrange [-2.95:-2.86] 
set ylabel "$ \\langle \\hat{H}  \\rangle $"
set title "$L=6,  \\beta  t=4$ "
set xlabel "\\Large $\\delta t_l $ "
f(x)  = a + b*x 
fit  [0.001:0.12] f(x)  "Langevin_ener.dat"  u 1:2:3 via a,b
plot "Langevin_ener.dat" u 1:2:3 w e lc rgb "black" lt 1 pt 7  lw 2 t "", f(x)  w l lt 1 lc rgb "blue"  t ""
unset multiplot
set output 
!pdflatex Langevin.tex
!open Langevin.pdf
