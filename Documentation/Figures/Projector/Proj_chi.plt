set size 0.8,0.8
set terminal epslatex color standalone
set output "Proj_chi.tex"
set key center bottom
#set key at graph 0.7,0.99

col1='#F7FBFF'
col2='#DEEBF7'
col3='#C6DBEF'
col4='#9ECAE1'
col5='#6BAED6'
col6='#4292C6'
col7='#2171B5'
col8='#084594'

set multiplot
#set origin 0.2,0.5
#set key at graph 0.5,-0.5
set pointsize 1.5
#set ytics 0.0 0.5 0.1
set yrange [0.0:0.5]
set xrange [0.0:0.21] 

set ylabel "$  \\chi_{c}  $"
set title "$L=6,  U/t=2$ "
set xlabel " $ \\frac{1}{\\beta t} $ "
#set x2label " $ \\frac{1}{2 \\theta t + \\beta t} $ "
plot "Finite_T/Chi.dat" u (1/$1):($1*$2):($1*$3) w e lc rgb col8 lt 1 pt 5 t "Finite T"
unset multiplot
set output 
!pdflatex Proj_chi.tex
!open Proj_chi.pdf
