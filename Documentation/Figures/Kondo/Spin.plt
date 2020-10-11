set size 0.85,0.8
set terminal epslatex color standalone
set output "Spin.tex"
set key center top
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
#set ytics -2.95,0.01, -2.85
#set yrange [-109:-108.4]
set xrange [0.0:3.0] 
#set logscale y
set ylabel "$ S_{ff}(\\pi,\\pi)$"
set title "$L=4, \\beta t = 5, J_k/t = 2 $ "
set xlabel "$ U_f/t $"

plot "Constraint.dat" u 1:4:5 w e lc rgb col8 lt 1 pt 5 t ""
#plot "Finite_T/Kin.dat" u (1/$1):2:3 w e lc rgb col8 lt 1 pt 5 t "Finite T", \
#     "Proj_tt1/Kin.dat" u (1/($1*2 + 1.0)):2:3 w e lc rgb col7  lt 1 pt 7 t "t-t'", \
#     "Proj_Kekule/Kin.dat" u (1/($1*2 + 1.0)):2:3 w e lc rgb col6 lt 1 pt 9 t "Kekule"
unset multiplot
set output 
!pdflatex Spin.tex
!open Spin.pdf
