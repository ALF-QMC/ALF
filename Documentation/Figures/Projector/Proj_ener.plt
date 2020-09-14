set size 0.8,0.8
set terminal epslatex color standalone
set output "Proj_ener.tex"
set key center top


#
#set style line 1 lt 1 lc rgb '#F7FBFF' # very light blue
#set style line 2 lt 1 lc rgb '#DEEBF7' # 
#set style line 3 lt 1 lc rgb '#C6DBEF' # 
#set style line 4 lt 1 lc rgb '#9ECAE1' # light blue
#set style line 5 lt 1 lc rgb '#6BAED6' # 
#set style line 6 lt 1 lc rgb '#4292C6' # medium blue
#set style line 7 lt 1 lc rgb '#2171B5' #
#set style line 8 lt 1 lc rgb '#084594' # dark blue

col1='#F7FBFF'
col2='#DEEBF7'
col3='#C6DBEF'
col4='#9ECAE1'
col5='#6BAED6'
col6='#4292C6'
col7='#2171B5'
col8='#084594'


#set key at graph 0.7,0.99
set multiplot
#set origin 0.2,0.5
#set key at graph 0.5,-0.5
set pointsize 1.5
#set ytics -2.95,0.01, -2.85
#set yrange [-80.5:-79.5]
set xrange [0.0:0.21] 
set ylabel "$ \\langle \\hat{H}  \\rangle $"
set title "$L=6,  U/t=2$ "
set xlabel " $ \\frac{1}{\\beta t},  \\frac{1}{2 \\theta t + \\beta t} $ "
#set x2label " $ \\frac{1}{2 \\theta t + \\beta t} $ "
plot "Finite_T/Ener.dat" u (1/$1):2:3 w e lc rgb col8 lt 1 pt 5 t "Finite T", \
     "Proj_tt1/Ener.dat" u (1/($1*2 + 1.0)):2:3 w e lc rgb col7 lt 1 pt 7 t "t-t'", \
     "Proj_Kekule/Ener.dat" u (1/($1*2 + 1.0)):2:3 w e lc rgb col6  lt 1 pt 9 t "Kekule"
unset multiplot

set output 
!pdflatex Proj_ener.tex
!open Proj_ener.pdf
