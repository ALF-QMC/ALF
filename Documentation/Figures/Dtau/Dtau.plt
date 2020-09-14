set size 1.1,1.0
set terminal epslatex color standalone
set output "Dtau.tex"
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


set multiplot
set key at graph 0.22,0.27
#set key bottom left 
set pointsize 1.5
set xtics 0.0, 0.01, 0.04
#set yrange [-80.5:-79.5]
set xrange [0.0:0.2*0.2] 
set ylabel "$ \\langle \\hat{H}  \\rangle $"
set xlabel "$ (\\Delta  \\tau t)^{2} $"
set title "Square, $L=6,  U/t=4$, $\\beta t = 5 $ "
#set x2label " $ \\frac{1}{2 \\theta t + \\beta t} $ "
plot "NoSymNoCheck/Ener.dat" u ($1*$1):2:3 w e lc rgb col8 lt 1 pt 5  t "Sym N Check N", \
     "NoSymCheck/Ener.dat"   u ($1*$1):2:3 w e lc rgb col7 lt 1 pt 7  t "Sym N Check Y", \
     "SymNoCheck/Ener.dat"   u ($1*$1):2:3 w e lc rgb col6 lt 1 pt 9  t "Sym Y Check N", \
     "SymCheck/Ener.dat"     u ($1*$1):2:3 w e lc rgb col5 lt 1 pt 11 t "Sym Y Check Y"

unset multiplot

set output 
!pdflatex Dtau.tex
!open Dtau.pdf
