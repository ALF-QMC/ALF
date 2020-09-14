set size 1.1,1.0
set terminal epslatex color standalone
set output "Dtau_1.tex"


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
#set key at graph 0.62,0.28
set key at graph 0.57,0.27
#set key  bottom left 
set logscale y
#set format y "$%.0t\\cdot10^{%+1T}$"
set format y "$10^{%+1T}$"
set pointsize 1.5
#set xtics 0.0, 0.01, 0.04
#set yrange [-80.5:-79.5]
set yrange [5*10e-8:1]
set xrange [0.0:20] 
set ylabel "$ |G(r=0,\\tau)| $"
set xlabel "$ \\tau $"
set title "Honeycomb, $L=6,  U/t=2$, $\\beta t = 40 $ "
#set x2label " $ \\frac{1}{2 \\theta t + \\beta t} $ "
plot "HoneycombFT_Beta40_Dtau0.2/Green_R0/g_R0" u 1:(abs($2)):3 w e lc rgb col8  lt 1 pt 5  t "Sym N Check Y, $\\Delta \\tau = 0.2 $" ,\
     "HoneycombFT_Beta40_Dtau0.1/Green_R0/g_R0" u 1:(abs($2)):3 w e lc rgb col7  lt 1 pt 7  t "Sym N Check Y, $\\Delta \\tau = 0.1 $", \
     "HoneycombTT_Beta40_Dtau0.2/Green_R0/g_R0" u 1:(abs($2)):3 w e lc rgb col6  lt 1 pt 9  t "Sym Y Check Y, $\\Delta \\tau = 0.2 $", \
     "HoneycombTT_Beta40_Dtau0.1/Green_R0/g_R0" u 1:(abs($2)):3 w e lc rgb col5  lt 1 pt 11 t "Sym Y Check Y, $\\Delta \\tau = 0.1 $"

unset multiplot

set output 
!pdflatex Dtau_1.tex
!open Dtau_1.pdf
