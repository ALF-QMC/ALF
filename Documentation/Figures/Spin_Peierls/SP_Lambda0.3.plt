set size 0.8,0.8
set terminal epslatex color standalone
set output "SP_Lambda0.3.tex"
set key top right


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
pi  =  3.1415927
set yrange [0.0:4.5]
set  ytics  0.0,1,4.5
set xrange [0.0:2.0*pi] 
set ylabel "$ S(q) $"
set title "$L=16,  U/t=1, \\omega_0 = 0.25, \\beta  = 8, \\lambda=0.3$ "
set xlabel " $q $ "
#set x2label " $ \\frac{1}{2 \\theta t + \\beta t} $ "
plot "Lx16_Ly1_Jx1.0_Jy0.000001_Om0.25_Lam0.3_Beta8/SpinZ_eqJK" u 1:2:3 w e lc rgb col8 lt 1 pt 5 t "ALF", \
     "Lx16_Ly1_Jx1.0_Jy0.000001_Om0.25_Lam0.3_Beta8/SpinZ_eqJK" u ($1+2*pi):2:3 w e lc rgb col8 lt 1 pt 5 t "", \
     "Lx16_Ly1_Jx1.0_Jy0.000001_Om0.25_Lam0.3_Beta8/spinztotJ_lambda0.3" u ($1):($2*4):($3*4) w e lc rgb col6 lt 1 pt 9 t "SSE"
unset multiplot

set output 
!pdflatex SP_Lambda0.3.tex
!open SP_Lambda0.3.pdf
