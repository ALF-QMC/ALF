set size 0.85,0.8
set terminal epslatex color standalone
set output "PAM.tex"
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
set key top right
set pointsize 1.5
#set ytics -2.95,0.01, -2.85
#set yrange [-109:-108.4]
set xrange [0.4:1.6] 
set ylabel "$ S^{ff}(\\pi,\\pi)  $"
set title " $ \\theta = 10,   U_f/t = 4.0 $ "
set xlabel "$ V/t $"

plot "Spin_PAM_L10.dat" u 1:2:3 w e lc rgb col8 lt 1 pt 5 t "L=10", \
     "Spin_PAM_L8.dat" u 1:2:3 w e lc rgb col7 lt 1 pt 7 t "L=8", \
     "Spin_PAM_L6.dat" u 1:2:3 w e lc rgb col6 lt 1 pt 9 t "L=6", \
     "Spin_PAM_L4.dat" u 1:2:3 w e lc rgb col5 lt 1 pt 11 t "L=4"
#plot "Finite_T/Kin.dat" u (1/$1):2:3 w e lc rgb col8 lt 1 pt 5 t "Finite T", \
#     "Proj_tt1/Kin.dat" u (1/($1*2 + 1.0)):2:3 w e lc rgb col7  lt 1 pt 7 t "t-t'", \
#     "Proj_Kekule/Kin.dat" u (1/($1*2 + 1.0)):2:3 w e lc rgb col6 lt 1 pt 9 t "Kekule"
unset multiplot
set output 
!pdflatex PAM.tex
!open PAM.pdf
