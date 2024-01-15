#set terminal postscript color eps 25 enhanced
#set output "R_Plot_Name.eps"
set pm3d 
set title "R_Plot_Title"
set size 1.0,0.8
set logscale cb
set cbrange [0.01:10]
set xrange [0.0:R_xrange] 
set yrange [R_Om_st:R_Om_en] 
#set ylabel " {/Symbol w}/t "
#set xlabel "q"
set tics out
#s= 2* 3.14159 *3/4
set pm3d corners2color c4
#set xtics ( '{/Symbol G}' 0, \
#            'K' 6, \
#            'M' 9, \
#            '{/Symbol G}'  18  )
plot "R_Plot_File" u 1:($2):($4)  w image t ""

#!epstopdf R_Plot_Name.eps
#!open R_Plot_Name.pdf
