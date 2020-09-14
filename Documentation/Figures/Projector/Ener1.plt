set yrange [-80.5:-79.5]
set xrange [0.0:0.3]
plot "Finite_T/Ener.dat" u (1/$1):2:3 w e,  "Proj_tt1/Ener.dat" u (1/($1*2)):2:3 w e , "Proj_Kekule/Ener.dat" u (1/($1*2)):2:3 w e
