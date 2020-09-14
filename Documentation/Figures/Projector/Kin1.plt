set xrange [0.0:0.5]
plot "Finite_T/Kin.dat" u (1/$1):2:3 w e,  "Proj_tt1/Kin.dat" u (1/($1*2)):2:3 w e, "Proj_Kekule/Kin.dat" u (1/($1*2)):2:3 w e
