set xrange [0.0:0.5]
plot "Finite_T/Pot.dat" u (1/$1):2:3 w e,  "Proj_tt1/Pot.dat" u (1/($1*2)):2:3 w e, "Proj_Kekule/Pot.dat" u (1/($1*2)):2:3 w e
