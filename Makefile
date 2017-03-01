Examples:
	(cd Libraries; make );\
	(cd Prog_8; make Examples )

Hub_Ising:
	(cd Libraries; make );\
	(cd Prog_8; make Hub_Ising )

Hub:
	(cd Libraries; make );\
	(cd Prog_8; make Hub )

Hub_Can:
	(cd Libraries; make );\
	(cd Prog_8; make Hub_Can )

SPT:
	(cd Libraries; make );\
	(cd Prog_8; make SPT )

Ising:
	(cd Libraries; make );\
	(cd Prog_8; make Ising )

Kondo_Honey:
	(cd Libraries; make );\
	(cd Prog_8; make Kondo_Honey )

clean: 	
	(cd Prog_8; make clean );\
	(cd Libraries; make clean )
