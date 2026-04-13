EXEC = main

OBJS = library_rel_ed.o module_for_slepc.o main.o

include ${SLEPC_DIR}/lib/slepc/conf/slepc_variables


${EXEC}: ${OBJS}
	-${FLINKER} -o ${EXEC} ${OBJS} ${SLEPC_LIB}
	${RM} ${OBJS} *.mod


main.o: module_for_slepc.o library_rel_ed.o
module_for_slepc.o: library_rel_ed.o


include ${SLEPC_DIR}/lib/slepc/conf/slepc_rules