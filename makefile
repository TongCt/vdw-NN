#Macro definition

FC = ifort
FFLAGS = -O2 -qopenmp
#end of Macro definition

OBJ_NN = main.o libax.o readv31.o readv21.o readv22.o readv23.o

main.exe: $(OBJ_NN)
	${FC} ${FFLAGS} $(OBJ_NN) -o main.exe 


clean:
	rm -f *.exe *.o 
