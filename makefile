CC = gcc
#FLAG = -O3 -Wall
#FLAG = -O3 -Wall -w
FLAG = -w
LIB = -lm
#CLLIB = -lblas -lf2c
#CLLIB = -llapack -lblas -lf2c

INCLUDE = def.h
SOURCE = main.c
OBJ = file.o err.o post.o theorem.o

run : $(OBJ) $(SOURCE) $(INCLUDE)
#	$(CC) $(FLAG) -o run $(SOURCE) $(OBJ) $(DEBUG) ${CLLIB} $(LIB)
	$(CC) $(FLAG) -o run $(SOURCE) $(OBJ) $(DEBUG) $(LIB)

file.o : file.c $(INCLUDE)
	$(CC) $(FLAG) -c file.c

err.o : err.c $(INCLUDE)
	$(CC) $(FLAG) -c err.c

post.o : post.c $(INCLUDE)
	$(CC) $(FLAG) -c post.c

theorem.o : theorem.c $(INCLUDE)
	$(CC) $(FLAG) -c theorem.c

clean:
	rm -f $(OBJ) *.o run* *.dat *.out
