CFLAGS = -Wall

main: main.c
	@gcc-7 main.c -c ${CFLAGS} -lm -lgsl -lgslcblas
	@gcc-7 main.o ${CFLAGS} -lm -lgsl -lgslcblas -o main.out 

	@echo "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
	@echo "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
	@echo "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
	@echo "-_-_          main.c compiled.            -_-_"
	@echo "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
	@echo "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
	@echo "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"

# CFLAGS = -Wall -O2
# LIBS   = -lgsl -lm -lgslcblas
# CC     =gcc

# main: main.o
# 	$(CC) main.o ${CFLAGS} $(LIBS) -o main.out 
# main.o: main.c
# 	$(CC) main.c -c ${CFLAGS} $(LIBS)
# clean:
# 	rm *.{o,out}	