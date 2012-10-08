# Target to build
TARGET = final
EXECS = ./$(TARGET)

all: $(TARGET)

# Generic compile rules
.c.o: 
	gcc -c -O -Wall $<

# Generic compile and link
%: %.c 
	gcc -Wall -O3 -o ./$@ $^ 

clean:
	rm -f $(EXECS) *.o *.a