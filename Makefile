CC       = g++
CFLAGS   = -Wall
LDFLAGS  = 
OBJFILES = linear_solver.o inv.o gbp.o
TARGET   = GABP

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~