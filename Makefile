TARGET = Test_dyn

SOURCES = \
	test/testDynamics.c \
	src/MatriceOps/MatriceOps.c \
	src/scaraDyn/scaraDynamics.c 

all: $(TARGET)

$(TARGET) : $(SOURCES)
	gcc -Iinclude -o $(TARGET) $(SOURCES)


