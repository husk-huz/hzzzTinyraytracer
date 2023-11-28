SYSCONF_LINK = g++
CPPFLAGS     =
LDFLAGS      = 
LIBS         = -lm
CPP_FLAGS = '-std=c++11'

DESTDIR = ./
TARGET  = main

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: $(DESTDIR)$(TARGET)

$(DESTDIR)$(TARGET): $(OBJECTS)
	$(SYSCONF_LINK) ${CPP_FLAGS} -Wall $(LDFLAGS) -o $(DESTDIR)$(TARGET) $(OBJECTS) $(LIBS)

$(OBJECTS): %.o: %.cpp
	$(SYSCONF_LINK) ${CPP_FLAGS} -Wall $(CPPFLAGS) -c $(CFLAGS) $< -o $@

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)
	-rm -f *.tga

