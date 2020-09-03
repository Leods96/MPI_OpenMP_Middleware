NAME = traffic
CC = mpicc
CXX = mpicxx
COMMONFLAGS = -Wall $(INCFLAGS) -fopenmp
CFLAGS = $(COMMONFLAGS)
CPPFLAGS = -std=c++11 $(COMMONFLAGS)
LDFLAGS =
SRCDIR = ./src
INCFLAGS = -I./include
BINDIR = /usr/local/bin
OBJDIR = ./build
CSOURCES = $(wildcard $(SRCDIR)/*.c)
CPPSOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(CSOURCES))
OBJECTS += $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(CPPSOURCES))

.PHONY = all clean install uninstall

all: $(ULFM_FILE) $(NAME)

$(NAME): $(OBJECTS)
	$(CXX) $^ -o $@ $(CPPFLAGS) $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $< -o $@ $(CPPFLAGS) $(LDFLAGS)

clean:
	rm -f $(NAME)
	rm -f $(OBJDIR)/*.o
install: all
	cp $(NAME) $(BINDIR)
uninstall:
	rm -f $(BINDIR)/$(NAME)