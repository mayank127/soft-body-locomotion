CC  = g++
LD	= ld
RM 	= rm
ECHO = /bin/echo
PRINTF	= printf
SED	= sed
CP = cp
MV = mv


PROJECT_ROOT=.
SRCDIR = $(PROJECT_ROOT)
OBJDIR = $(PROJECT_ROOT)/obj

INCLUDES = -I external/VegaFEM-v2.0/libraries/include/
LIBS =  -lcorotationalLinearFEM -lvolumetricMesh -lsparseMatrix -lpolarDecomposition -lminivector -lrenderVolumetricMesh -lopenGLHelper -lGL -lglut -L external/VegaFEM-v2.0/libraries/lib/
TARGET = soft-body

SRCS := $(wildcard $(SRCDIR)/*.cpp)
INCS := $(wildcard $(SRCDIR)/*.h)
OBJS := $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

.PHONY: all setup clean distclean vega_lib

all: setup $(TARGET)

setup: vega_lib volume_render
	@$(ECHO) "Setting up compilation.."
	@mkdir -p obj

$(TARGET): $(OBJS)
	@$(ECHO) "Building executable..."
	@$(CC) -o $@  -O3 $(OBJS) $(INCLUDES) $(LIBS)

-include $(OBJS:.o=.d)

vega_lib: external/VegaFEM-v2.0.zip
	@cd external; unzip -qo VegaFEM-v2.0.zip
	@cd external/VegaFEM-v2.0/; ./build

volume_render: vega_lib external/VegaFEM-v2.0/libraries/renderVolumetricMesh/renderVolumetricMesh.h
	@cd external/VegaFEM-v2.0/libraries/renderVolumetricMesh/; make

$(OBJS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@$(PRINTF) "Compiling $(notdir $<)\n"
	@$(CC) $(INCLUDES) $(LIBS) -c $< -o $@ -MD

clean:
	@$(ECHO) -n "Cleaning up..."
	@$(RM) -rf $(OBJDIR) *~ $(SRCDIR)/*~ $(TARGET)
	@$(ECHO) "Done"

distclean:
	@$(ECHO) -n "Cleaning up.."
	@$(RM) -rf $(OBJDIR) *~  $(TARGET) external/VegaFEM-v2.0/
	@$(ECHO) "Done"