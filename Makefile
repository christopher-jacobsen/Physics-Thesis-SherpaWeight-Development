# makefile for SherpaWeight and SherpaME

BUILD_DIR = build

# ROOT configuration
ROOT_CFLAGS  = $(shell root-config --cflags)
ROOT_LDFLAGS = $(shell root-config --ldflags)
ROOT_LIBS    = $(shell root-config --libs)

# Sherpa configuration
SHERPA_CPPFLAGS = $(shell sherpa-config --cppflags)
SHERPA_LDFLAGS  = $(shell sherpa-config --ldflags)
SHERPA_LIBS     = $(shell sherpa-config --libs)

# compiler and linker setup
CXX = clang++
CPP_FLAGS = $(ROOT_CFLAGS) $(SHERPA_CPPFLAGS) -std=c++11 -ISource/Common
LD = $(CXX)
LD_FLAGS = $(ROOT_LDFLAGS) $(SHERPA_LDFLAGS) $(ROOT_LIBS) $(SHERPA_LIBS)

#CPP_FILES := $(wildcard src/*.cpp)
#OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))



SHERPA_WEIGHT_SOURCE = $(wildcard Source/SherpaWeight/*.cpp) $(wildcard Source/Common/*.cpp)
SHERPA_WEIGHT_DEPS   = $(SHERPA_WEIGHT_SOURCE) $(wildcard Source/SherpaWeight/*.h) $(wildcard Source/Common/*.h)

SHERPA_ME_SOURCE = $(wildcard Source/SherpaME/*.cpp) $(wildcard Source/Common/*.cpp)
SHERPA_ME_DEPS   = $(SHERPA_ME_SOURCE) $(wildcard Source/SherpaME/*.h) $(wildcard Source/Common/*.h)


all: $(BUILD_DIR)/SherpaWeight $(BUILD_DIR)/SherpaME

$(BUILD_DIR)/SherpaWeight: $(SHERPA_WEIGHT_DEPS)
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CPP_FLAGS) $(LD_FLAGS) $(SHERPA_WEIGHT_SOURCE) -o $@


$(BUILD_DIR)/SherpaME: $(SHERPA_ME_DEPS)
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CPP_FLAGS) $(LD_FLAGS) $(SHERPA_ME_SOURCE) -o $@

clean:
	rm -rf $(BUILD_DIR)
