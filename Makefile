# Target executable 
TARGET := CCC 

# Compiler
CXX := mpicxx 

# HDF5 library (use by default)
use_hdf := 1

# Folders
OBJ_DIR := obj
BIN_DIR := bin
SRC_DIR := src
INC_DIR := include
INCLUDE := -I$(INC_DIR)

# Flags and libraries
CXXFLAGS := -std=c++11 -O3
LDFLAGS  := -lHYPRE 

ifneq ($(use_hdf), 0)
  CXXFLAGS += -DUSE_HDF
  LDFLAGS  += -lhdf5
endif

# Files
SRC := $(wildcard $(SRC_DIR)/*.cc)
OBJ := $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(SRC))

all: build $(BIN_DIR)/$(TARGET)

$(BIN_DIR)/$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

build:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)

clean: 
	rm -rf $(OBJ_DIR) $(BIN_DIR)

cleaner:
	rm $(SRC_DIR)/*~
	rm $(INC_DIR)/*~
