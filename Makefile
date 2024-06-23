CC := g++
CFLAGS := -std=c++17 -Werror -g
LDFLAGS := -lSDL2 -lfftw3
BIN_DIR := bin
SRC_DIR := src
OBJ_DIR := obj

# Source files
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC))

# Target
TARGET := $(BIN_DIR)/main

# Build
all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# Clean

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean