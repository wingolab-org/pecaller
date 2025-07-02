# Makefile for pecaller2 project

# Configuration
CC := gcc
CFLAGS := -std=c99 -g -O3 -Wall -Wextra -pedantic
CFLAGS_PROD := -std=c99 -g -O3 -Wall
LIBS := -lm -lz -lpthread
PREFIX ?= /usr/local

# Directories
SRC_DIR := c
BUILD_DIR := build
PERL_DIR := perl
BIN_DIR := $(PREFIX)/bin

# Find all C source files in c/ directory (not subdirectories)
C_SOURCES := $(wildcard $(SRC_DIR)/*.c)
# Extract program names from source files
PROGRAMS := $(patsubst $(SRC_DIR)/%.c,%,$(C_SOURCES))
# Build targets in build directory
TARGETS := $(patsubst %,$(BUILD_DIR)/%,$(PROGRAMS))

# Find all Perl scripts
PERL_SCRIPTS := $(wildcard $(PERL_DIR)/*.pl)

# Default target
.PHONY: all
all: $(TARGETS)

# Development build with extra warnings
.PHONY: dev
dev: CFLAGS := -std=c99 -g -O2 -Wall -Wextra -Wpedantic -Wshadow -Wcast-align -Wcast-qual -Wformat=2
dev: all

# Create build directory
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Build individual programs
$(BUILD_DIR)/%: $(SRC_DIR)/%.c | $(BUILD_DIR)
	@echo "Building $*..."
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

# Clean build artifacts
.PHONY: clean
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR)
	@find . -name "*.gc*" -exec rm -f {} +
	@rm -rf `find . -name "*.dSYM" -print`

# Install binaries and Perl scripts
.PHONY: install
install: all
	@echo "Installing to $(PREFIX)..."
	@install -d $(BIN_DIR)
	@install -m 755 $(TARGETS) $(BIN_DIR)/
	@install -m 755 $(PERL_SCRIPTS) $(BIN_DIR)/

# Uninstall
.PHONY: uninstall
uninstall:
	@echo "Uninstalling from $(PREFIX)..."
	@rm -f $(patsubst %,$(BIN_DIR)/%,$(PROGRAMS))
	@rm -f $(patsubst $(PERL_DIR)/%,$(BIN_DIR)/%,$(PERL_SCRIPTS))

# Show help
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  all        - Build all programs (default)"
	@echo "  dev        - Development build with extra warnings"
	@echo "  clean      - Remove build artifacts"
	@echo "  install    - Install binaries and Perl scripts to $(PREFIX)/bin"
	@echo "  uninstall  - Remove installed files"
	@echo "  help       - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  PREFIX     - Installation prefix (default: /usr/local)"
	@echo "  CC         - C compiler (default: gcc)"
	@echo ""
	@echo "Programs to be built:"
	@printf "  %s\n" $(PROGRAMS)

# Debug: show variables
.PHONY: debug
debug:
	@echo "C_SOURCES: $(C_SOURCES)"
	@echo "PROGRAMS: $(PROGRAMS)"
	@echo "TARGETS: $(TARGETS)"
	@echo "PERL_SCRIPTS: $(PERL_SCRIPTS)"
