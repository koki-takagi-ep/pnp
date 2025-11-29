# Makefile for 1D Poisson-Nernst-Planck Solver
#
# Usage:
#   make          - Build the solver
#   make run      - Build and run with default parameters
#   make plot     - Generate plots from results
#   make clean    - Remove build artifacts
#   make all      - Build, run, and plot

# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic
INCLUDES = -I./include

# Directories
SRC_DIR = src
BUILD_DIR = build
RESULTS_DIR = results

# Source files
SOURCES = $(SRC_DIR)/pnp_solver.cpp $(SRC_DIR)/main.cpp
OBJECTS = $(BUILD_DIR)/pnp_solver.o $(BUILD_DIR)/main.o
TARGET = $(BUILD_DIR)/pnp_solver

# Default target
.PHONY: default
default: $(TARGET)

# Build target
$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^
	@echo "Build complete: $(TARGET)"

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

# Run the solver with default parameters
.PHONY: run
run: $(TARGET)
	@mkdir -p $(RESULTS_DIR)
	./$(TARGET) --output $(RESULTS_DIR)/pnp_results.dat

# Run parametric study
.PHONY: parametric
parametric: $(TARGET)
	@mkdir -p $(RESULTS_DIR)
	@echo "Running parametric study with different surface potentials..."
	./$(TARGET) --phi0 25  --output $(RESULTS_DIR)/pnp_phi25mV.dat
	./$(TARGET) --phi0 50  --output $(RESULTS_DIR)/pnp_phi50mV.dat
	./$(TARGET) --phi0 100 --output $(RESULTS_DIR)/pnp_phi100mV.dat
	./$(TARGET) --phi0 200 --output $(RESULTS_DIR)/pnp_phi200mV.dat
	@echo ""
	@echo "Running parametric study with different concentrations..."
	./$(TARGET) --c0 0.1 --output $(RESULTS_DIR)/pnp_c0_0.1M.dat
	./$(TARGET) --c0 1.0 --output $(RESULTS_DIR)/pnp_c0_1.0M.dat
	./$(TARGET) --c0 5.0 --output $(RESULTS_DIR)/pnp_c0_5.0M.dat
	@echo "Parametric study complete."

# Generate plots
.PHONY: plot
plot:
	@if [ -f $(RESULTS_DIR)/pnp_results.dat ]; then \
		python3 scripts/plot_results.py; \
	else \
		echo "No results found. Run 'make run' first."; \
	fi

# Generate all plots including parametric study
.PHONY: plot-all
plot-all:
	python3 scripts/plot_results.py
	python3 scripts/plot_parametric.py

# Build, run, and plot
.PHONY: all
all: $(TARGET) run plot

# Full pipeline with parametric study
.PHONY: full
full: $(TARGET) parametric plot-all

# Clean build artifacts
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(RESULTS_DIR)/*.dat

# Clean everything including results and plots
.PHONY: distclean
distclean: clean
	rm -rf $(RESULTS_DIR)

# Help
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  make          - Build the solver"
	@echo "  make run      - Build and run with default parameters"
	@echo "  make parametric - Run parametric study"
	@echo "  make plot     - Generate plots from results"
	@echo "  make plot-all - Generate all plots including parametric"
	@echo "  make all      - Build, run, and plot"
	@echo "  make full     - Full pipeline with parametric study"
	@echo "  make clean    - Remove build artifacts"
	@echo "  make distclean- Remove all generated files"
	@echo "  make help     - Show this help message"

# Dependencies
$(BUILD_DIR)/pnp_solver.o: $(SRC_DIR)/pnp_solver.cpp include/pnp_solver.hpp
$(BUILD_DIR)/main.o: $(SRC_DIR)/main.cpp include/pnp_solver.hpp
