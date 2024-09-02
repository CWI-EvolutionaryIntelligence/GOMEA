CXXL=$(CXX)
ifneq ($(OS),Windows_NT)
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		CXXL=clang++
	endif
endif

default: install

install-deps:
	pip3 -q install build --user
	pip3 -q install cython==3.0.10 --user
	pip3 -q install numpy==2.0.1 --user

install-meson-deps:
	pip3 -q install meson==1.5.1 --user

build-sdist: install-deps
	CXX=$(CXXL) python3 -m build --sdist

build-wheel: install-deps
	CXX=$(CXXL) python3 -m build --wheel

build-all: build-sdist build-wheel

dev: install-deps install-meson-deps
	CXX=$(CXXL) meson setup build --python.install-env auto
	CXX=$(CXXL) meson compile -C build
	python -m pip install --no-build-isolation --editable .

build-debug: install-deps install-meson-deps
	CXX=$(CXXL) meson setup build-debug --python.install-env auto --buildtype="debug"
	CXX=$(CXXL) meson compile -C build-debug

install-debug: install-deps install-meson-deps build-debug
	CXX=$(CXXL) python -m pip install --no-binary :all: --debug --no-build-isolation --editable .

debug: uninstall build-debug
	mkdir -p debug/gomea/
	cp -r build-debug/* debug/gomea/
	cp -r gomea/__init__.py debug/gomea/

install: install-deps build-wheel
	pip3 install dist/*.whl --user

reinstall: install-deps build-wheel
	pip3 install dist/*.whl --user --force-reinstall

src-install: install-deps build-sdist
	pip3 install dist/*.tar.gz --user

uninstall:
	pip3 uninstall -y gomea

cibuildwheel: clean
	sudo pipx run cibuildwheel     

doc:
	sphinx-build -b html docs/source/ docs/build/html

clean:
	rm -f *.so
	rm -rf gomea.egg-info/
	rm -rf build/
	rm -rf build_cpp/
	rm -rf dist/
	rm -rf build-debug/
	rm -rf debug/gomea/
	rm -f gomea/*.cpp
	rm -f gomea/*.h
	rm -f gomea/*.so
	rm -rf cython_debug/
	rm -rf gomea/__pycache__/


######################################################################################################


CXX=g++
CXX_INC=-I./ -Igomea/ -Igomea/lib/Eigen/ -Igomea/lib/cxxopts-3.1.1/include/
CXXFLAGS=-g -fopenmp -Wall -std=c++17 -DCPP_STANDALONE $(CXX_INC)  #-pg ## For profiling
SRCDIR=gomea/src
OBJDIR=build_cpp/obj
OBJDIR_DISCRETE=build_cpp/obj_discrete
OBJDIR_RV=build_cpp/obj_realvalued
BINDIR=build_cpp

# list of all source files
SRCS=$(wildcard $(SRCDIR)/fitness/*.cpp $(SRCDIR)/common/*.cpp $(SRCDIR)/utils/*.cpp)
SRCS_DISCRETE=$(wildcard $(SRCDIR)/fitness/benchmarks-discrete/*.cpp $(SRCDIR)/discrete/*.cpp)
SRCS_RV=$(wildcard $(SRCDIR)/fitness/benchmarks-rv/*.cpp $(SRCDIR)/real_valued/*.cpp)

# generate a list of object files based on source files
OBJS=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
OBJS_DISCRETE=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR_DISCRETE)/%.o,$(SRCS_DISCRETE))
OBJS_RV=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR_RV)/%.o,$(SRCS_RV))

# the final executable file
TARGET_DISCRETE=$(BINDIR)/DiscreteGOMEA
TARGET_RV=$(BINDIR)/RealValuedGOMEA

#### Object files ####
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR_DISCRETE)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR_RV)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET_DISCRETE): $(OBJS_DISCRETE) $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TARGET_RV): $(OBJS_RV) $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

cpp: $(TARGET_DISCRETE) $(TARGET_RV)
