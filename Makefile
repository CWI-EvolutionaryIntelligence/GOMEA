default: install

install-deps:
	pip3 -q install build --user
	pip3 -q install cython==3.0.10 --user
	pip3 -q install numpy==2.0.1 --user

install-meson-deps:
	pip3 -q install meson==1.5.1 --user

build-sdist: install-deps
	python3 -m build --sdist

build-wheel: install-deps
	python3 -m build --wheel

build-all: build-sdist build-wheel

dev: install-deps install-meson-deps
	meson setup build --python.install-env auto
	meson compile -C build
	python -m pip install --no-build-isolation --editable .

build-debug: install-deps install-meson-deps
	meson setup build-debug --python.install-env auto --buildtype="debug"
	meson compile -C build-debug

install-debug: install-deps install-meson-deps build-debug
	python -m pip install --no-binary :all: --debug --no-build-isolation --editable .

debug: uninstall build-debug
	mkdir -p debug/gomea/
	cp -r build-debug/* debug/gomea/
	cp -r gomea/__init__.py debug/gomea/

#cpp:
#	@mkdir -p build
#	g++ -g -Wall -std=c++17 -DCPP_STANDALONE -I./ -Igomea/ -IEigen/ gomea/src/discrete/*.cpp gomea/src/fitness/*.cpp gomea/src/common/*.cpp gomea/src/utils/*.cpp -o build/DiscreteGOMEA

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
CXXFLAGS=-O2 -fopenmp -Wall -std=c++17 -DCPP_STANDALONE $(CXX_INC)
#CXXFLAGS_O=-O2 -fopenmp -Wall -std=c++17 -DCPP_STANDALONE $(CXX_INC) #-pg ## For profiling
SRCDIR=gomea/src
OBJDIR_DISCRETE=build_cpp/obj_discrete
BINDIR=build_cpp

# list of all source files
SRCS_DISCRETE=$(wildcard $(SRCDIR)/fitness/*.cpp $(SRCDIR)/fitness/benchmarks-discrete/*.cpp $(SRCDIR)/common/*.cpp $(SRCDIR)/utils/*.cpp $(SRCDIR)/discrete/*.cpp)
SRCS_RV=$(wildcard $(SRCDIR)/fitness/*.cpp $(SRCDIR)/fitness/benchmarks-discrete/*.cpp $(SRCDIR)/common/*.cpp $(SRCDIR)/utils/*.cpp $(SRCDIR)/discrete/*.cpp)

# generate a list of object files based on source files
OBJS_DISCRETE=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR_DISCRETE)/%.o,$(SRCS_DISCRETE))

# the final executable file
TARGET_DISCRETE=$(BINDIR)/DiscreteGOMEA

#### Same, but for discrete version ####
$(OBJDIR_DISCRETE)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET_DISCRETE): $(OBJS_DISCRETE)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compile discrete version of code, in debug mode
cpp: $(TARGET_DISCRETE)
