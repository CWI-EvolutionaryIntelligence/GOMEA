default: reinstall

here: 
	python3 setup.py build_ext --inplace

build:
	pip3 -q install build --user
	pip3 -q install cython --user
	python3 -m build --sdist
	python3 -m build --wheel

src-install:
	pip3 -q install build --user
	pip3 -q install cython --user
	python3 -m build --sdist
	pip3 install dist/*.tar.gz --user

cpp:
	@mkdir -p build
	g++ -g -Wall -std=c++17 -DCPP_STANDALONE -I./ -Igomea/ -IEigen/ gomea/src/discrete/*.cpp gomea/src/fitness/*.cpp gomea/src/common/*.cpp gomea/src/utils/*.cpp -o build/DiscreteGOMEA

debug:
	python3 setup.py build_ext --inplace --debug

install:
	pip3 -q install build --user
	pip3 -q install cython --user
	python3 -m build --wheel
	pip3 install dist/*.whl --user

reinstall:
	python3 -m build --wheel
	pip3 install dist/*.whl --user --force-reinstall

doc:
	sphinx-build -b html docs/source/ docs/build/html

clean:
	rm -f *.so
	rm -rf gomea.egg-info/
	rm -rf build/
	rm -rf dist/
	rm -f gomea/*.cpp
	rm -f gomea/*.h
	rm -f gomea/*.so
	rm -rf cython_debug/
	rm -rf gomea/__pycache__/
