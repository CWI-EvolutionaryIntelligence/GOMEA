default: install-meson

here: 
	python3 setup.py build_ext --inplace

install-deps:
	pip3 -q install build --user
	pip3 -q install cython --user
	pip3 -q install numpy --user

build: install-deps
	python3 -m build --sdist
	python3 -m build --wheel

build-sdist: install-deps
	python3 -m build --sdist

build-wheel: install-deps
	python3 -m build --wheel

src-install: install-deps
	python3 -m build --sdist
	pip3 install dist/*.tar.gz --user

cpp:
	@mkdir -p build
	g++ -g -Wall -std=c++17 -DCPP_STANDALONE -I./ -Igomea/ -IEigen/ gomea/src/discrete/*.cpp gomea/src/fitness/*.cpp gomea/src/common/*.cpp gomea/src/utils/*.cpp -o build/DiscreteGOMEA

debug: install-deps
	python3 setup.py build_ext --inplace --debug

install: install-deps
	python3 -m build --wheel
	pip3 install dist/*.whl --user

meson-build:
	meson setup build --python.install-env auto
	meson compile -C build

meson-install: meson-build
	python -m pip install --no-build-isolation --editable .

docker-wheel: clean
	sudo pipx run cibuildwheel     

reinstall: install-deps
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
