default: install

install-deps:
	pip3 -q install build --user
	pip3 -q install cython==0.29.32 --user
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

debug: install-deps install-meson-deps
	meson setup debug --python.install-env auto --buildtype="debug"
	meson compile -C debug 
	python -m pip install --no-build-isolation --editable .

cpp:
	@mkdir -p build
	g++ -g -Wall -std=c++17 -DCPP_STANDALONE -I./ -Igomea/ -IEigen/ gomea/src/discrete/*.cpp gomea/src/fitness/*.cpp gomea/src/common/*.cpp gomea/src/utils/*.cpp -o build/DiscreteGOMEA

install: install-deps build-wheel
	pip3 install dist/*.whl --user

reinstall: install-deps build-wheel
	pip3 install dist/*.whl --user --force-reinstall

src-install: install-deps build-sdist
	pip3 install dist/*.tar.gz --user

cibuildwheel: clean
	sudo pipx run cibuildwheel     

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
