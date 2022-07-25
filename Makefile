all:
	python setup.py build_ext --inplace -j4

debug:
	python setup.py build_ext --inplace --debug -j4

install:
	python setup.py install --user

pip-dist:
	python setup.py sdist bdist_wheel
	auditwheel repair $(`ls dist/*-cp38-cp38-linux_x86_64.whl`)
	rm dist/*-cp38-cp38-linux_86_64.whl
	mv wheelhouse/* dist/

clean:
	rm -f *.so
	rm -rf build/
	rm -rf dist/
	rm -f modules/*/src/cython/*.cpp
	rm -f modules/*/src/cython/*.c
	rm -f modules/*/src/cython/*.h
