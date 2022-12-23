default: 
	python setup.py build_ext --inplace

RVGOMEA-cpp: GOMEAlib-py
	g++ -g -Wall -std=c++17 $(wildcard src/gomea/cpp/*.cpp) -o build/RealValuedGOMEA -I/usr/include/python3.8/ -Iinclude/ -Isrc/ -L/usr/lib64/ -L./ -lpython3.8 -l:RealValuedGOMEA.cpython-38-x86_64-linux-gnu.so

debug:
	python setup.py build_ext --inplace --debug -j4

install:
	python setup.py bdist_wheel
	pip install dist/*.whl --user --force-reinstall

pip-dist:
	python setup.py sdist bdist_wheel
	auditwheel repair $(`ls dist/*-cp38-cp38-linux_x86_64.whl`)
	rm dist/*-cp38-cp38-linux_86_64.whl
	mv wheelhouse/* dist/

clean:
	rm -f *.so
	rm -rf gomea.egg-info/
	rm -rf build/
	rm -rf dist/
	rm -f gomea/*.cpp
	rm -f gomea/*.h
	rm -f gomea/*.so
	rm -rf gomea/__pycache__/
