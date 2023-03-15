default: 
	python setup.py build_ext --inplace

cpp: install 
	g++ -g -Wall -std=c++17 $(python-config --embed --cflags) -I./ $(python-config --embed --includes) gomea/src/discrete/main.cpp -o build/DiscreteGOMEA -L~/.local/lib/python3.10/site-packages/gomea/ $(python-config --embed --ldflags) -l:discrete.cpython-310-x86_64-linux-gnu.so

debug:
	python setup.py build_ext --inplace --debug

install:
	pip -q install build --user
	pip -q install cython --user
	python -m build --wheel
	pip install dist/*.whl --user

reinstall:
	python setup.py bdist_wheel
	pip install dist/*.whl --user --force-reinstall

clean:
	rm -f *.so
	rm -rf gomea.egg-info/
	rm -rf build/
	rm -rf dist/
	rm -f gomea/*.cpp
	rm -f gomea/*.h
	rm -f gomea/*.so
	rm -rf gomea/__pycache__/
