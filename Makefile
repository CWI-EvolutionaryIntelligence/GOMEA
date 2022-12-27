default: 
	python setup.py build_ext --inplace

debug:
	python setup.py build_ext --inplace --debug -j4

install:
	python setup.py bdist_wheel
	pip install dist/*.whl --user --force-reinstall

clean:
	rm -f *.so
	rm -f *.dat
	rm -rf gomea.egg-info/
	rm -rf build/
	rm -rf dist/
	rm -f gomea/*.cpp
	rm -f gomea/*.h
	rm -f gomea/*.so
	rm -rf gomea/__pycache__/
