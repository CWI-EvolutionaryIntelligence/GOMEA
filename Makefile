all:
	python setup.py build_ext --inplace

install:
	python setup.py install --user
	
	#python setup.py sdist bdist_wheel
	#auditwheel repair $(`ls dist/*-cp38-cp38-linux_x86_64.whl`)
	#rm dist/*-cp38-cp38-linux_86_64.whl
	#mv wheelhouse/* dist/

clean:
	rm -rf build/
	rm -rf dist/
	rm -f src/cython/*.cpp
	rm -f src/cython/*.c
