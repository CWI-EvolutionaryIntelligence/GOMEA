
CXX=g++
CXXFLAGS=-O2
EXE=GOMEA

$(EXE): build/main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

build/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

clean:
	rm build/*.o
	rm GOMEA
