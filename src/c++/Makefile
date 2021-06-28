FLAGS =-Wall -Wextra -std=c++11 -O3 -funroll-loops -finline -march=native -larmadillo -fdiagnostics-color -fPIC

speed: Main

debug: FLAGS +=-g -pg -time
debug: Main

# -MM for dependencies
	
Graph.o: Graph.cpp Graph.h Histogram.o RandomNumbers.o
	g++ $(FLAGS) -c $<	
		
graphs = ForestFire.o Erdos.o StochKronecker.o BarabasiAlbert.o Arxiv.o MetaGraph.o WattsStrogatz.o

$(graphs): %.o: %.cpp %.h Graph.o
	g++ $(FLAGS) -c $< -o $@
	
Main.o: Main.cpp Main.h LearningNetwork.hpp Clustering.hpp
	g++ $(FLAGS) -c $<
	
LearningNetwork.hpp: Scientist.o ForestFire.o Erdos.o StochKronecker.o Histogram.o Graph.o

Clustering.hpp: LearningNetwork.hpp
		
%.o: %.cpp %.h
	g++ $(FLAGS) -c $<
	
Main: Main.o Graph.o Histogram.o RandomNumbers.o Scientist.o $(graphs)
	g++ -o Echse Main.o Graph.o Histogram.o RandomNumbers.o Scientist.o $(graphs) $(FLAGS)