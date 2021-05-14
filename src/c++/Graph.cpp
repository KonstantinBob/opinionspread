/*
 * Graph.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: konbob
 */

#include "Graph.h"
#include "Histogram.h"
#include "RandomNumbers.h"

#include <armadillo>

#include <iostream>
#include <vector>
#include <exception>
#include <cmath>
#include <queue>
#include <fstream>
#include <random>
#include <algorithm>
#include <iterator>
#include <set>
#include <functional>

/*
 * Returns vector of nodenumbers of adjacent nodes that are connected by outgoing edges.
 * Complexity: O(1)
 */
std::vector<long int> Graph::getAdjacentOutNodes(const unsigned int& NodeNumber) const {
	// check if node number is valid
	if (NodeNumber >= getSize()) {
		std::cout << "node number out of range at getAdjacentOutNodes! Exit. " << std::endl;
		exit(1);
	} else {
		try {
			return listOfOutNeighborlists.at(NodeNumber);
		} catch (std::exception& e) {
			std::cout << "Could not access neigbhor list at getAdjacentOutNodes! Exit. " << std::endl;
			exit(1);
		}

	}
}

/*
 * Returns vector of nodenumbers of adjacent nodes that are connected by ingoing edges.
 * Complexity: O(1)
 */
std::vector<long int> Graph::getAdjacentInNodes(const unsigned int& NodeNumber) const {
	// check if node number is valid
	if (NodeNumber >= getSize()) {
		std::cout << "node number out of range at getAdjacentInNodes! Exit. " << std::endl;
		exit(1);
	} else {
		try {
			return listOfInNeighborlists.at(NodeNumber);
		} catch (std::exception& e) {
			std::cout << "Could not access neigbhor list at getAdjacentInNodes! Exit. " << std::endl;
			exit(1);
		}
	}
}

/*
 * Returns vector of nodenumbers of adjacent nodes that are connected by ingoing and outgoing edges.
 * Complexity: O(1)
 */
std::vector<long int> Graph::getAdjacentAllNodes(const unsigned int& NodeNumber) const {
	// check if node number is valid
	if (NodeNumber >= getSize()) {
		std::cout << "node number out of range at getAdjacentAllNodes! Exit. " << std::endl;
		exit(1);
	} else {
		try {
			std::vector<long int> neighborhoodIn = listOfInNeighborlists.at(NodeNumber);
			std::vector<long int> neighborhoodOut = listOfOutNeighborlists.at(NodeNumber);

			if (neighborhoodOut.size() > 0 && neighborhoodIn.size() > 0) {
				neighborhoodOut.insert(neighborhoodOut.end(), neighborhoodIn.begin(), neighborhoodIn.end());
				return neighborhoodOut;
			} else if (neighborhoodOut.size() > 0) {
				return neighborhoodOut;
			} else if (neighborhoodIn.size() > 0) {
				return neighborhoodIn;
			} else {
				std::cout << "Disconnected node (" << NodeNumber << ") detected! Exit. " << std::endl;
				exit(1);
			}

		} catch (std::exception& e) {
			std::cout << "Could not access neigbhor list at getAdjacentAllNodes! Exit. \n " << e.what() << std::endl;

			exit(1);
		}
	}
}

/*
 * Returns the number of nodes as the size of the adjacency matrix
 * Complexity: O(1)
 */
unsigned long int Graph::getNumberOfNodes() const {
	return getSize();
}
/*
 * Returns the number of edges in the graph by iterating over the adjacency matrix.
 * NOTE: Undirected Edges will be counted twice!
 * Complexity: O(1)
 */
unsigned long int Graph::getNumberOfEdges() const {
	//Undirected edges will be counted twice!
	return adjacencyMatrix.n_nonzero;
}

/*
 * Returns a histogram of edges by iterating over the adjacency matrix and summing the entries
 */
Histogram Graph::getDegreeDistribution() const {
	// vector for the counts
	std::vector<double> counts;
	counts.assign(getSize(), 0);

	// loop over the matrix
	for (unsigned long int nodeIndex = 0; nodeIndex < getSize(); ++nodeIndex) {
		counts.at(nodeIndex) = adjacencyMatrix.col(nodeIndex).n_nonzero + adjacencyMatrix.row(nodeIndex).n_nonzero;
	}

	double binwidth = 1.0;
	Histogram degreeHisto(counts, binwidth);
	return degreeHisto;

}

/*
 * Returns a list of distances in shortest path sense to all nodes from source node given by NodeNumber.
 * A breadth first algorithm is used.
 */
std::vector<double> Graph::getShortestPathDistances(const unsigned int& SourceNodeNumber) const {
	if (SourceNodeNumber < getNumberOfNodes()) {
		// construct the lists for distances and visited nodes
		std::vector<double> distance(getNumberOfNodes(), std::numeric_limits<double>::infinity());
		std::vector<bool> visited(getNumberOfNodes(), false);

		// construct empty queue. Nodes will be managed by node number.
		std::queue<unsigned int> nodeQueue;

		// initialize the source node correctly
		distance.at(SourceNodeNumber) = 0.0;
		visited.at(SourceNodeNumber) = true;
		nodeQueue.push(SourceNodeNumber);

		// perform breadth first search
		while (!nodeQueue.empty()) {
			// dequeue node
			unsigned int currentNodeNumber = nodeQueue.front();
			nodeQueue.pop();

			// check out the neighbors
			for (auto & neighbor : getAdjacentAllNodes(currentNodeNumber)) {
				if (visited.at(neighbor) == false) {
					visited.at(neighbor) = true;
					distance.at(neighbor) = distance.at(currentNodeNumber) + 1;
					nodeQueue.push(neighbor);
				}
			}
		}

		// clean up infinities and zeros by removing them
		distance.erase(std::remove(distance.begin(), distance.end(), std::numeric_limits<double>::infinity()), distance.end());
		distance.erase(std::remove(distance.begin(), distance.end(), 0.0), distance.end());

		// return distances
		return distance;

	} else {
		std::cout << "node number to big. Exiting." << std::endl;
		exit(1);
	}

}
/*
 * Returns an approximation of the effective diameter by performing numberOfTrails breadth first searches.
 * The 90% percentile of the total distance histogram and its variance estimated by bootstraping are reported.
 * return values
 * vec.at(0): estimate ( = mean of distribution ) of effective diameter
 * vec.at(1): standard deviation of distribution of effective diameter
 */
std::vector<double> Graph::getApproximateEffectiveDiameter(const long int& numberOfTrails, const long int& numberOfResamples) const {

	if (numberOfTrails > 0 && numberOfResamples > 0) {

		// ========== get raw data from graph =======================

		// declaration outside the loop
		unsigned int currentNodeNumber = 0;
		std::vector<double> sampledDistances;
		std::vector<double> totalDistances;
		std::vector<double> results(2, 0.0);

		// increase capacity at the beginning to avoid copying
		totalDistances.reserve(numberOfTrails * getNumberOfNodes());

		// get distances from random nodes
		for (long int trailCount = 0; trailCount < numberOfTrails; ++trailCount) {
			// pick random node
			currentNodeNumber = RandomNumbers::getRandomLongInt(0, getNumberOfNodes());

			// perform breadth first search and store results
			sampledDistances = getShortestPathDistances(currentNodeNumber);
			totalDistances.insert(totalDistances.end(), std::make_move_iterator(sampledDistances.begin()), std::make_move_iterator(sampledDistances.end()));
		}

		// ========== do bootstrap resampling to get the distribution of the effective diameter estimate ===========

		// declaration outside the loop
		std::vector<double> bootstrapedDistances;
		std::vector<double> bootstrapedEstimates;
		long int randomIndex;

		// increase capacity at the beginning to avoid copying
		bootstrapedDistances.reserve(totalDistances.size());
		bootstrapedEstimates.reserve(numberOfResamples);

		// resample the distances
		for (long int bootstrapCount = 0; bootstrapCount < numberOfResamples; ++bootstrapCount) {
			// clear the bootstraped sample
			bootstrapedDistances.clear();

			// sample from the real data set
			for (unsigned int sampleCounter = 0; sampleCounter < totalDistances.size(); ++sampleCounter) {
				// get random index
				randomIndex = RandomNumbers::getRandomLongInt(0, totalDistances.size());

				// put entry into sample
				bootstrapedDistances.push_back(totalDistances.at(randomIndex));
			}

			// calculate the 90% percentile
			Histogram diameterDistribution(bootstrapedDistances, 1.0);
			bootstrapedEstimates.push_back(diameterDistribution.getPercentile(0.9));

		}

		// =============== calculate mean and variance estimators of the effective diameter estimate ===========================

		Histogram histo(bootstrapedEstimates, 1e-1);

		results.at(0) = histo.getMeanEstimator();
		results.at(1) = sqrt(histo.getVarianceEstimator());

		return results;

	} else {
		std::cerr << " Invalid arguments at getApproximateEffectiveDiameter. Exit." << std::endl;
		exit(1);
	}
}

/*
 * Returns the effective diameter of the graph.
 * ( That is the 90% percentile of the distance histogram)
 * For construction of the distance matrix the Floyd–Warshall algorithm is used.
 * NOTE: Disconnected nodes will be treated as if they were a distance of 0 apart.
 */
double Graph::getEffectiveDiameter() const {

	// get the size of the matrix
	long int matrixSize = getSize();

	// construct an empty distance matrix
	arma::mat distanceMatrix(matrixSize, matrixSize);

	// fill it for Floyd–Warshall algorithm
	for (long int k = 0; k < matrixSize; ++k) {
		for (long int i = 0; i < matrixSize; ++i) {
			//if (i!=k) {
			double entry = adjacencyMatrix(k, i);
			if (entry != 0) {
				distanceMatrix(k, i) = entry;
			} else {
				distanceMatrix(k, i) = std::numeric_limits<double>::infinity();
			}
			//}
		}
	}

	// do the loops of Floyd–Warshall algorithm
	for (int k = 0; k < matrixSize; ++k) {
		for (int i = 0; i < matrixSize; ++i) {
			for (int j = 0; j < matrixSize; ++j) {
				if (distanceMatrix(i, j) > distanceMatrix(i, k) + distanceMatrix(k, j)) {
					distanceMatrix(i, j) = distanceMatrix(i, k) + distanceMatrix(k, j);
				}
			}
		}
	}

	// create vector for the histogram
	std::vector<double> distanceList;
	// increase capacity
	distanceList.reserve(matrixSize * matrixSize);

	// loop over for the histogram
	for (long int k = 0; k < matrixSize; ++k) {
		for (long int i = 0; i < matrixSize; ++i) {
			// remove zeros or infinities
			if (distanceMatrix(i, k) != 0.0 && distanceMatrix(i, k) != std::numeric_limits<double>::infinity()) {
				distanceList.push_back(distanceMatrix(i, k));
			}
		}
	}

	Histogram diameterDistribution(distanceList, 1.0);

	return diameterDistribution.getPercentile(0.9);
}

/*
 * Writes a file in dot language to draw the graph by dot, neato or some other program from graphviz package.
 */
void Graph::writeDotFile(const std::string& filename) const {

	std::ofstream file;
	try {
		// Open file
		file.open(filename, std::ofstream::out);

		// write preamble
		file << "Digraph {" << std::endl;
		file << "node [fillcolor=\"#000080\"]" << "node [color=\"#EEEEEE\"]\n"
		//<< "node [label=\"\"]\n"
				<< "node [shape=\"circle\"]\n"
				//<< "node [style=filled]\n"
				<< "node [width=0.1]\n" << "edge [color=\"#ff6600\"]\n" << "edge [arrowsize=0.6]\n" << "overlap=false\n" << "rankdir=\"LR\";" << std::endl;

		//write data to file
		for (unsigned int nodeNumber = 0; nodeNumber < getNumberOfNodes(); ++nodeNumber) {
			for (unsigned int adjacentNodeNumber = 0; adjacentNodeNumber < getNumberOfNodes(); ++adjacentNodeNumber) {
				// if nodes are connected write edge to file. Self loops are ignored
				if (nodeNumber != adjacentNodeNumber && adjacencyMatrix(nodeNumber, adjacentNodeNumber) != 0) {
					file << nodeNumber << " -> " << adjacentNodeNumber << std::endl;
				}
			}

		}

		// for case of single node
		if (getNumberOfNodes() == 1) {
			file << "0" << std::endl;
		}

		// write closing
		file << "}" << std::endl;

		// close file
		file.close();
		std::cout << "Finished writing graph to " << filename << std::endl;

	} catch (...) {
		std::cout << "Could not open or write to file " << filename << std::endl;
	}
}

/*
 * Reads graph from a file created by snap.
 */
void Graph::readGraphFromSnapFile(const std::string& filename) {
	try {

		// file object for reading
		std::ifstream file;

		// open the file
		file.open(filename, std::ifstream::in);

		// skip the first two lines
		std::string line;
		std::getline(file, line);
		std::getline(file, line);

		// ===== parse the third line to get the number of nodes =====
		// get the line
		std::getline(file, line);

		// split it into parts at delimeter " "
		std::stringstream stringstream;
		stringstream.str(line);
		std::string substring;
		std::vector<std::string> parsedEntries;
		while (std::getline(stringstream, substring, ' ')) {
			parsedEntries.push_back(substring);
		}

		// use the third entry of the vector for number of nodes
		long int readOutSize = std::stol(parsedEntries.at(2));

		// skip the fourth line
		std::getline(file, line);

		// check if graph is valid
		if (readOutSize > 0) {
			size = readOutSize;

			// create empty adjacency matrix
			adjacencyMatrix = arma::SpMat<short>(size, size);

			// === read the following lines ================
			// get the node numbers
			long int a, b;
			while (file >> a >> b) {

				// add entry
				createEdge(a, b);
			}

			// close file
			file.close();
		} else {
			std::cerr << " Problem with matrix size." << std::endl;
			exit(1);
		}

	} catch (std::exception& e) {
		std::cerr << " Problem occurred in readFromGraphFile. Is the file missing?" << std::endl;
		std::cerr << e.what() << std::endl;
		exit(1);
	}
}

/*
 * Writes a file in dot language to draw the graph by dot, neato or some other program from graphviz package.
 * Every node can have its own color, encoded in argument vector color Data.
 */
void Graph::writeDotFileWithColorData(const std::string& filename, const std::vector<std::string>& colorData, const std::string& appendix) const {
	std::ofstream file;
	try {
		// Open file
		file.open(filename, std::ofstream::out);

		file << "Digraph {" << std::endl;
		file << "node [label=\"\"]\n" << "node [style=filled]\n" << "node [shape=\"circle\"]\n" << "overlap=false\n" << "edge [color=\"#ff6600\"]\n" << "rankdir=\"LR\";"
				<< "edge [arrowsize=0.6]\n" << std::endl;

		//write data to file
		for (unsigned int nodeNumber = 0; nodeNumber < getNumberOfNodes(); ++nodeNumber) {
			for (unsigned int adjacentNodeNumber = 0; adjacentNodeNumber < getNumberOfNodes(); ++adjacentNodeNumber) {
				// if nodes are connected write edge to file. Self loops are ignored
				if (nodeNumber != adjacentNodeNumber && adjacencyMatrix(nodeNumber, adjacentNodeNumber) != 0) {
					file << nodeNumber << " -> " << adjacentNodeNumber << std::endl;
				}
			}
			// write color value for node
			file << nodeNumber << " [color=\" " << colorData.at(nodeNumber) << " \"]" << std::endl;

		}

		// write appendix
		file << appendix << std::endl;

		// write closing
		file << "}" << std::endl;

		// close file
		file.close();
		std::cerr << "Finished writing graph to " << filename << std::endl;

	} catch (...) {
		std::cerr << "Could not open or write to file " << filename << " or index out of bounds." << std::endl;
	}
}

/*
 * Picks up to <numberNeighbors> random numbers of list <Neighborhood>.
 * Returns a vector of node numbers.
 * non const: Vector in argument will be shuffled!
 */
std::vector<long int> Graph::getRandomNeighbors(const std::vector<long int>& Neighborhood, const long int& numberNeighbors, const std::set<long int>& visitedNodeSet) const {

	if (numberNeighbors > 0) {

		// storage of nodes
		std::vector<long int> pickedNeighbors;
		pickedNeighbors.reserve(Neighborhood.size());

		// add only those that have not already been visited
		for (const long int & nodeNumber : Neighborhood) {
			if (visitedNodeSet.count(nodeNumber) == 0) {
				pickedNeighbors.push_back(nodeNumber);
			}
		}

		// pick from neighborhood without replacement by shuffling the neighborhood
		std::random_shuffle(pickedNeighbors.begin(), pickedNeighbors.end()); // FIXME standard generator

		//calculate the maximum number
		long int maxNeighbors;
		if (Neighborhood.size() <= (unsigned) numberNeighbors) {
			maxNeighbors = Neighborhood.size();
		} else {
			maxNeighbors = numberNeighbors;
		}

		// remove all others
		for (long int counter = pickedNeighbors.size(); counter > maxNeighbors; --counter) {
			pickedNeighbors.pop_back();
		}
		pickedNeighbors.shrink_to_fit();

		return pickedNeighbors;

	} else {
		std::cerr << "number of Neighbors <= 0. Exit." << std::endl;
		exit(1);
	}
}

/*
 * Create an edge between the two given nodes
 */
void Graph::createEdge(const long int& sourceNodeNumber, const long int& targetNodeNumber) {
	adjacencyMatrix(sourceNodeNumber, targetNodeNumber) = 1;
}

/*
 * Removes the  edge between the two given nodes
 */
void Graph::removeEdge(const long int& sourceNodeNumber, const long int& targetNodeNumber) {
	adjacencyMatrix(sourceNodeNumber, targetNodeNumber) = 0;
}

std::vector<long int> Graph::findDisconnectedNodes() const {
	// list of disconnected nodes
	std::vector<long int> disconnectedNodes;

	// constant and flag
	long int numberOfNodes = getSize();
	bool isDisconnected;

	// loop over adjacency matrix
	for (long int currentNode = 0; currentNode < numberOfNodes; currentNode++) {
		isDisconnected = true;
		for (long int nodeToCompare = 0; nodeToCompare < numberOfNodes; nodeToCompare++) {
			if (adjacencyMatrix(currentNode, nodeToCompare) == 1.0) {
				isDisconnected = false;
				break;
			}
			if (adjacencyMatrix(nodeToCompare, currentNode) == 1.0) {
				isDisconnected = false;
				break;
			}

		}
		if (isDisconnected) {
			disconnectedNodes.push_back(currentNode);
		}
	}

	return disconnectedNodes;
}

void Graph::printAdjacencyMatrix(const std::string& message) const {
	std::cerr << message << std::endl;

	// loop through row of adjacency matrix
	for (unsigned int columnIndex = 0; columnIndex < getSize(); ++columnIndex) {
		std::cerr << "[";
		for (unsigned int rowIndex = 0; rowIndex < getSize(); rowIndex++) {
			if (adjacencyMatrix(rowIndex, columnIndex) == 0.0) {
				std::cerr << " 0 ";
			}
			if (adjacencyMatrix(rowIndex, columnIndex) == 1.0) {
				std::cerr << " 1 ";
			}
		}
		std::cerr << "]" << std::endl;
	}
}

/*
 * Create the lists of outgoing  and ingoing neighbors. To prevent disconnected nodes, add a random edge in case.
 */
void Graph::createNeighborlists(const bool& fixUnconnected) {
	// declare and reserve space for listOfNeighborlists
	listOfInNeighborlists = std::vector<std::vector<long int>>(getSize(), std::vector<long int>(0));
	listOfOutNeighborlists = std::vector<std::vector<long int>>(getSize(), std::vector<long int>(0));

	// loop over all nodes
	for (unsigned long int currentNode = 0; currentNode < getSize(); ++currentNode) {
		// ====== ingoing neighbors =========
		// prepare vector for in nodes
		std::vector<long int> AdjacentInNodeNumbers = { };

		// loop through columns of adjacency matrix
		for (unsigned int index = 0; index < getSize(); index++) {
			if (adjacencyMatrix(index, currentNode) != 0.0 && index != currentNode) {
				AdjacentInNodeNumbers.push_back(index);
			}
		}

		listOfInNeighborlists.at(currentNode) = AdjacentInNodeNumbers;

		// ====== outgoing neighbors =========
		// prepare vector for in nodes
		std::vector<long int> AdjacentOutNodeNumbers = { };

		// loop through rows of adjacency matrix
		for (unsigned int index = 0; index < getSize(); index++) {
			if (adjacencyMatrix(currentNode, index) != 0.0 && index != currentNode) {
				AdjacentOutNodeNumbers.push_back(index);
			}
		}

		listOfOutNeighborlists.at(currentNode) = AdjacentOutNodeNumbers;

		// ====== fix unconnected nodes in random graphs =========
		if (fixUnconnected == true) {
			//check if node is disconnected
			if (AdjacentInNodeNumbers.size() == 0 && AdjacentOutNodeNumbers.size() == 0) {

				// add random edge
				long int randomNode = RandomNumbers::getRandomLongInt(0, getSize());
				createEdge(currentNode, randomNode);

				// write message
				std::cerr << "#I had to create an edge to prevent disconnected nodes in " << typeIdentifier << "!\n" << std::endl;

				// update neighborhoods
				updateOutNeighborlists(currentNode, randomNode);
				updateInNeighborlists(randomNode, currentNode);

			}
		}

	}
}

void Graph::updateInNeighborlists(const unsigned int& NodeNumber, const unsigned int& newNeighbor) {
	// get old neighborlist
	std::vector<long int> neighboorhood = listOfInNeighborlists.at(NodeNumber);

	// add new neighbor
	neighboorhood.push_back(newNeighbor);

	// store new neighboorhood
	listOfInNeighborlists.at(NodeNumber) = neighboorhood;
}

/*
 * Returns the distribution of cover times.
 * This is done by sending out a swarm surfers that move along every edge and measuring its cover time.
 */
Histogram Graph::getCoverTime(const long int& numberOfTrails) const {
	// store the results
	std::vector<double> results;
	results.reserve(numberOfTrails);

	for (long int trail = 0; trail < numberOfTrails; ++trail) {

		// prepare list of already visited nodes
		std::vector<bool> visitedList(getNumberOfNodes(), false);

		// queue of current surfer locations
		std::queue<long int> currentSurferLocations;

		// queue of next generation surfer locations
		std::queue<long int> nextSurferLocations;

		// start at random point in graph
		long int startNode = RandomNumbers::getRandomLongInt(0, getNumberOfNodes());

		// add it to the current locations
		currentSurferLocations.push(startNode);

		// counter for generations
		long int generations = 0;

		// counter for visited nodes
		long int visitedNodes = 0;

		// flag for breaking the loop
		bool breakLoop = false;

		// now proceed the next generations
		do {
			// ======== work on current nodes
			do {
				// get next node
				long int currentNode = currentSurferLocations.front();
				currentSurferLocations.pop(); // remove it

				// check if it has already been visited
				if (visitedList.at(currentNode) == false) {
					// then add the neighbors
					for (long int neighbor : getAdjacentAllNodes(currentNode)) {
						nextSurferLocations.push(neighbor);
					}

					// mark current node as visited
					visitedList.at(currentNode) = true;
					++visitedNodes;
				}
				++generations;

			} while (!currentSurferLocations.empty());

			// ===== swap current and next generation
			// swap them
			currentSurferLocations = nextSurferLocations;

			// make it empty
			std::queue<long int> empty;
			std::swap(nextSurferLocations, empty);

			// ==== check break flag ====
			if (nextSurferLocations.empty() && currentSurferLocations.empty()) {
				breakLoop = true;
			}
			if ((double) visitedNodes > 0.9 * getNumberOfNodes()) {
				breakLoop = true;
			}

		} while (!breakLoop);

		// check if the results are reliable
		if ((double) visitedNodes > 0.75 * getNumberOfNodes()) {
			results.push_back(generations);
		}

	}

	// return histogram
	if (results.size() == 0) {
		std::cerr << " No trail successful. Exit" << std::endl;
		exit(1);
	}

	Histogram histo(results, 1.0);
	return histo;
}

void Graph::updateOutNeighborlists(const unsigned int& NodeNumber, const unsigned int& newNeighbor) {
	// get old neighborlist
	std::vector<long int> neighboorhood = listOfOutNeighborlists.at(NodeNumber);

	// add new neighbor
	neighboorhood.push_back(newNeighbor);

	// store new neighboorhood
	listOfOutNeighborlists.at(NodeNumber) = neighboorhood;
}

Graph::~Graph() {

}

unsigned long int Graph::getSize() const {
	return adjacencyMatrix.n_rows;
}

long int Graph::getNumberOfOutNodes(const unsigned int& NodeNumber) const {
	return listOfOutNeighborlists.at(NodeNumber).size();
}

const std::string& Graph::getTypeIdentifier() const {
	return typeIdentifier;
}

long int Graph::getNumberOfInNodes(const unsigned int& NodeNumber) const {
	return listOfInNeighborlists.at(NodeNumber).size();
}

/*
 * Count the number of nodes that are not symmetrically connected, i.e. do not have both in- and outgoing edges.
 */
unsigned long int Graph::countUnsymmetricalConnectedNodes() const {
	// prepare counter
	unsigned long int counter = 0;

	//loop over nodes
	for (unsigned int nodeNumber = 0; nodeNumber < getNumberOfNodes(); ++nodeNumber) {

		// check if there are outgoing egdes
		if (getAdjacentOutNodes(nodeNumber).size() > 0) {

			// compare to the ingoing edges
			if (getAdjacentOutNodes(nodeNumber).size() != getAdjacentInNodes(nodeNumber).size()) {
				counter++;
			}
		} else {
			// If there is no outgoing edge there must be a least one ingoing. Therefore unsymmetrically.
			counter++;
		}

	}

	return counter;
}

/*
 * Writes a dot file for a graph with given weights.
 * The index in the array corresponds to the node number.
 */
void Graph::writeDotFileWithWeights(const std::string& filename, const std::vector<double>& weights) const {

	std::ofstream file;
	try {
		// Open file
		file.open(filename, std::ofstream::out);

		file << "Graph {" << std::endl;
		file << "node [shape=\"circle\"]\n" << "overlap=false\n" << " splines=false \n" << "edge [color=\"#ff6600\"]\n" << "rankdir=\"LR\";" << "edge [arrowsize=0.6]\n"
				<< std::endl;

		//write data to file
		for (unsigned int nodeNumber = 0; nodeNumber < getNumberOfNodes(); ++nodeNumber) {

			// write color value for node
			file << nodeNumber << " [width=\" " << 0.05 * weights.at(nodeNumber) << " \"]" << std::endl;

			for (unsigned int adjacentNodeNumber = 0; adjacentNodeNumber < getNumberOfNodes(); ++adjacentNodeNumber) {
				// if nodes are connected write edge to file. Self loops are ignored
				if (nodeNumber != adjacentNodeNumber && adjacencyMatrix(nodeNumber, adjacentNodeNumber) != 0) {
					file << nodeNumber << " -- " << adjacentNodeNumber << std::endl;
				}
			}

		}

		// write closing
		file << "}" << std::endl;

		// close file
		file.close();
		std::cerr << "Finished writing graph to " << filename << std::endl;

	} catch (...) {
		std::cerr << "Could not open or write to file " << filename << " or index out of bounds." << std::endl;
	}
}

const arma::SpMat<short>& Graph::getAdjacencyMatrix() const {
	return adjacencyMatrix;
}

/*
 * Computes the fraction of unsymmetrical edges.
 */

double Graph::getSkewness() const {

	double counterUpper = 0.0;
	double counterLower = 0.0;

	for (long int columnIndex = 0; columnIndex < size; columnIndex++) {
		for (long int rowIndex = columnIndex; rowIndex < size; rowIndex++) {

			// count if edge is in upper triangle
			if (adjacencyMatrix(columnIndex, rowIndex) == 1) {
				++counterUpper;
			}

			// count if edge is in lower triangle
			if (adjacencyMatrix(rowIndex, columnIndex) == 1) {
				++counterLower;
			}

		}
	}


	return abs(counterLower-counterUpper) / (counterLower + counterUpper);
}
