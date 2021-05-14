/*
 * Graph.h
 *
 *  Created on: Sep 14, 2016
 *      Author: konbob
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <armadillo>
#include "Histogram.h"
#include <vector>
#include <string>
#include <functional>
#include <set>

class Graph {
public:
	// creation methods
	virtual void resize(const long int& numberOfNodes) = 0;
	void readGraphFromSnapFile(const std::string& filename);

	// analysis methods
	double getEffectiveDiameter() const;
	std::vector<double> getApproximateEffectiveDiameter(const long int& numberOfTrails, const long int& numberOfResamples) const;
	unsigned long int getNumberOfEdges() const;
	unsigned long int getNumberOfNodes() const;
	Histogram getDegreeDistribution() const;
	Histogram getCoverTime(const long int& numberOfTrails) const;
	unsigned long int countUnsymmetricalConnectedNodes() const;
	double getSkewness() const;

	// getter
	unsigned long int getSize() const;
	const std::string& getTypeIdentifier() const;
	std::vector<long int> getAdjacentOutNodes(const unsigned int& NodeNumber) const;
	std::vector<long int> getAdjacentInNodes(const unsigned int& NodeNumber) const;
	std::vector<long int> getAdjacentAllNodes(const unsigned int& NodeNumber) const;
	std::vector<double> getShortestPathDistances(const unsigned int& NodeNumber) const;
	std::vector<long int> getRandomNeighbors(const std::vector<long int>& Neighborhood, const long int& numberNeighbors, const std::set<long int>& visitedNodeSet) const;
	long int getNumberOfInNodes(const unsigned int& NodeNumber) const;
	long int getNumberOfOutNodes(const unsigned int& NodeNumber) const;
	std::vector<long int> findDisconnectedNodes() const;
	const arma::SpMat<short>& getAdjacencyMatrix() const;

	// output methods
	void writeDotFile(const std::string& filename) const;
	void writeDotFileWithColorData(const std::string& filename, const std::vector<std::string>& colorData, const std::string& appendix) const;
	void writeDotFileWithWeights(const std::string& filename, const std::vector<double>& weights) const;
	void printAdjacencyMatrix(const std::string& message) const;

	// virtual destructor
	virtual ~Graph();


protected:
	// utility methods
	void createEdge(const long int& sourceNodeNumber,const long int& targetNodeNumber);
	void removeEdge(const long int& sourceNodeNumber,const long int& targetNodeNumber);
	void createNeighborlists(const bool& fixUnconnected);
	void updateInNeighborlists(const unsigned int& NodeNumber, const unsigned int& newNeighbor);
	void updateOutNeighborlists(const unsigned int& NodeNumber, const unsigned int& newNeighbor);

	// member
	arma::SpMat<short> adjacencyMatrix;
	std::vector<std::vector<long int>> listOfOutNeighborlists;
	std::vector<std::vector<long int>> listOfInNeighborlists;
	unsigned long int size;
	std::string typeIdentifier;
};

#endif /* GRAPH_H_ */
