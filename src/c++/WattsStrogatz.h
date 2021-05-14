/*
 * WattsStrogatz.h
 *
 *  Created on: Jul 19, 2017
 *      Author: konbob
 */

#ifndef WATTSSTROGATZ_H_
#define WATTSSTROGATZ_H_

#include "Graph.h"

class WattsStrogatz: public Graph {
public:
	WattsStrogatz(const long int& k, const double& pRewire, const long int& numberOfNodes);
	void resize(const long int& numberOfNodes);

private:
	long int k;
	double pRewire;

	void createEdges(const long int& numberOfNodes);
};

#endif /* WATTSSTROGATZ_H_ */
