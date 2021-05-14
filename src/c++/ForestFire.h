/*
 * ForestFire.h
 *
 *  Created on: Nov 30, 2016
 *      Author: konbob
 */

#ifndef FORESTFIRE_H_
#define FORESTFIRE_H_

#include "Graph.h"

class ForestFire: public Graph {
public:
	ForestFire(const double& forwardBuringProbability, const double& backwardBurningProbability, const long int& numberOfNodes);
	virtual void resize(const long int& numberOfNodes);

private:
	double forwardBuringProbability;
	double backwardBurningProbability;

	void createEdges(const long int& numberOfNodes);
};

#endif /* FORESTFIRE_H_ */
