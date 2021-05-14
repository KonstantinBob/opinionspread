/*
 * BarabasiAlbert.h
 *
 *  Created on: Jul 17, 2017
 *      Author: konbob
 */

#ifndef BARABASIALBERT_H_
#define BARABASIALBERT_H_

#include "Graph.h"

class BarabasiAlbert: public Graph {
public:
	BarabasiAlbert(const long int & mZero, const long int & m, const long int& numberOfNodes);
	void resize(const long int& numberOfNodes);

private:
	long int mZero;
	long int m;

	void createEdges(const long int& numberOfNodes);
};

#endif /* BARABASIALBERT_H_ */
