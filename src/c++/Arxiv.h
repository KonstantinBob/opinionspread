/*
 * Arxiv.h
 *
 *  Created on: Dec 6, 2016
 *      Author: konbob
 */

#ifndef ARXIV_H_
#define ARXIV_H_

#include "Graph.h"


class Arxiv: public Graph {
public:
	Arxiv(const std::string& filename);
	void resize(const long int& numberOfNodes);

};

#endif /* ARXIV_H_ */
