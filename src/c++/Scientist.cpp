/*
 * Scientist.cpp
 *
 *  Created on: Oct 10, 2016
 *      Author: konbob
 */

#include "Scientist.h"
#include "RandomNumbers.h"

#include <cmath>
#include <array>
#include <fstream>
#include <iostream>

/*
 * Constructor for scientist.
 */
Scientist::Scientist(long int ID, double selfConfidence):
ID(ID),
selfConfidence(selfConfidence)
{
	// choose random opinion
	opinion=RandomNumbers::getRandomBool();

}

/*
 * Getter for id
 */
long int Scientist::getId() const {
	return ID;
}

/*
 * Getter for self confidence
 */
double Scientist::getSelfConfidence() const {
	return selfConfidence;
}

/*
 * Setter for opinion.
 */
void Scientist::setOpinion(bool opinion) {
	this->opinion = opinion;
}

/*
 * Getter for opinion.
 */
bool Scientist::getOpinion() const {
	return this->opinion;
}
/*
 * Mutates the theory.
 */
void Scientist::doMutation(const double& expGap) {
	// throw coin for up or down
	if (RandomNumbers::getRandomBool()==false) {
		// new opinion is down. Make Metropolis check:
		if (RandomNumbers::getRandomReal() < expGap) {
			opinion=false;
		}
	} else {
		// new opinion is up.
		// Accept it instantly.
		opinion= true;
	}
}

/*
 * Similarity measure between two scientists. Return value must be in the range [0,1].
 */
double Scientist::getSimilarity(const Scientist& scientist1, const Scientist& scientist2){
	if (scientist1.getOpinion() == scientist2.getOpinion()){
		return 1;
	}
	return 0;
}


