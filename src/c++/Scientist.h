/*
 * Scientist.h
 *
 *  Created on: Oct 10, 2016
 *      Author: konbob
 */

#ifndef SCIENTIST_H_
#define SCIENTIST_H_

#include <array>
#include <string>

class Scientist {

private:
	long int ID;
	double selfConfidence;
	bool opinion;

public:
	// constructor
	Scientist(long int ID,double selfConfidence);

	//getter
	long int getId() const;
	double getSelfConfidence() const;
	bool getOpinion() const;

	// setter
	void setOpinion(bool opinion);

	// methods
	void doMutation(const double& expGap);
	static double getSimilarity(const Scientist& scientist1, const Scientist& scientist2);

};

#endif /* SCIENTIST_H_ */
