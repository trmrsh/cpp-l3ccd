#include <cstdlib>
#include <cmath>
#include "trm/l3ccd.h"

/* Computes the start probability of a register in which the multiplication probability
 * is assumed to increase from start to end logarithmically.
 * \param nstep    the number of steps
 * \param increase the factor increase from first to last multiplication step, measured in the directions
 * that the electrons travel. Must be greater than 0.
 * \param gain the gain factor desired, must be greater than or equal to 1.
 * \return the multiplication probability of the first step.
 * \exception Throws exceptions in inputs are out of range or if the derived probability exceeds 1 at any point.
 */

double L3ccd::start_prob(int nstep, double increase, double gain){

    if(increase <= 0.)
	throw L3ccd_Error("L3ccd::start_prob: increase <= 0");

    if(gain < 1.)
	throw L3ccd_Error("L3ccd::start_prob: gain < 1");

    if(gain == 1.) return 0.;

    // Slightly tricky computation to normalise the probabilities
    double p1 = 0., l1 = 0., p2=0.0001;
    const double LGAIN_TEST = log(gain);
    double l2 = LGAIN_TEST - 1.;

    // First a value for p2 that gives a gain (=l2) larger than wanted
    while(l2 < LGAIN_TEST){
	l2 = 0.;
	for(int i=0; i<nstep; i++){
	    double pv = exp(log(p2)  + log(increase)*i/(nstep-1));
	    l2 += log(1.+pv);
	}
	if(l2 < LGAIN_TEST) p2 *= 2;
    }

    // Now refine the number with a binary chop
    while(l2-l1 > 1.e-10){
	double lm = 0., pm = (p1+p2)/2.;
	for(int i=0; i<nstep; i++){
	    double pv = exp(log(pm)  + log(increase)*i/(nstep-1));
	    lm += log(1.+pv);
	}
	if(lm < LGAIN_TEST){
	    p1 = pm;
	    l1 = lm;
	}else{
	    p2 = pm;
	    l2 = lm;
	}
    }

    double pstart = (p1+p2)/2.;
    if(pstart > 1. || pstart*increase > 1.)
	throw L3ccd_Error("L3ccd::start_prob: multiplication probability has to exceed 1 to match gain");

    return pstart;
}

