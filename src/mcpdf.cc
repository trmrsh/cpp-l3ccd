/*

!!begin
!!title   mcpdf -- Monte Carlo generator of L3CCD output events
!!author  T.R. Marsh
!!created 12 Dec 2006
!!created 29 Dec 2006
!!descr   Monte Carlo generator of L3CCD output events
!!root    mcpdf
!!index   mcpdf
!!class   Programs
!!css     style.css
!!head1   mcpdf -- Monte Carlo generator of L3CCD output events
 

Program to generate the PDF of an avalanche register by brute force Monte
Carlo approach.  i.e. it follows each electron through every stage. It is
too slow to be of much practical use but provide an independent check of more 
sophisticated methods. It allows for CICs generated within the register (i.e. 
output electrons even in the absence of any input electrons) and readout noise added after it. Since
CICs generated in the avalanche register do not run through the full set of
stages they have a different probability distribution from the output
generated from genuine input electrons.

!!head2 Invocation

mcpdf nstep nmax gain increase pcic nelec poisson read bias cfac nbin x1 x2 nevent
seed nrandom

!!head2 Arguments

!!table 

!!arg{nstep}{Number of avalanche multiplication steps.}

!!arg{nmax}{The maximum output number to consider.}  

!!arg{gain}{Mean gain from start to end of the avalanche register (sets the
multiplication probability per stage).}

!!arg{increase}{Factor increase (can be < 1 to give a decrease) in
multiplication probability from input to output of register. This is to assess
the influence of variable multiplication probabilities on the statistics. One
might by default anticipate a higher probability near the output from where
the clocking voltages are driven. The probabilities increase logarithmically
in the direction of the electron flow.  }

!!arg{pcic}{CIC probability per stage,
i.e. the chance that an electron is created spontaneously.}

!!arg{nelec}{Number of electrons at input} 

!!arg{poisson}{Assume 'nelec' represents a mean with a Poisson distribution 
or an exact number (nearest integer will be taken).}  

!!arg{read}{read noise in electrons assumed to be added on ouput. 'read' is 
the RMS of this, which is assumed to be gaussian.}  

!!arg{bias}{Bias level, i.e. what zero will become in the final output}

!!arg{cfac}{Factor to convert electrons at the end of the register into counts (ADU)}

!!arg{nbin}{Number of output bins}

!!arg{x1}{Left-edge of first output pixel. Its up to the user to ensure that it overlaps with nmax OK}

!!arg{x2}{Right-edge of last output pixel.}

!!arg{nevent}{Number of events to generate.} 

!!arg{seed}{seed integer for random number generator.} 

!!arg{nrandom}{Which of four random number generators to pick (1 to 4)}

!!table

!!end

*/

#include <cstdlib>
#include <cfloat>
#include <climits>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/array1d.h"
#include "trm/format.h"
#include "trm/input.h"
#include "trm/l3ccd.h"

int main (int argc, char *argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, L3ccd::L3CCD_ENV, L3ccd::L3CCD_DIR);

    // Define inputs
    input.sign_in("nstep",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nmax",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("gain",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("increase", Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("pcic",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nelec",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("poisson",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("read",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("bias",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("cfac",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nbin",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",       Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x2",       Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nevent",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("seed",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nrandom",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    int nstep;
    input.get_value("nstep", nstep, 591, 1, 10000, "number of multiplication steps");
    double gain;
    input.get_value("gain", gain, 100., 1., 1000000., "mean gain");
    double increase;
    input.get_value("increase", increase, 1., 0.001, 1000., "probability increase factor");

    double pstart = L3ccd::start_prob(nstep, increase, gain);

    double pcic;
    input.get_value("pcic", pcic, 0., 0., 1., "CIC multiplication probability per stage");
    double nelec;
    input.get_value("nelec", nelec, 1., 0., 1000., "number of input electrons");
    bool poisson;
    input.get_value("poisson", poisson, true, "account for poisson distribution on input?");
    double read;
    input.get_value("read", read, 0., 0., 1000000., "read noise (electrons)");
    double bias;
    input.get_value("bias", bias, 0., -DBL_MAX, DBL_MAX, "output bias level");
    double cfac;
    input.get_value("cfac", cfac, 1., DBL_MIN, DBL_MAX, "electrons at output/ADU");
    int nbin;
    input.get_value("nbin", nbin, 1000, 1, 100000000, "number of output bins");
    double x1;
    input.get_value("x1", x1, 0.,  -DBL_MAX, DBL_MAX, "left-edge of first bin");
    double x2;
    input.get_value("x2", x2, 100., -DBL_MAX, DBL_MAX, "right-edge of last bin");
    int nevent;
    input.get_value("nevent", nevent, 1000, 1, INT_MAX, "number of events to generate");
    Subs::INT4 seed;
    input.get_value("seed", seed, Subs::INT4(-79791), Subs::INT4(INT_MIN), Subs::INT4(INT_MAX), 
		    "seed integer for random number generator");
    if(seed > 0) seed = -seed;
    int nrandom;
    input.get_value("nrandom", nrandom, 3, 1, 4, "which random number generator to use");

    input.save();

    // Calculate the variance
    double vr = 0., gr = 1., mcic = 0., vcic = 0.;
    for(int i=1; i<=nstep; i++){

	// Compute multiplication probability
	double pmult  = exp(log(pstart)  + log(increase)*(nstep-i)/(nstep-1));

	mcic += gr;
	vcic += pcic*vr      + pcic*(1-pcic)*gr*gr;

	vr    = (1+pmult)*vr + pmult*(1-pmult)*gr*gr;
	gr   *= 1+pmult;
    }
    mcic *= pcic;

    double mui = poisson ? nelec : Subs::nint(nelec);

    // print out the parameters
    Subs::Format form(12);
    std::cout << "#\n";
    std::cout << "# Generated by mcpdf with the following inputs:\n";
    std::cout << "#\n";
    std::cout << "# nstep    = " << nstep  << std::endl;
    std::cout << "# gain     = " << form(gain)    << std::endl;
    std::cout << "# increase = " << form(increase)<< std::endl;
    std::cout << "# pcic     = " << form(pcic)    << std::endl;
    std::cout << "# nelec    = " << form(nelec)   << std::endl;
    std::cout << "# poisson  = " << poisson << std::endl;
    std::cout << "# read     = " << form(read)    << std::endl;
    std::cout << "# bias     = " << form(bias)    << std::endl;
    std::cout << "# cfac     = " << form(cfac)    << std::endl;
    std::cout << "# nbin     = " << nbin        << std::endl;
    std::cout << "# x1       = " << form(x1)    << std::endl;
    std::cout << "# x2       = " << form(x2)    << std::endl;
    std::cout << "# nevent   = " << nevent  << std::endl;
    std::cout << "# seed     = " << seed  << std::endl;
    std::cout << "# nrandom  = " << nrandom << std::endl;
    std::cout << "#\n";
    std::cout << "# Derived quantities:\n";
    std::cout << "#\n";
    std::cout << "# Start multiplication probability   = " << form(pstart)          << std::endl;
    std::cout << "# End multiplication probability     = " << form(increase*pstart) << std::endl;
    std::cout << "# Directly measured gain             = " << form(gr/cfac)           << std::endl;
    std::cout << "# Mean for 1 electron input, no CICs = " << form(bias+gr/cfac)        << std::endl;
    std::cout << "# RMS for 1 electron input, no CICs  = " << form(sqrt(vr)/cfac)        << std::endl;
    std::cout << "# Mean for 0 electron input          = " << form(bias+mcic/cfac)      << std::endl;
    std::cout << "# Overall mean                       = " << form(bias + (mcic + gain*mui)/cfac) << std::endl;
    if(poisson)
	std::cout << "# Overall RMS                        = " << form(sqrt(read*read + vcic + (gain*gain + vr)*mui)/cfac)      << std::endl;
    else
	std::cout << "# Overall RMS                        = " << form(sqrt(read*read + vcic + vr*mui)/cfac)      << std::endl;

    std::cout << "#\n";

    // Function pointers to select between possibilities
    float (*poiss)(float, Subs::INT4&);
    double (*gauss)(Subs::INT4&);
    double (*ran)(Subs::INT4&);
    switch(nrandom){
	case 1:
	    poiss   = &Subs::poisson1;
	    gauss   = &Subs::gauss1;
	    ran     = &Subs::ran1;
	    break;
	case 2:
	    poiss   = &Subs::poisson2;
	    gauss   = &Subs::gauss2;
	    ran     = &Subs::ran2;
	    break;
	case 3:
	    poiss   = &Subs::poisson3;
	    gauss   = &Subs::gauss3;
	    ran     = &Subs::ran3;
	    break;
	case 4:
	    poiss   = &Subs::poisson4;
	    gauss   = &Subs::gauss4;
	    ran     = &Subs::ran4;
	    break;
	default:
	    throw L3ccd::L3ccd_Error("nrandom out of range");
    }

    Subs::Array1D<int> hbuff(nbin);
    hbuff = 0;
    
    int n_out_of_range = 0;
    for(int i=0; i<nevent; i++){

	// Initial number of electrons
	int ne;
	if(poisson)
	    ne = int(poiss(float(nelec), seed));
	else
	    ne = int(Subs::nint(nelec));

	// run through every stage
	for(int j=0; j<nstep; j++){

	    // Compute multiplication probability for step ns
	    double pmult  = exp(log(pstart)  + log(increase)*j/(nstep-1));

	    // go through each input electron for this stage
	    int nout = 0;
	    for(int k=0; k<ne; k++){
		if(ran(seed) < pmult)
		    nout += 2;
		else
		    nout++;
	    }

	    // CIC
	    if(ran(seed) < pcic) nout++;

	    // Feed the output of one stage to the input of the next
	    ne = nout;

	}

	// Add readout noise, digitise and bin in one line	
	int nb = int(Subs::nint(nbin*(Subs::nint(bias + (ne + read*gauss(seed))/cfac)-x1)/(x2-x1)-0.5));
	if(nb >= 0 && nb < nbin)
	    hbuff[nb]++;
	else
	    n_out_of_range++;

    }

    std::cerr << "Number out of range = " << n_out_of_range << std::endl;
    
    // Print out the results
    for(int i=0; i<nbin; i++)
	std::cout << form(x1+(x2-x1)*(i+0.5)/nbin) << " " << form(hbuff[i]/double(nevent)) << " " << form(sqrt(hbuff[i])/nevent) << std::endl;

  }
  catch(const std::string& err){
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }
}

