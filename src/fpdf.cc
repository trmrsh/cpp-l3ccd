/*

!!begin
!!title   fpdf - fast computation PDF of L3CCD output distribution
!!author  T.R. Marsh
!!created 11 Dec 2006
!!revised 29 Dec 2006
!!descr   fast computation of PDF of L3CCD output distribution
!!root    fpdf
!!index   fpdf
!!class   Programs
!!css     style.css
!!head1   fpdf -- fast computation of PDF of L3CCD output distribution
 
Program to generate the PDF of an avalanche register. It can cope with
a specific number of input electrons or a Poisson distribution and
allows for CICs generated within the register (i.e. output electrons 
even in the absence of any input electrons) and readout noise added after
it. Since CICs generated in the avalanche register do not run through the full
set of stages they have a different probability distribution from the
output generated from genuine input electrons. The program makes use of FFT methods
enabled by the special case of constant avalanche probability and is much faster
than alternative methods.

The readout noise is added by first casting the probabilities onto a finely spaced array
(factor of 10 more fine) and then convolving in the readout noise. This reduces digitisation
effects.

The program finally digitises the results to give values that can be compared directly
to those from ultracam's 'hist'. Note that parameters 'nmax' and 'x1' and 'x2' should not
really be independent but allow for cases where x2 has been taken to unrealistically high
values such that an equivalent 'nmax' would be far too large.

!!head2 Invocation

fpdf nstep nmax ffac gain increase pcic nelec poisson read bias cfac nbin x1 x2



!!head2 Arguments

!!table 

!!arg{nstep}{Number of avalanche multiplication steps.}

!!arg{nmax}{The maximum output number to consider in electrons.}  

!!arg{ffac}{FFT convolutions are used which wrap round, forcing arrays of at least 2
times nmax. This parameter allows you to specify the parameter
precisely. The reason for making it larger than 2 is that you can
rarely choose nmax to be anywhere near the true highest possible value
(=2^nstep) and so you have to rely on the decrease of probability at
high output values (goes as exp(-n/gain) for a single electron input
for example where gain is the mean gain that the program will
report). The actual length for the FFTs will be the smallest power of
2 greater than nstep*ffac; at any one time two such arrays in double
precision will be allocated. The best way to determine the best value
is to experiment by successively increasing it in multiples of 2. At the 
most two FFTs are computed.}

!!arg{gain}{Mean gain from start to end of the avalanche register (sets the multiplication probability per stage).}

!!arg{increase}{Factor increase (can be < 1 to give a decrease) in multiplication probability from 
input to output of register. This is to assess
the influence of variable multiplication probabilities on the statistics. One might by default anticipate a higher
probability near the output from where the clocking voltages are driven. The probabilities increase logarithmically
in the direction of the electron flow.
}

!!arg{pcic}{CIC probability per stage,
i.e. the chance that an electron is created spontaneously on any given step (constant).}

!!arg{nelec}{Number of electrons at input. Has two meanings as determined by the next variable.} 

!!arg{poisson}{Assume 'nelec' represents a mean with a Poisson distribution 
or an exact number (nearest integer will be taken).}  

!!arg{read}{read noise in electrons assumed to be added on ouput. 'read' is 
the RMS of this, which is assumed to be gaussian.}  

!!arg{bias}{Bias level, i.e. what zero will become in the final output}

!!arg{cfac}{Factor to convert electrons at the end of the register into counts (ADU)}

!!arg{nbin}{Number of output bins}

!!arg{x1}{Left-edge of first output pixel. Its up to the user to ensure that it overlaps with nmax OK}

!!arg{x2}{Right-edge of last output pixel.}

!!table

!!end

*/

#include <cfloat>
#include <cstdlib>
#include <iostream>
#include "trm_subs.h"
#include "trm_array1d.h"
#include "trm_constants.h"
#include "trm_format.h"
#include "trm_input.h"
#include "trm_l3ccd.h"

int main (int argc, char *argv[]){

  try{

    // Maximum number of RMS for read noise
    const double GFAC = 9.5;

    // Construct Input object
    Subs::Input input(argc, argv, L3ccd::L3CCD_ENV, L3ccd::L3CCD_DIR);

    // Define inputs
    input.sign_in("nstep",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nmax",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ffac",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("gain",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("increase", Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("pcic",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nelec",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("poisson",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("read",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("bias",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("cfac",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nbin",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",       Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x2",       Subs::Input::GLOBAL, Subs::Input::PROMPT);

    int nstep;
    input.get_value("nstep", nstep, 591, 1, 10000, "number of multiplication steps");
    int nmax;
    input.get_value("nmax", nmax, 100000, 1, 100000000, "maximum number of output electrons");
    double ffac;
    input.get_value("ffac", ffac, 2., 2., 1000., "extra length for Fourier convolutions");
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
    std::cout << "# Generated by fpdf with the following inputs:\n";
    std::cout << "#\n";
    std::cout << "# nstep    = " << nstep  << std::endl;
    std::cout << "# nmax     = " << nmax   << std::endl;
    std::cout << "# ffac     = " << form(ffac)    << std::endl;
    std::cout << "# gain     = " << form(gain)    << std::endl;
    std::cout << "# increase = " << form(increase)<< std::endl;
    std::cout << "# pcic     = " << form(pcic)    << std::endl;
    std::cout << "# nelec    = " << nelec   << std::endl;
    std::cout << "# poisson  = " << poisson << std::endl;
    std::cout << "# read     = " << form(read)    << std::endl;
    std::cout << "# bias     = " << form(bias)    << std::endl;
    std::cout << "# cfac     = " << form(cfac)    << std::endl;
    std::cout << "# nbin     = " << nbin        << std::endl;
    std::cout << "# x1       = " << form(x1)    << std::endl;
    std::cout << "# x2       = " << form(x2)    << std::endl;
    std::cout << "#\n";
    std::cout << "# Derived quantities:\n";
    std::cout << "#\n";
    std::cout << "# Start multiplication probability   = " << form(pstart)            << std::endl;
    std::cout << "# End multiplication probability     = " << form(increase*pstart)   << std::endl;
    std::cout << "# Directly measured gain             = " << form(gr/cfac)           << std::endl;
    std::cout << "# Mean for 1 electron input, no CICs = " << form(bias+gr/cfac)      << std::endl;
    std::cout << "# RMS for 1 electron input, no CICs  = " << form(sqrt(vr)/cfac)     << std::endl;
    std::cout << "# Mean for 0 electron input          = " << form(bias+mcic/cfac)     << std::endl;
    std::cout << "# Overall mean                       = " << form(bias + (mcic+gain*mui)/cfac) << std::endl;
    if(poisson)
	std::cout << "# Overall RMS                        = " << form(sqrt(read*read + vcic + (gain*gain + vr)*mui)/cfac)      << std::endl;
    else
	std::cout << "# Overall RMS                        = " << form(sqrt(read*read + vcic + vr*mui)/cfac)      << std::endl;

    std::cout << "#\n";

    // Compute size for FFTs
    const int NFFT    = int(pow(2.,int(log(ffac*(nmax+1))/log(2.))+1));
    std::cerr << "Number of points for the FFTs = " << NFFT << std::endl;

    // Grab all memory now so we know quickly whether it will fail
    double *prob   = new double[NFFT];
    double *prcic  = new double[NFFT];
    
    // Generate FFT equivalent to a single electron input (i.e. with P(1) = 1)
    // Note prob[0] never modified after here because it must always = 1 
    prob[0] =  1;
    prob[1] = -1;
    for(int i=2; i<NFFT; i += 2){
      double theta = Constants::TWOPI*i/2/NFFT;
      prob[i]   = cos(theta);
      prob[i+1] = sin(theta);
    }

    // Initialise the CIC array
    const double CPCIC = 1.-pcic;
    prcic[0] = 1;
    prcic[1] = 1;
    for(int i=2; i<NFFT; i += 2){
      prcic[i]   = 1.;
      prcic[i+1] = 0.;
    }

    double preal, pimag, real, imag;
    for(int ns=1; ns<=nstep; ns++){

	// Accumulate CIC FFT product
	prcic[1] *= (CPCIC + pcic*prob[1]);
	for(int i=2; i<NFFT; i += 2){
	    preal      = CPCIC + pcic*prob[i];
	    pimag      = pcic*prob[i+1]; 
	    real       = preal*prcic[i] - pimag*prcic[i+1];
	    imag       = pimag*prcic[i] + preal*prcic[i+1];
	    prcic[i]   = real;
	    prcic[i+1] = imag;
	}

	// Compute multiplication probability for step ns (note reverse counting is used here)
	double pmult  = exp(log(pstart)  + log(increase)*(nstep-ns)/(nstep-1));
	double cpmult = 1.-pmult;

	// Apply recurrence
	prob[1] = cpmult*prob[1] + pmult*Subs::sqr(prob[1]);
	for(int i=2; i<NFFT; i += 2){
	    real      = prob[i]*prob[i] - prob[i+1]*prob[i+1];
	    imag      = 2*prob[i]*prob[i+1];
	    prob[i]   = cpmult*prob[i]   + pmult*real;
	    prob[i+1] = cpmult*prob[i+1] + pmult*imag;
	}

    }    
    
    // Now account for the number and distribution of input electrons
    if(poisson || nelec == 0){

	// For a mean NE input electrons with a Poisson distribution, the output FFT is the probability
	// generating function for the Poisson distribution applied to the FFT for 1 electron
	prob[1] = exp(nelec*(prob[1]-1));
	for(int i=2; i<NFFT; i += 2){
	    double a     = exp(nelec*(prob[i]-1));
	    double theta = nelec*prob[i+1];
	    prob[i]      = a*cos(theta);
	    prob[i+1]    = a*sin(theta);
	}
	
    }else{
	
	// For exactly NE input electrons, the FFT is the FFT for 1 electron to the power NE
	int ne = int(Subs::nint(nelec));
	prob[1] = pow(prob[1], ne);
	for(int i=2; i<NFFT; i += 2){
	    double a     = pow(sqrt(prob[i]*prob[i] + prob[i+1]*prob[i+1]), ne);
	    double theta = ne*atan2(prob[i+1], prob[i]);
	    prob[i]     = a*cos(theta);
	    prob[i+1]   = a*sin(theta);
	}
	
    }
    
    // Convolve in the CICs
    prob[1] *= prcic[1];
    for(int i=2; i<NFFT; i += 2){
	real      = prob[i]*prcic[i] - prob[i+1]*prcic[i+1];
	imag      = prob[i+1]*prcic[i] + prob[i]*prcic[i+1];
	prob[i]   = real;
	prob[i+1] = imag;
    }
    
    // Free up CIC buffer
    delete[] prcic;

    // We now have the FFT of the final output PDF, so invert and normalise
    Subs::fftr(prob, NFFT, -1);

    for(int i=0; i<NFFT; i++) prob[i] *= (2./NFFT);

    // Array for binned data
    Subs::Array1D<double> ybin(nbin);
    ybin = 0.;

    if(read > 0.){

      // Add read noise if any specified. Do so by casting into a finer
      // array so that read noise can be added without digitisation noise.

      // Factor to make fine array for adding read noise
      const int FFAC = 10;

      // Number for fine array
      const int NFINE  = int(pow(2.,ceil(log(FFAC*(nmax + 1 + 2*GFAC*read))/log(2.))));

      // For DFT convolution
      const int NRFFT  = 2*NFINE; 

      std::cerr << "Number of points for read noise FFTs = " << NRFFT << std::endl;

      // Registration of coarse relative fine array
      const int OFFSET = int((NFINE - FFAC*nmax)/2.);

      // Load up fine array
      double *fine = new double[NRFFT];

      for(int i=0; i<NRFFT; i++) fine[i] = 0;

      for(int i=0; i<nmax+1; i++)
	fine[OFFSET+FFAC*i] = prob[i];
      
      // Clean up dynamically allocated memory
      delete[] prob;

      // Now load up gaussian for convolution
      double *gauss = new double[NRFFT];

      // Generate gaussian read noise RMS = read
      for(int i=0; i<NFFT; i++) gauss[i] = 0.;

      const int NRMAX =  int(FFAC*GFAC*read) + 1 > NFINE ? NFINE : int(FFAC*GFAC*read) + 1;
      double sum = 1.;
      gauss[0] = 1.;
      for(int i=1; i<NRMAX; i++){
	gauss[i] = exp(-Subs::sqr(i/(FFAC*read))/2.);
	sum += 2*gauss[i];
      }
	
      // Normalise to unit probability
      gauss[0] /= sum;
      for(int i=1; i<NRMAX; i++){
	gauss[i] /= sum;
	gauss[NRFFT-i] = gauss[i];
      }
	
      // Carry out convolution
      Subs::fftr(fine,  NRFFT, 1);
      Subs::fftr(gauss, NRFFT, 1);
	
      fine[1] *= gauss[1];
      for(int i=2; i<NRFFT; i += 2){
	double real = fine[i]*gauss[i] - fine[i+1]*gauss[i+1];
	double imag = fine[i+1]*gauss[i] + fine[i]*gauss[i+1];
	fine[i]     = real;
	fine[i+1]   = imag;
      }
      delete[] gauss;

      // Inverse FFT
      Subs::fftr(fine, NRFFT, -1);

      // Digitise

      for(int i=0; i<NFINE; i++){

	// Digitise & bin in one line
	int nb = int(Subs::nint(nbin*(Subs::nint(bias + double(i-OFFSET)/FFAC/cfac)-x1)/(x2-x1)-0.5));
	if(nb >= 0 && nb < nbin)
		ybin[nb] += 2*fine[i]/NRFFT;
      }
      delete[] fine;

    }else{
    
      // No read noise case. Digitise immediately
      for(int i=0; i<nmax+1; i++){

	// Digitise & bin in one line
	int nb = int(Subs::nint(nbin*(Subs::nint(bias + i/cfac)-x1)/(x2-x1)-0.5));
	if(nb >= 0 && nb < nbin) ybin[nb] += prob[i];
      }

      // Clean up dynamically allocated memory
      delete[] prob;

    }

    // Print out the results
    for(int i=0; i<nbin; i++)
	std::cout << form(x1+(x2-x1)*(i+0.5)/nbin) << " " << form(ybin[i]) << std::endl;    

  }
  catch(const std::string& err){
      std::cerr << err << std::endl;
      exit(EXIT_FAILURE);
  }
}

