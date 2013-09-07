#ifndef TRM_L3CCD
#define TRM_L3CCD

#include <string>
using namespace std;

//! L3CCD software namespace

/*
 * This is the namespace for all the \c l3ccd software mainly to do with
 * simulation of their statistics.
 */

namespace L3ccd {

  //! Default directory for command defaults
  const char L3CCD_DIR[] = ".l3ccd";

  //! Environment variable for switching directory for command defaults
  const char L3CCD_ENV[] = "L3CCD_ENV";

  //! An exception class.

  /** L3ccd::L3ccd_Error is the error class for the L3ccd programs.
   * It is inherited from the standard string class.
   */
  class L3ccd_Error : public std::string {
  public:

    //! Default constructor
    L3ccd_Error() : std::string() {}

    //! Constructor storing a message
    L3ccd_Error(const std::string& err) : std::string(err) {} 
  };

  //! Computes start probability
  double start_prob(int nstep, double increase, double gain);

}

#endif
