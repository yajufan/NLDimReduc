//
// File:        nldr_LTSA.hxx
//


#ifndef included_nldr_LTSA
#define included_nldr_LTSA

#ifndef included_Sapphire_config
#include "Sapphire_config.hxx"
#endif

#ifndef included_nldr_Isomap
#include "nldr_Isomap.hxx"
#endif

#ifndef included_nldr_NLDimRed
#include "nldr_NLDimRed.hxx"
#endif

#ifndef included_nns_NNSearch
#include "nns_NNSearch.hxx"
#endif

#ifndef included_di_InstanceArray
#include "di_InstanceArray.hxx"
#endif

#ifndef included_tbox_DenseMatrix
#include "tbox_DenseMatrix.hxx"
#endif

#ifndef included_iostream
#include <iostream>
#define included_iostream
#endif

#ifndef included_string
#include <string>
#define included_string
#endif

#ifndef included_dsyevd
#include <dsyevd.h> 
#define included_dsyevd
#endif

#ifndef included_dgeev
#include <dgeev.h> 
#define included_dgeev
#endif

#define _DXX_


class nldr_LTSA: public nldr_NLDimRed
{
public:

   /** 
   * Default Constructor
   **/
   nldr_LTSA();

   /** 
   * Copy Constructor
   **/
   nldr_LTSA(const nldr_LTSA &rhs);

   /** 
   * Constructor created with an appropriate nearest neighbor
   * algorithm. The nearest neighbor object is created in the calling
   * program and passed as a pointer to the base class.
   **/
   nldr_LTSA(nns_NNSearch<double> *nns_algo);


  /** 
   * Assignment operator
   **/
   nldr_LTSA &operator=(const nldr_LTSA &rhs);

  /** 
   * Destructor
   **/
   virtual ~nldr_LTSA();

   /**
   * Set the method to be used for finding the nearest neighbors. This
   * should be created in the calling program and passed to this class as a
   * pointer to the base nearest-neighbor class.
   *
   * If this object is not set by the user, the default is the 1-nearest
   * neighbor implemented using a naive search. This is set in the constructor.
   **/
   void setNearestNbrAlgo(nns_NNSearch<double> *nns_algo);


   /**
   * Set the number of dimentions used for embedding
   * Default is 1
   **/
   void setNumDimentions(int nDim);


   int getNumDimentions();


   /**
   *  Applies the non-linear dimension reduction model to input, filling
   *  output with result. 
   *
   *   @param input Reference to input instance array 
   *   @param output  Reference to output instance array
   **/
   void apply(const di_InstanceArray &Input,
              di_InstanceArray &Output);


private:


   void LTSA_algorithm(tbox_DenseMatrix<double> &Input,
					  tbox_DenseMatrix<double> &NewCoord);

   /** 
   * Find neighbor's indices
   * using nearest neighbor search.
   **/
   void findNeighborIdx(int centerIdx,
						vector<unsigned int> &neighborIdx);

   int eigvalues(tbox_DenseMatrix<double> &A,
						  tbox_DenseMatrix<double> &Evects,
                          double *Evals);


   /** The algorithm to be used for nearest neighbor search **/
   nns_NNSearch<double> *d_nns_algo;
   int d_nDim;
};


#ifndef LACKS_NAMESPACE
}
#endif

#endif






