//
// File:        nldr_pca.hxx
// Package:     Sapphire Non-Linear Dimension Reduction
// Copyright:   
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Non-Linear Dimension Reduction Class for Isomap algorithm
// Authors:     Chandrika Kamath and Ya-Ju Fan, LLNL
//


#ifndef included_nldr_pca
#define included_nldr_pca

#ifndef included_Sapphire_config
#include "Sapphire_config.hxx"
#endif

#ifndef included_nldr_NLDimRed
#include "nldr_NLDimRed.hxx"
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

#ifndef LACKS_NAMESPACE
using namespace std;
namespace Sapphire {
#endif
#define _DXX_


/**
 *
 * Class for the Isomap algorithm for nonlinear dimension reduction. 
 *
 * TODO: Yaru- pls add any references used in the implementation. You will
 * probably need to include other header files, such as the one for
 * eigenvalue decomp. You will also need to add functions to get and set
 * different variables used in the code or to provide the user relevant
 * information e.g. the eigenvalues.
 *
 */

class nldr_pca: public nldr_NLDimRed
{
public:

   /** 
   * Default Constructor
   **/
   nldr_pca();

   /** 
   * Copy Constructor
   **/
   nldr_pca(const nldr_pca &rhs);

  /** 
   * Assignment operator
   **/
   nldr_pca &operator=(const nldr_pca &rhs);

  /** 
   * Destructor
   **/
   virtual ~nldr_pca();

   void setNumDimentions(int nDim);

   int getNumDimentions();

   void setOptDimFlag(int optDimFlag);

   int getOptDimFlag();

   void apply(const di_InstanceArray &Input,
              di_InstanceArray &Output);

   void apply(const di_InstanceArray &Input,
		 	  vector<double> &eigenvalues,
			  tbox_DenseMatrix<double> &EigVec,
              di_InstanceArray &Output);


private:

   /**
   * PCA
   **/
   void pca(tbox_DenseMatrix<double> &X,
			int &initDim,
			vector<double> &eigVal,
			tbox_DenseMatrix<double> &EigVec,
			tbox_DenseMatrix<double> &Y);

   /**
   * Covariance Matrix
   **/
   void covarMx(tbox_DenseMatrix<double> &X,
				tbox_DenseMatrix<double> &CovarM);


   /**
   * Subtract off the mean
   **/
   void zeroMean(tbox_DenseMatrix<double> &Y);


   int d_nDim;
   int d_optDim_flag;

};


#ifndef LACKS_NAMESPACE
}
#endif

#endif





