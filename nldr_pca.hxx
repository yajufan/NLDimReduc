//
// File:        nldr_pca.hxx
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

#define _DXX_



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






