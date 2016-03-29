//
// File:        nldr_LLEmbed.hxx
//


#ifndef included_nldr_LLEmbed
#define included_nldr_LLEmbed

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


#define _DXX_




class nldr_LLEmbed: public nldr_NLDimRed
{
public:

   /** 
   * Default Constructor
   **/
   nldr_LLEmbed();

   /** 
   * Copy Constructor
   **/
   nldr_LLEmbed(const nldr_LLEmbed &rhs);

   /** 
   * Constructor created with an appropriate nearest neighbor
   * algorithm. The nearest neighbor object is created in the calling
   * program and passed as a pointer to the base class.
   **/
   nldr_LLEmbed(nns_NNSearch<double> *nns_algo);


  /** 
   * Assignment operator
   **/
   nldr_LLEmbed &operator=(const nldr_LLEmbed &rhs);

  /** 
   * Destructor
   **/
   virtual ~nldr_LLEmbed();

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
   * Set the delta value, parameter used for conditioning a nearly-sigular matrix
   * delta^2 << 1
   **/
   void setDelta(double delta);


   double getDelta();

   void setLSerrorFlag(int LSerror_flag);

   int getLSerrorFlag();

   void setnFeature4DimEstLimits(int nFeature4DimEstLimits);

   int getnFeature4DimEstLimits();

   /**
   *  Applies the non-linear dimension reduction model to input, filling
   *  output with result. 
   *
   *   @param input Reference to input instance array 
   *   @param output  Reference to output instance array
   **/
   void apply(const di_InstanceArray &Input,
              di_InstanceArray &Output);

   void apply(const di_InstanceArray &Input,
              di_InstanceArray &Output,
			  tbox_DenseMatrix<double> &EigVec,
			  vector<double> &eigVal,
			  vector<double> &leastSquareErrors);

   void apply(const di_InstanceArray &Input_instArray,
              di_InstanceArray &Output,
		 	  tbox_DenseMatrix<double> &EigVec,
			  vector<double> &eigenValues,
			  vector<double> &leastSquareErrors,
			  vector<double> &residualVariances,
			  vector<double> &ReconError);


private:

   /**
   * Compute the matrix of optimal reconstruction weights 
   **/
   void reconstrctWeights(tbox_DenseMatrix<double> &Input,
		 				  tbox_DenseMatrix<double> &OmegaMx,
						  vector<double> &leastSquareErrors);
   /**
   * Compute the matrix of optimal reconstruction weights 
   * without computing least square errors
   **/
   void reconstrctWeights(tbox_DenseMatrix<double> &Input,
		 				  tbox_DenseMatrix<double> &OmegaMx);

   /** 
   * Find neighbor's indices
   * using nearest neighbor search.
   **/
   void findNeighborIdx(int centerIdx,
						vector<unsigned int> &neighborIdx);

   /**
   * Given a sample index and its neighbors' indices,
   * solving for its optimal weights on neighbors.
   **/
   void leastSquareProblem(int centerIdx,
		  			       vector<unsigned int> &neighborIdx,
						   tbox_DenseMatrix<double> &Input,
						   vector<double> &omega);

   /**
   * Compute the center that is described by the neighbors and weights
   **/
   void reconstructedCenter(tbox_DenseMatrix<double> &Input,
							vector<unsigned int> &neighborIdx,
							vector<double> &omega,
							vector<double> &reconC);

   /**
   * Compute least square error,
   * A step after finding the optimal weights on neighbors (omega)
   **/
   void leastSquareError(int centerIdx,
						 vector<unsigned int> &neighborIdx,
						 tbox_DenseMatrix<double> &Input,
						 vector<double> &omega,
						 double &LSerror);

  /**
   * Compute dot product of two row vectors s and t
   * of matrix "input".
   **/
/*
   double dotProduct(const di_InstanceArray &Input,
	  			     int s,
				     int t);
*/
   /**
   * Compute dot product of two STL vectors.
   **/
   double dotProduct(vector<double> &A,
					 vector<double> &B);

   void localGram(int centerIdx,
				  tbox_DenseMatrix<double> &Input,
				  vector<unsigned int> &neighborIdx,
				  tbox_DenseMatrix<double> &LocalGramMx);
   /**
   * Compute the neighborhood correlation matrix
   * for data index "centerIdx".
   * Note: The NCMatrix may be nearly singular.
   **/
   void nbrhdCorr(int centerIdx,
				  tbox_DenseMatrix<double> &Input,
				  vector<unsigned int> &neighborIdx,
				  tbox_DenseMatrix<double> &NCMatrix);


   void condition(tbox_DenseMatrix<double> &InputM,
				 tbox_DenseMatrix<double> &OutputM);
/*
   double trace(tbox_DenseMatrix<double> &InputM);
*/
   int inv(tbox_DenseMatrix<double> &InputM,
		   tbox_DenseMatrix<double> &OutputM);

   void lapackErr(int info);


/*   void getCenterMinusNbr(tbox_DenseMatrix<double> &Input,
						  int centerIdx,
						  vector<unsigned int> &neighborIdx,
						  vector<vector<double> > &centerMinusNbr);
*/
   void getCenterMinusNbr(tbox_DenseMatrix<double> &Input,
						  int centerIdx,
						  vector<unsigned int> &neighborIdx,
						  tbox_DenseMatrix<double> &centerMinusNbr);
/*
   void getCenterDotNbr(const di_InstanceArray &Input,
						int centerIdx,
						vector<unsigned int> &neighborIdx,
						vector<double> &centerDotNbr);
*/
   void computeWeights(tbox_DenseMatrix<double> &invNCMatrix,
				 	   vector<double> &omega);
/*
   void computeWeights(const di_InstanceArray &Input,
				 	   int centerIdx,
				 	   vector<unsigned int> &neighborIdx,
				 	   tbox_DenseMatrix<double> &invNCMatrix,
				 	   vector<double> &omega);
*/


   void matrixForEmbedding(tbox_DenseMatrix<double> &OmegaMx,
	  					   tbox_DenseMatrix<double> &EmbdCostMx);


   /**
   * Compute new coordinates 
   **/
   void coordinates(int nDim,
		 		    tbox_DenseMatrix<double> &EigVec, 
					tbox_DenseMatrix<double> &NewCoord);


   /**
   * Compute residual variances of a list of possible dimensions
   **/
   void intrinsicDim(tbox_DenseMatrix<double> &Input,
		 			 tbox_DenseMatrix<double> &OmegaMx,
					 tbox_DenseMatrix<double> &EigVec, 
					 vector<double> &residualVariances);

   /**
   * Use reconstruction weight from Y on X
   * Calculate the least square error
   **/
   void reconstructionError(tbox_DenseMatrix<double> &Input,
						    tbox_DenseMatrix<double> &EigVec,
		 				    vector<double> &ReconError);

   /** The algorithm to be used for nearest neighbor search **/
   nns_NNSearch<double> *d_nns_algo;
   int d_nDim;
   double d_delta; // A value << 1.
   int d_LSerror_flag;
   int d_nFeature4DimEstLimits;
};


#ifndef LACKS_NAMESPACE
}
#endif

#endif






