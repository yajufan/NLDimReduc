//
// File:        nldr_LapEigMap.hxx
//


#ifndef included_nldr_LapEigMap
#define included_nldr_LapEigMap

#ifndef included_Sapphire_config
#include "Sapphire_config.hxx"
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


class nldr_LapEigMap: public nldr_NLDimRed
{
public:

   /** 
   * Default Constructor
   **/
   nldr_LapEigMap();

   /** 
   * Copy Constructor
   **/
   nldr_LapEigMap(const nldr_LapEigMap &rhs);

   /** 
   * Constructor created with an appropriate nearest neighbor
   * algorithm. The nearest neighbor object is created in the calling
   * program and passed as a pointer to the base class.
   **/
   nldr_LapEigMap(nns_NNSearch<double> *nns_algo);


  /** 
   * Assignment operator
   **/
   nldr_LapEigMap &operator=(const nldr_LapEigMap &rhs);

  /** 
   * Destructor
   **/
   virtual ~nldr_LapEigMap();

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

   /**
   * Get the number of dimentions used for embedding
   **/
   int getNumDimentions();

   /**
   * Set the parameter sigma in Heat Kernel for computing weights
   * Default is 2.0.
   **/
   void setHeatKernelParam(double sigma);

   double getHeatKernelParam();

   /**
   * Set the type of weighting methods
   * 1: simple (considering only the connection to neighbors 0/1)
   * 2: Heat kernel
   * Default is 1.
   **/
   void setWeightType(int weightType);

   int getWeightType();


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
			  vector<double> &objectiveValues,
			  int &estDim);

 /*  void applySNN(const di_InstanceArray &Input,
              di_InstanceArray &Output,
			  tbox_DenseMatrix<double> &EigVec,
			  vector<double> &eigVal,
			  vector<double> &objectiveValues,
			  int &estDim);
*/

private:

   /** 
   * Compute Adjacency Matrix
   **/
   void adjacency(tbox_DenseMatrix<double> &AdjacencyMx);


   /**
   * Make sure the adjacency matrix is symmetric
   * Need this especially when using k-nearest-neighbors
   * If component s is among neighbors of t or
   * t is among neighbors of s,  
   * s and t are connected.
   **/
   void symmetrize(tbox_DenseMatrix<double> &AdjacencyMx);

   /**
   * Compute the matrix of optimal reconstruction weights 
   **/
   void weights(tbox_DenseMatrix<double> &Input,
				tbox_DenseMatrix<double> &AdjacencyMx,
				tbox_DenseMatrix<double> &WeightMx);

   void neighborIdxMx(tbox_DenseMatrix<double> &NeighborhoodTable,
					  tbox_DenseMatrix<double> &AdjacencyMx);

/*   void sharedNNstrength(tbox_DenseMatrix<double> &NeighborhoodTable,
						 tbox_DenseMatrix<double> &AdjacencyMx,
						 tbox_DenseMatrix<double> &SNNstrength);
*/
   double heatKernel(int s, 
			   	     int t,
					 tbox_DenseMatrix<double> &Input);


   /**
   * Compute the diagonal weight matrix
   **/
   void diagWeights(tbox_DenseMatrix<double> &WeightMx,
		  	 	    tbox_DenseMatrix<double> &DiagWMx);


   /**
   * Compute the Laplacian matrix
   **/
   void laplacianMatrix(tbox_DenseMatrix<double> &WeightMx,
						tbox_DenseMatrix<double> &DiagWMx,
						tbox_DenseMatrix<double> &LapMx);

   /**
   * Solve the generalized eigenvalue problem
   * Obtaining the right eigenvector v(j):
   * A*v(j) = lambada(j)*B*v(j)
   * For Laplacian Eigenmaps,
   * A = LapMx
   * B = DiagWMx.
   * eigVal_i = -1 if unable to calculate the i-th eigenvalue.
   **/
   int genEigProb(tbox_DenseMatrix<double> &LapMx,
				  tbox_DenseMatrix<double> &DiagWMx,
				  tbox_DenseMatrix<double> &EigVec,
				  vector<double> &eigVal);


/* !!!!!!!!!!!! Should delete this function!!!!!!! */
   /**
   * Order eigen values and their corresponding eigen vectors.
   **/
/*   void orderEigVals(tbox_DenseMatrix<double> &EigVec_temp,
				     vector<double> &eigVal_temp,
				     tbox_DenseMatrix<double> &EigVec,
				     vector<double> &eigVal);

*/

   /**
   * Compute new coordinates 
   **/
   void coordinates(int nDim,
					vector<double> &eigenValues,
					tbox_DenseMatrix<double> &EigVec, 
					tbox_DenseMatrix<double> &NewCoord);

/*   void coordinates(int nDim,
					tbox_DenseMatrix<double> &EigVec, 
					tbox_DenseMatrix<double> &NewCoord);
*/
   /**
   * Compute objective values for a list of possible dimensions
   **/
   void intrinsicDim(int &nFeatures,
					 tbox_DenseMatrix<double> &WeightMx,
					 tbox_DenseMatrix<double> &EigVec, 
					 vector<double> &objectiveValues,
					 int &estDim);
   /**
   * Compute the averaged value of the objective function
   **/
   double objective(tbox_DenseMatrix<double> &WeightMx,
			    	tbox_DenseMatrix<double> &NewCoord);

   int minIndex(vector<double> &obj);

   /** 
   * The algorithm to be used for nearest neighbor search 
   **/
   nns_NNSearch<double> *d_nns_algo;

   unsigned int d_nDim;
   
   /** 
   * Heat kernel parameter, default is 1.0.
   **/
   double d_sigma;

   /**
   * Method used to weight neighbors - 1: simple 2: Heat kernel.
   **/ 
   int d_weightType; 
};


#ifndef LACKS_NAMESPACE
}
#endif

#endif






