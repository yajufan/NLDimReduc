//
// File:        nldr_Isomap.hxx
//


#ifndef included_nldr_Isomap
#define included_nldr_Isomap

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


class nldr_Isomap: public nldr_NLDimRed
{
public:

   /** 
   * Default Constructor
   **/
   nldr_Isomap();

   /** 
   * Copy Constructor
   **/
   nldr_Isomap(const nldr_Isomap &rhs);

   /** 
   * Constructor created with an appropriate nearest neighbor
   * algorithm. The nearest neighbor object is created in the calling
   * program and passed as a pointer to the base class.
   **/
   nldr_Isomap(nns_NNSearch<double> *nns_algo);


  /** 
   * Assignment operator
   **/
   nldr_Isomap &operator=(const nldr_Isomap &rhs);

  /** 
   * Destructor
   **/
   virtual ~nldr_Isomap();


   /**
   * Set the method to be used for finding the nearest neighbors. This
   * should be created in the calling program and passed to this class as a
   * pointer to the base nearest-neighbor class.
   *
   * If this object is not set by the user, the default is the 1-nearest
   * neighbor implemented using a naive search. This is set in the constructor.
   **/

   void setNearestNbrAlgo(nns_NNSearch<double> *nns_algo);


   void setNumDimentions(int nDim);


   int getNumDimentions();


   /**
   *  Applies the non-linear dimension reduction model to input, filling
   *  output with result. 
   *
   *   @param input Reference to input instance array 
   *   @param output  Reference to output instance array
   **/
   void apply(const di_InstanceArray &input,
              di_InstanceArray &output);

   /**
   * Compute Isomap results
   * Without computing the residual variances since it takes a long time
   **/
   void apply(const di_InstanceArray &input,
                        di_InstanceArray &output,
						vector<int> &keptIdx,
						vector<double> &eigenValues,
						tbox_DenseMatrix<double> &EigVec,
						tbox_DenseMatrix<double> &GeoDist);
 
   /**
   * Compute Isomap results
   * New coordinates, eigen values, eigen vectors and residual variances. 
   **/
   void apply(const di_InstanceArray &input,
                        di_InstanceArray &output,
						vector<int> &keptIdx,
						vector<double> &residualVariances,
						vector<double> &eigenValues,
						tbox_DenseMatrix<double> &EigVec,
						tbox_DenseMatrix<double> &GeoDist);

   void applySNN(const di_InstanceArray &input_instArray,
                        di_InstanceArray &output,
						vector<int> &keptIdx,
						vector<double> &residualVariances,
						vector<double> &eigenValues,						
						tbox_DenseMatrix<double> &EigVec,
						tbox_DenseMatrix<double> &GeoDist);

private:

   /** 
   * Compute Adjacency Distance Matrix
   * Input: nRow - Number of rows (samples)
   * Output: Distances - Adjacency distance matrix 
   *         contains distances to neighbors
   * The adjacency distance matrix is stored in a STL vector 
   * for row-based matrix access.
   **/
   void adjacencyDist(vector<double> &Distances,
					  const int &nRow);

   /** 
   * Symmetrize Adjacency Matrix
   **/
   void symmetrize(vector<double> &Distances,
				   const int &nRow);


   /**
   * Compute a matrix, each column contains nearest neighbor indices 
   * ordered according to their closeness to that specific sample.
   **/
   void neighborIdxMx(tbox_DenseMatrix<double> &NeighborhoodTable,
								   tbox_DenseMatrix<double> &AdjacencyMx);


   void sharedNNstrength(tbox_DenseMatrix<double> &NeighborhoodTable,
						 tbox_DenseMatrix<double> &AdjacencyMx,
						 tbox_DenseMatrix<double> &SNNstrength);


   /**
   * Floyd's Algorithm
   **/
   void floyd(vector<double> &Distances,
			  const int &nRow);

   /**
   * Remove outliers
   **/
   void removeOutlier(vector<double> &inDistances,
					  const int &nRow,
					  vector<double> &outDistances,
					  int &nRowC,
					  vector<int> &keptIdx);

   /**
   * Out put a double-centered tbox_DenseMatrix
   **/
   void doubleCentering(vector<double> &D_cleaned,
					    int nRowC,
					    tbox_DenseMatrix<double> &m4eigDecom);

   /**
   * Input a matrix (m4eigDecom) in the form of tbox_DenseMatrix
   * which is a column-based matrix
   * used in LAPACK.
   * Output eigen values and left eigen vectors.
   **/
/*   void eigDecom(tbox_DenseMatrix<double> &MxForEigDecom,
	   			 const int &nRow,
				 double *EigVal,
				 double *LeftEigVec);
*/
   /**
   * Compute new coordinates in the form of STL vector 
   * 
   * Note that "eigVal" contains eigenvalues in increasing order
   **/
   void coordinates(int nDim,
					const int nFeature,
					double *EigVal,
					tbox_DenseMatrix<double> &LeftEigVec, 
					vector<double> &NewCoord);

   double min(double value1,
			  double value2);

   double mean(vector<double> &values);

   /**
   * Calculate n-by-n distance matrix
   **/
   void DistMatrix(vector<double> &DataMatrix,
	    			 int nRow,
				  	 int nCol,
				  	 vector<double> &DistMatrix);

   /**
   * Calculate Residual Variance of two variables, A and B.
   **/
   double residualVar(vector<double> &attributeA,
				 	  vector<double> &attributeB);

   /** The algorithm to be used for nearest neighbor search **/
   nns_NNSearch<double> *d_nns_algo;
   int d_nDim;
};


#ifndef LACKS_NAMESPACE
}
#endif

#endif






