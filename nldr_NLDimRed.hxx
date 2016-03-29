//
// File:        nldr_NLDimRed.hxx
// Package:     Sapphire Non-Linear Dimension Reduction
// Copyright:   
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Non-Linear Dimension Reduction Base Class
// Authors:     Chandrika Kamath and Ya-Ju Fan, LLNL
// 	Copied wrapResults() from dr_Transform 
//		- written by Erick Cantu-Paz and Cyrus Harrison
//


#ifndef included_nldr_NLDimRed
#define included_nldr_NLDimRed

#ifndef included_Sapphire_config
#include "Sapphire_config.hxx"
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
 * Abstract base class for nonlinear dimension reduction algorithms. 
 *
 * These techniques cannot easily embed a new data point into the reduced
 * dimensional space created using the original dataset. Therefore, the
 * build and apply approach used in the dr\_DimRed class and its derived
 * classes cannot be used here. Hence, the use of a separate set of
 * classes.
 *
 */

class nldr_NLDimRed
{
public:

   /**
      Applies the non-linear dimension reduction model to input, filling output
      with result. Subclasses implement this method with specialized algorithms.

      @param input  Reference to input instance array
      @param output Reference to output instance array
    */
   virtual void apply(const di_InstanceArray &input, 
                      di_InstanceArray &output)=0;


protected:

   /**
   * Transform InstanceArray to DenseMatrix
   **/
   void InstanceArray2DenseMatrix(const di_InstanceArray &Input,
								  tbox_DenseMatrix<double> &InputDMx);

   double residualVar(vector<double> &X,
				    vector<double> &Y);

   /**
      Wraps a matrix result, into an instance array, also adding and ignored
      features from the input array.
    */
   void wrapResults(const di_InstanceArray &input,
                    tbox_DenseMatrix<double> &mat,
                    di_InstanceArray &output);

   void orderEigVals(tbox_DenseMatrix<double> &EigVec_temp,
					 vector<double> &eigVal_temp,
				  	 tbox_DenseMatrix<double> &EigVec,
				  	 vector<double> &eigVal);


   void vector2tboxDM(vector<double> &rowBasedMatrix,
	     			  tbox_DenseMatrix<double> &columnBasedMatrix,
					  int nRow,
					  int nCol);

   void array2tboxDM(double *rowBasedMatrix,
	     			  tbox_DenseMatrix<double> &columnBasedMatrix,
					  int nRow,
					  int nCol);

   void tboxDM2vector(tbox_DenseMatrix<double> &columnBasedMatrix,
					  vector<double> &output_vector);

   void array2vector(double *input_array,
	     			 vector<double> &output_vector,
					 int nRow);

   double min(double value1,
			  double value2);

   int min(int value1,
		   int value2);

   double max(double value1,
			  double value2);

   double mean(vector<double> &values);

};


#ifndef LACKS_NAMESPACE
}
#endif

#endif






