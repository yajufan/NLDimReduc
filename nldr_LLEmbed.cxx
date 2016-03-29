//
// File:        nldr_LLEmbed.cxx
// Package:     Sapphire Non-Linear Dimension Reduction
// Copyright:   
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Non-Linear Dimension Reduction Class for LLEmbed algorithm
// Authors:     Chandrika Kamath and Ya-Ju Fan, LLNL
//


#ifndef included_nldr_LLEmbed
#include "nldr_LLEmbed.hxx"
#endif

#ifndef included_dgetri
#include <dgetri.h>
#define included_dgetri
#endif

#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>


#ifndef LACKS_NAMESPACE
using namespace std;
namespace Sapphire {
#endif
#define _DXX_

class Timer {
private:

    timeval startTime,endTime;
    long seconds, nseconds;

public:
    void tic(){
        gettimeofday(&startTime, NULL);
    }

    void toc(){
        gettimeofday(&endTime, NULL);
    }

    double duration(){
		double du;
        seconds  = endTime.tv_sec  - startTime.tv_sec;
        nseconds = endTime.tv_usec - startTime.tv_usec;
		du = seconds*1000 + nseconds/1000.0;
        //du = seconds + nseconds/1000000.0;

        return du;
    }
};

/**
 *
 * Class for the LLEmbed algorithm for nonlinear dimension reduction. 
 *
 * TODO: Yaru- pls add any references used in the implementation. You will
 * probably need to include other header files, such as the one for
 * eigenvalue decomp. You will also need to add functions to get and set
 * different variables used in the code or to provide the user relevant
 * information e.g. the eigenvalues.
 *
 */

/** 
* Default Constructor
**/
nldr_LLEmbed::nldr_LLEmbed()
{
	d_nDim = 1;
	d_delta = 0.001;
	d_LSerror_flag = 0;
	d_nFeature4DimEstLimits = 100;
}

/** 
* Copy Constructor
**/
nldr_LLEmbed::nldr_LLEmbed(const nldr_LLEmbed &rhs)
{
	*this = rhs;
}

/** 
* Constructor created with an appropriate nearest neighbor
* algorithm. The nearest neighbor object is created in the calling
* program and passed as a pointer to the base class.
**/
nldr_LLEmbed::nldr_LLEmbed(nns_NNSearch<double> *nns_algo)
{
	d_nDim = 1;
	d_delta = 0.001;
	d_LSerror_flag = 0;
	d_nns_algo = nns_algo;
	d_nFeature4DimEstLimits = 100;
}


/** 
* Assignment operator
**/
nldr_LLEmbed& nldr_LLEmbed::operator=(const nldr_LLEmbed &rhs)
{
	*this = rhs;
}

/** 
* Destructor
**/
nldr_LLEmbed::~nldr_LLEmbed()
{

}

/**
* Set the method to be used for finding the nearest neighbors. This
* should be created in the calling program and passed to this class as a
* pointer to the base nearest-neighbor class.
*
* If this object is not set by the user, the default is the 1-nearest
* neighbor implemented using a naive search. This is set in the constructor.
**/


void nldr_LLEmbed::setNearestNbrAlgo(nns_NNSearch<double> *nns_algo)
{
	d_nns_algo = nns_algo;
}


void nldr_LLEmbed::setNumDimentions(int nDim)
{
	d_nDim = nDim;
}


int nldr_LLEmbed::getNumDimentions()
{
	return d_nDim;
}

void nldr_LLEmbed::setDelta(double delta)
{
	d_delta = delta;
}

double nldr_LLEmbed::getDelta()
{
	return d_delta;
}

void nldr_LLEmbed::setLSerrorFlag(int LSerror_flag)
{
	d_LSerror_flag = LSerror_flag;
}

int nldr_LLEmbed::getLSerrorFlag()
{
	return d_LSerror_flag;
}

void nldr_LLEmbed::setnFeature4DimEstLimits(int nFeature4DimEstLimits)
{
	d_nFeature4DimEstLimits = nFeature4DimEstLimits;
}

int nldr_LLEmbed::getnFeature4DimEstLimits()
{
	return d_nFeature4DimEstLimits;
}

/**
*  Applies the non-linear dimension reduction model to Input, filling
*  Output with result. 
*
*   @param Input Reference to Input instance array 
*   @param Output  Reference to Output instance array
**/
void nldr_LLEmbed::apply(const di_InstanceArray &Input_instArray,
                         di_InstanceArray &Output)
{
	int nRow,nCol;
	tbox_DenseMatrix<double> Input,Input_temp;
	tbox_DenseMatrix<double> OmegaMx;
	tbox_DenseMatrix<double> EmbdCostMx;
	tbox_DenseMatrix<double> NewCoord;
	double *eigVal;
	tbox_DenseMatrix<double> EigVec;
	vector<double> eigenValues;

	InstanceArray2DenseMatrix(Input_instArray,Input);

	nCol = Input.getNcols();
	nRow = Input.getNrows();	

	if(d_nDim < 1)
	{
		d_nDim = nCol;
	}
	else if(d_nDim > nCol)
	{
		d_nDim = nCol;
	}

	cout<<"Number of samples: "<<nRow<<endl;
	cout<<"Number of features: "<<nCol<<endl;

	//-----------------------
	// Solving for weights
	//-----------------------
	cout<<"Solving for weights... "<<endl;
/*	if(d_LSerror_flag == 1)
	{
		reconstrctWeights(Input,OmegaMx,leastSquareErrors);
	}
	else
	{*/
		reconstrctWeights(Input,OmegaMx);	
/*	}*/
	//-----------------------
	// Matrix for embedding
	//-----------------------
	matrixForEmbedding(OmegaMx,EmbdCostMx);

	//----------------------------
	// Eigenvalue Decomposition
	//----------------------------
	cout<<"Eigenvalue Decomposition... "<<endl;
	EigVec.resize(nRow,nRow);
	tbox_DenseMatrix<double> EigVec_temp;
	EigVec_temp.resize(nRow,nRow);
	vector<double> eigVal_temp;

	eigVal = new double[nRow];
	EmbdCostMx.eig(EigVec_temp,eigVal);
	
	array2vector(eigVal,eigVal_temp,nRow);
	delete [] eigVal;

	orderEigVals(EigVec_temp,eigVal_temp,EigVec,eigenValues);

	//------------------------------
	// New coordinates
	//------------------------------
	if(d_nDim > nCol)
	{
		// Give a warning.
		cout<<"The number of dimension required exceeds totol number of dimensions available."<<endl;
		cout<<"Set d_nDim = nCol. "<<endl;
		d_nDim = nCol;
	}

	coordinates(d_nDim,EigVec,NewCoord);

	wrapResults(Input_instArray,NewCoord,Output);
}

void nldr_LLEmbed::apply(const di_InstanceArray &Input_instArray,
                         di_InstanceArray &Output,
						 tbox_DenseMatrix<double> &EigVec,
						 vector<double> &eigenValues,
						 vector<double> &leastSquareErrors)
{
	int nRow,nCol;
	tbox_DenseMatrix<double> Input,Input_temp;
	tbox_DenseMatrix<double> OmegaMx;
	tbox_DenseMatrix<double> EmbdCostMx;
	tbox_DenseMatrix<double> NewCoord;
	double *eigVal;

	InstanceArray2DenseMatrix(Input_instArray,Input);

	nCol = Input.getNcols();
	nRow = Input.getNrows();	

	if(d_nDim < 1)
	{
		d_nDim = nCol;
	}
	else if(d_nDim > nCol)
	{
		d_nDim = nCol;
	}

	cout<<"Number of samples: "<<nRow<<endl;
	cout<<"Number of features: "<<nCol<<endl;

	//-----------------------
	// Solving for weights
	//-----------------------
	cout<<"Solving for weights... "<<endl;
	if(d_LSerror_flag == 1)
	{
		reconstrctWeights(Input,OmegaMx,leastSquareErrors);
	}
	else
	{
		reconstrctWeights(Input,OmegaMx);	
	}
	//-----------------------
	// Matrix for embedding
	//-----------------------
	matrixForEmbedding(OmegaMx,EmbdCostMx);

	//----------------------------
	// Eigenvalue Decomposition
	//----------------------------
	cout<<"Eigenvalue Decomposition... "<<endl;
	EigVec.resize(nRow,nRow);
	tbox_DenseMatrix<double> EigVec_temp;
	EigVec_temp.resize(nRow,nRow);
	vector<double> eigVal_temp;

	eigVal = new double[nRow];
	EmbdCostMx.eig(EigVec_temp,eigVal);
	
	array2vector(eigVal,eigVal_temp,nRow);
	delete [] eigVal;

	orderEigVals(EigVec_temp,eigVal_temp,EigVec,eigenValues);

	//------------------------------
	// New coordinates
	//------------------------------
	if(d_nDim > nCol)
	{
		// Give a warning.
		cout<<"The number of dimension required exceeds totol number of dimensions available."<<endl;
		cout<<"Set d_nDim = nCol. "<<endl;
		d_nDim = nCol;
	}

	coordinates(d_nDim,EigVec,NewCoord);

	wrapResults(Input_instArray,NewCoord,Output);

}

/**
*   @param Input Reference to Input instance array 
*   @param Output  Reference to Output instance array
*	Note: Having reconstruction error (ReconError)
* 		  will take a much longer time for the function
*		  to complete.
**/
void nldr_LLEmbed::apply(const di_InstanceArray &Input_instArray,
                         di_InstanceArray &Output,
						 tbox_DenseMatrix<double> &EigVec,
						 vector<double> &eigenValues,
						 vector<double> &leastSquareErrors,
						 vector<double> &residualVariances,
						 vector<double> &ReconError)
{
	int nRow,nCol;
	tbox_DenseMatrix<double> Input,Input_temp;
	tbox_DenseMatrix<double> OmegaMx;
	tbox_DenseMatrix<double> EmbdCostMx;
	tbox_DenseMatrix<double> NewCoord;
	double *eigVal;

	InstanceArray2DenseMatrix(Input_instArray,Input);

	nCol = Input.getNcols();
	nRow = Input.getNrows();	

	if(d_nDim < 1)
	{
		d_nDim = nCol;
	}
	else if(d_nDim > nCol)
	{
		d_nDim = nCol;
	}

	cout<<"Number of samples: "<<nRow<<endl;
	cout<<"Number of features: "<<nCol<<endl;

	//-----------------------
	// Solving for weights
	//-----------------------
	cout<<"Solving for weights... "<<endl;
	if(d_LSerror_flag == 1)
	{
		reconstrctWeights(Input,OmegaMx,leastSquareErrors);
	}
	else
	{
		reconstrctWeights(Input,OmegaMx);	
	}
	//-----------------------
	// Matrix for embedding
	//-----------------------
	matrixForEmbedding(OmegaMx,EmbdCostMx);

	//----------------------------
	// Eigenvalue Decomposition
	//----------------------------
	cout<<"Eigenvalue Decomposition... "<<endl;
	EigVec.resize(nRow,nRow);
	tbox_DenseMatrix<double> EigVec_temp;
	EigVec_temp.resize(nRow,nRow);
	vector<double> eigVal_temp;

	eigVal = new double[nRow];
	EmbdCostMx.eig(EigVec_temp,eigVal);
	
	array2vector(eigVal,eigVal_temp,nRow);
	delete [] eigVal;

	orderEigVals(EigVec_temp,eigVal_temp,EigVec,eigenValues);

	//------------------------------
	// New coordinates
	//------------------------------
	if(d_nDim > nCol)
	{
		// Give a warning.
		cout<<"The number of dimension required exceeds totol number of dimensions available."<<endl;
		cout<<"Set d_nDim = nCol. "<<endl;
		d_nDim = nCol;
	}
	
/*	cout<<"Residual Variances..."<<endl;
	intrinsicDim(Input,OmegaMx,EigVec,residualVariances);
*/
	coordinates(d_nDim,EigVec,NewCoord);

	cout<<"Reconstruction Error..."<<endl;
	reconstructionError(Input,EigVec,ReconError);

	wrapResults(Input_instArray,NewCoord,Output);
}

/**
* Compute the matrix of optimal reconstruction weights 
**/
void nldr_LLEmbed::reconstrctWeights(tbox_DenseMatrix<double> &Input,
									 tbox_DenseMatrix<double> &OmegaMx,
									 vector<double> &LSerrors)
{
	int i,j,nRow,numNbr;
	vector<unsigned int> neighborIdx;
	vector<double> omega;
	double error;
	
	nRow = Input.getNrows();

	OmegaMx.resize(nRow,nRow);
	// Initialize
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			OmegaMx.setItem(i,j,0.0);
		}
	}

	// Weights on neighbors
	LSerrors.clear();
	for(i=0;i<nRow;i++)
	{
		findNeighborIdx(i,neighborIdx);

		leastSquareProblem(i,neighborIdx,Input,omega);

		leastSquareError(i,neighborIdx,Input,omega,error);
		LSerrors.push_back(error);

		numNbr = neighborIdx.size();
		for(j=0;j<numNbr;j++)
		{
			OmegaMx(i,neighborIdx[j]) = omega[j];
		}
	}
}

/**
* Compute the matrix of optimal reconstruction weights 
**/
void nldr_LLEmbed::reconstrctWeights(tbox_DenseMatrix<double> &Input,
									 tbox_DenseMatrix<double> &OmegaMx)
{
	int i,j,nRow,numNbr;
	vector<unsigned int> neighborIdx;
	vector<double> omega;
	
	//---------------------------
	// Temp
	//---------------------------
	int numberOfNeighbors;
	//---------------------------


	nRow = Input.getNrows();

	OmegaMx.resize(nRow,nRow);
	// Initialize
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			OmegaMx.setItem(i,j,0.0);
		}
	}

	numberOfNeighbors = 0;

	// Weights on neighbors
	for(i=0;i<nRow;i++)
	{

		if(nRow > 500 & (i%500) < 0.1)
		{		
			cout<<"Solving LS Problem for "<<i<< " of "<<nRow<<" components..."<<endl;
		}

		findNeighborIdx(i,neighborIdx);

		leastSquareProblem(i,neighborIdx,Input,omega);

		numNbr = neighborIdx.size();
		for(j=0;j<numNbr;j++)
		{
			OmegaMx(i,neighborIdx[j]) = omega[j];
		}

				numberOfNeighbors += numNbr;
	}

	cout<<"Average number of neighbors:"<<numberOfNeighbors/nRow<<endl;

}

/** 
* Find neighbor's indices
* using nearest neighbor search.
**/
void nldr_LLEmbed::findNeighborIdx(int centerIdx,
							       vector<unsigned int> &neighborIdx)
{
	int i,nbrSize;
	vector<unsigned int> neighborIdxTemp;
	
	neighborIdxTemp.clear();
	neighborIdx.clear();

	// (+1) For using index to search 
	// need to include "centerIdx" itself
    // d_nns_algo->setNumNbrs(d_nns_algo->getNumNbrs()+1);
	// Can't increment it if having a for loop outside this!
	//========================================================

	//----------------------------------------
	// Find neighbors
	//----------------------------------------
	
	if (d_nns_algo->search(centerIdx,neighborIdxTemp))
	{
		nbrSize = neighborIdxTemp.size();
			
/*		if (nbrSize > 1)
		{
			// Take out centerIdx itself from neighbor ID list
			for(i=0;i<nbrSize;i++)
			{
				if(neighborIdxTemp[i] != centerIdx)
				{
					neighborIdx.push_back(neighborIdxTemp[i]);
				}
			}
		}*/
		if (nbrSize > 0)
		{
			for(i=0;i<nbrSize;i++)
			{
				neighborIdx.push_back(neighborIdxTemp[i]);
			}
		}
		else
		{
			cout<<"No neighbors exists for component "<<centerIdx<<"."<<endl;
			exit(-1);
		}
	}
	else
	{
		cout<<"No search of neighbors has been made."<<endl;
		cout<<"Invalid search type."<<endl;
		exit(-1);
	}
}



/**
* Given a sample index and its neighbors' indices,
* solving for its optimal weights on neighbors.
**/
void nldr_LLEmbed::leastSquareProblem(int centerIdx,
									  vector<unsigned int> &neighborIdx,
								 	  tbox_DenseMatrix<double> &Input,
								  	  vector<double> &omega)
{
	int i,numNbr;
	tbox_DenseMatrix<double> LocalCovMx,LocalCovMxCond,invLocalCovMx;
	int j;

	//---------------------------------------
	// Local covariance matrix (Gram matrix)
	//---------------------------------------
	numNbr = neighborIdx.size();
	localGram(centerIdx,Input,neighborIdx,LocalCovMx);


	if(numNbr > d_nDim)
	{
		condition(LocalCovMx,LocalCovMxCond);
		//cout<<"k>D; regularization will be used."<<endl;
	}
	else
	{
		int n_temp,m_temp;
		n_temp = LocalCovMx.getNrows();
		m_temp = LocalCovMx.getNcols();
		LocalCovMxCond.resize(n_temp,m_temp);

		for(j=0;j<m_temp;j++)
		{
			for(i=0;i<n_temp;i++)
			{
				LocalCovMxCond(i,j) = LocalCovMx(i,j);
			}
		}
	}
	
	//--------------------------------------------------
	// Apply matrix inverse function in tbox_DenseMatrix
	//--------------------------------------------------
	//LocalCovMxCond.inverse(invLocalCovMx);	

	//----------------------------------------------
	// Apply matrix inverse function in nldr_LLEmbed
	//----------------------------------------------
	//cout<<"Inverting..."<<endl;
	//timer.tic();
	if(inv(LocalCovMxCond,invLocalCovMx) != 0)
	{
		cout<<"Failed to Invert."<<endl;
	}
	//timer.toc();
	//cout << timer.duration() << endl;

/*
	cout<<"Inverted matrix: "<<endl;
	for(int i=0;i<numNbr;i++)
	{
		for(int j=0;j<numNbr;j++)
		{
			cout<<invLocalCovMx.getItem(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/	
	//cout<<"Computing weights for "<<centerIdx<<endl;
	//timer.tic();
	computeWeights(invLocalCovMx,omega);
	//timer.toc();
	//cout << timer.duration() << endl;

/*
	cout<<"Original center: "<<endl;
	for(int j=0;j<nCol;j++)
	{
		cout<<Input[centerIdx][j]<<' ';
	}
	cout<<endl;
*/
}

/**
* Compute the center that is described by the neighbors and weights
**/
void nldr_LLEmbed::reconstructedCenter(tbox_DenseMatrix<double> &Input,
									   vector<unsigned int> &neighborIdx,
									   vector<double> &omega,
									   vector<double> &reconC)
{	
	double temp;
	int i,j,nCol,numNbr;	

	reconC.clear();
	nCol = Input.getNcols();

	numNbr = neighborIdx.size();

	for(j=0;j<nCol;j++)
	{
		temp = 0.0;
		for(i=0;i<numNbr;i++)
		{
			temp += Input(neighborIdx[i],j)*omega[i];
		}
		reconC.push_back(temp);
		
	}

}

/**
* Compute least square error,
* A step after finding the optimal weights on neighbors (omega)
**/
void nldr_LLEmbed::leastSquareError(int centerIdx,
									vector<unsigned int> &neighborIdx,
									tbox_DenseMatrix<double> &Input,
									vector<double> &omega,
									double &LSerror)
{
	int j;
	int nCol;
	vector<double> reconC;

	nCol = Input.getNcols();

	reconstructedCenter(Input,neighborIdx,omega,reconC);

	LSerror= 0.0;
	for(j=0;j<nCol;j++)
	{
		LSerror += pow(Input(centerIdx,j) - reconC[j],2);	
	}

}





/**
* Compute dot product of two row vectors s and t
* of matrix "Input".
**/
/*
double nldr_LLEmbed::dotProduct(const di_InstanceArray &Input,
				  				int s,
				  				int t)	
{
	int i,nCol;	
	double temp;

	nCol = Input.getNumFeatures();

	temp = 0;
	for(i=0;i<nCol;i++)
	{
		temp += Input[s][i]*Input[t][i];
	}
	return temp;
}
*/

double nldr_LLEmbed::dotProduct(vector<double> &A,
				  				vector<double> &B)
{
	int i,nCol,nTemp;	
	double temp;

	nCol = A.size();
	nTemp = B.size();
	
	if(nCol != nTemp)
	{
		cout<<"Unable to compute dot product.";
		cout<<"Sizes of two vectors are not equal."<<endl;
		exit(-1);
	}

	temp = 0;
	for(i=0;i<nCol;i++)
	{
		temp += A[i]*B[i];
	}
	return temp;
}


void nldr_LLEmbed::localGram(int centerIdx,
							 tbox_DenseMatrix<double> &Input,
							 vector<unsigned int> &neighborIdx,
							 tbox_DenseMatrix<double> &LocalGramMx)
{
	unsigned int i,j,t,numNbr;
	double temp;
	//vector<vector<double> > centerMinusNbr;
	tbox_DenseMatrix<double> centerMinusNbr;
	
	vector<double> A,B;
	vector<double>::iterator it;
	

	numNbr = neighborIdx.size();

	//cout<<"C - N"<<endl;
	getCenterMinusNbr(Input,centerIdx,neighborIdx,centerMinusNbr);
	//cout<<"C - N done"<<endl;
	unsigned int nCol = Input.getNcols();
/*
	cout<<"centerMinusNbr:"<<endl;
	for(i=0;i<numNbr;i++)
	{
		for(j=0;j<nCol;j++)
		{
			cout<<centerMinusNbr[i][j]<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/

	if(centerMinusNbr.getNrows() != nCol)
	{
		cout<<"Number of columns does not match at localGram()!"<<endl;
	}

	LocalGramMx.resize(numNbr,numNbr);
	for(j=0;j<numNbr;j++)
	{
		//B = centerMinusNbr[j];

		for(i=j;i<numNbr;i++)
		{
			//it = centerMinusNbr.begin()+i*nCol;
			//A = centerMinusNbr[i];
			//it = centerMinusNbr.begin()+j*nCol;
			temp = 0.0;
			//LocalGramMx.setItem(i,j,dotProduct(A,B));
			for(t=0;t<nCol;t++)
			{
				temp = temp + centerMinusNbr(t,j)*centerMinusNbr(t,i);
			}
			LocalGramMx(i,j) = temp;
			LocalGramMx(j,i) = LocalGramMx(i,j);
		}
	}
}


/**
* Compute the local covariance matrix C
* for data index "centerIdx".
* C_ij = z_i^T z_j where z's are neighbors of sample "centerIdx".
* Note: The LocalCovMx may be nearly singular.
**/
/*
void nldr_LLEmbed::localCov(int centerIdx,
							 const di_InstanceArray &Input,
							 vector<unsigned int> &neighborIdx,
							 tbox_DenseMatrix<double> &LocalCovMx)
{
	int i,j,nbrSize;

	cout<<"Data:"<<endl;
	for(i=0;i<Input.getNumInstances();i++)
	{
		for(j=0;j<Input.getNumFeatures();j++)
		{
			cout<<Input[i][j]<<' ';
		}
		cout<<endl;
	}
	cout<<endl;

	// Search for neighbors of "centerIdx"
	findNeighborIdx(centerIdx,neighborIdx);
	nbrSize = neighborIdx.size(); // Size of neighbors

	cout<<"Center:    "<<centerIdx<<endl;
	cout<<"Neighbors: ";
	for(i=0;i<nbrSize;i++)
	{
		cout<<neighborIdx[i]<<' ';
	}
	cout<<endl;

	// Compute neighborhood correlation matrix
	LocalCovMx.resize(nbrSize,nbrSize);
	for(i=0;i<nbrSize;i++)
	{
		for(j=0;j<nbrSize;j++)
		{
			//cout<<neighborIdx[i]<<","<<neighborIdx[j]<<endl;
			LocalCovMx.setItem(i,j,dotProduct(Input,neighborIdx[i],neighborIdx[j]));
			cout<< LocalCovMx.getItem(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;

}
*/
/**
* Conditioning a nearly singular matrix.
* Add a small multiple of the identity matrix to "InputM"
**/
void nldr_LLEmbed::condition(tbox_DenseMatrix<double> &InputM,
					        tbox_DenseMatrix<double> &OutputM)
{
	int i,n;
	double add;
	

	n = InputM.getNcols();
	OutputM = InputM;

	add = d_delta*InputM.trace()/n;

	//cout<<"Add = "<<add<<endl;

	for(i=0;i<n;i++)
	{
		OutputM.setItem(i,i,InputM.getItem(i,i)+add);
	}
}

/*
double nldr_LLEmbed::trace(tbox_DenseMatrix<double> &InputM)
{
	int i,n,m;
	double temp;

	n = InputM.getNrows();
	m = InputM.getNcols();
		
	if(n!=m)
	{
		cout<<"Unable to compute trace. ";
		cout<<"Not a square matrix."<<endl;
		exit(-1);
	}

	temp = 0.0;
	for(i=0;i<n;i++)
	{
		temp += InputM.getItem(i,i);
	}

	return temp;
}
*/

/**
* Invert a matrix using LAPACK function.
* Will return a value to see if the inverse is computed successfully.
**/
int nldr_LLEmbed::inv(tbox_DenseMatrix<double> &InputM,
					  tbox_DenseMatrix<double> &OutputM)
{
	int i,j;
	int n = InputM.getNcols();


	int LDA = n;
	int *ipiv = new int[n];
	double *work = new double[n*n];;
	int lwork = n*n;
	int info;
//	int err = 0;

	double *A = new double[n*n];

	OutputM.resize(n,n);

	// Assign matrix values to A
	for(j=0;j<n;j++)
	{
		for(i=0;i<n;i++)
		{
			A[i*n+j] = InputM.getItem(i,j);
		}
	}


	dgetrf_(n,n,A,LDA,ipiv,info);

	if(info != 0)
	{
		lapackErr(info);
		goto bailout;
	}

	dgetri_(n,A,LDA,ipiv,work,lwork,info);

	if(info != 0)
	{
		lapackErr(info);
		goto bailout;
	}

	// Assign inverted matix values	to Output
	//OutputM.resize(n,n);
	for(j=0;j<n;j++)
	{
		for(i=0;i<n;i++)
		{
			OutputM.setItem(i,j,A[j*n+i]);
		}
	}
	delete [] ipiv;
	delete [] work;
	delete [] A;

	return 0;

bailout:
	delete [] ipiv;
	delete [] work;
	delete [] A;

	return 1;
}


void nldr_LLEmbed::lapackErr(int info)
{
	cout<<"info = "<<info<<endl;
	if(info > 0)
	{
		cout<<"U("<<info<<","<<info<<") is exactly zero."<<endl;
		cout<<"The matrix is singular and its inverse could not be computed."<<endl;
	}
	else
	{
		cout<<"The "<<-info<<"-th argument had an illegal value."<<endl;
	}
}
/*
void nldr_LLEmbed::getCenterMinusNbr(tbox_DenseMatrix<double> &Input,
							  	     int centerIdx,
							  	     vector<unsigned int> &neighborIdx,
								     vector<vector<double> > &centerMinusNbr)
{

	int k,j,nCol,numNbr;

	nCol = Input.getNcols();
	numNbr = neighborIdx.size();

	centerMinusNbr.resize(numNbr);
	for(k=0;k<numNbr;k++)
	{
		for(j=0;j<nCol;j++)
		{
			centerMinusNbr[k].push_back(Input(centerIdx,j)-Input(neighborIdx[k],j));
		}
	}
}
*/
void nldr_LLEmbed::getCenterMinusNbr(tbox_DenseMatrix<double> &Input,
							  	     int centerIdx,
							  	     vector<unsigned int> &neighborIdx,
								     tbox_DenseMatrix<double> &centerMinusNbr)
{

	int k,j,nCol,numNbr;

	nCol = Input.getNcols();
	numNbr = neighborIdx.size();

	centerMinusNbr.resize(nCol,numNbr);
	for(k=0;k<numNbr;k++)
	{
		for(j=0;j<nCol;j++)
		{
			centerMinusNbr(j,k) = Input(centerIdx,j)-Input(neighborIdx[k],j);
		}
	}
}


/*
void nldr_LLEmbed::getCenterDotNbr(const di_InstanceArray &Input,
							  	   int centerIdx,
							  	   vector<unsigned int> &neighborIdx,
								   vector<double> &centerDotNbr)
{
	int i,numNbr;

	numNbr = neighborIdx.size();
	centerDotNbr.clear();	
	for(i=0;i<numNbr;i++)
	{
		centerDotNbr.push_back(dotProduct(Input,centerIdx,neighborIdx[i]));		
	}
}
*/

void nldr_LLEmbed::computeWeights(tbox_DenseMatrix<double> &invLocalCovMx,
								  vector<double> &omega)
{
	int i,j,numNbr;
	double sumRow,sumAll;

	numNbr = invLocalCovMx.getNrows();
	

	sumAll = 0.0;
	for(j=0;j<numNbr;j++)
	{
		for(i=0;i<numNbr;i++)
		{
			sumAll += invLocalCovMx.getItem(i,j);
		}
	}

	omega.clear();
	for(i=0;i<numNbr;i++)
	{
		sumRow = 0.0;
		for(j=0;j<numNbr;j++)
		{
			sumRow += invLocalCovMx.getItem(i,j);
		}	
		omega.push_back(sumRow/sumAll);
	}

	double sumO;
	sumO = 0.0;

	for(i=0;i<numNbr;i++)
	{
		sumO += omega[i];
	}
	//cout<<"Sum of omega = "<<sumO<<endl;
}
/*
void nldr_LLEmbed::computeWeights(const di_InstanceArray &Input,
							  int centerIdx,
							  vector<unsigned int> &neighborIdx,
							  tbox_DenseMatrix<double> &invLocalCovMx,
							  vector<double> &omega)
{
	int i,j,numNbr;
	double tempO;
	double lambda;
	double alpha,beta;
	vector<double> centerDotNbr;
	
	numNbr = neighborIdx.size();

	getCenterDotNbr(Input,centerIdx,neighborIdx,centerDotNbr);
	
	//-------------------
	// Compute lambda
	//-------------------
	alpha = 0.0;
	beta= 0.0;
	for(j=0;j<numNbr;j++)
	{
		for(i=0;i<numNbr;i++)
		{		
			alpha += invLocalCovMx.getItem(i,j)*centerDotNbr[j];
			beta += invLocalCovMx.getItem(i,j);
		}
	}

	cout<<"Center Dot Nbr: "<<endl;
	for(j=0;j<numNbr;j++)
	{
		cout<< centerDotNbr[j]<<' ';
	}
	cout<<endl;

	alpha = 1-alpha;
	lambda = alpha/beta;

	//---------------------
	// Compute omega
	//---------------------
	omega.clear();
	

	for(j=0;j<numNbr;j++)
	{
		tempO = 0.0;
		for(i=0;i<numNbr;i++)
		{
			tempO += invLocalCovMx.getItem(j,i)*(centerDotNbr[i] + lambda);
		}
		omega.push_back(tempO);
	}

	//------------------------------
	// Reconstruct the center point
	//------------------------------

	double sumOmega = 0.0;
	for(i=0;i<omega.size();i++)
	{
		sumOmega += omega[i];
	}
	cout<<"Sum of omega: "<< sumOmega <<endl;

	vector<double> Y;
	double tempY;
	Y.clear();
	
	for(j=0;j<Input.getNumFeatures();j++)
	{
		tempY = 0.0;
		for(i=0;i<numNbr;i++)
		{
			tempY += omega[i]*Input[neighborIdx[i]][j];
		}
		Y.push_back(tempY);
	}
	for(j=0;j<Input.getNumFeatures();j++)
	{
		cout<<Y[j]<<' ';
	}
	cout<<endl;

	for(j=0;j<Input.getNumFeatures();j++)
	{
		cout<<Input[centerIdx][j]<<' ';
	}
	cout<<endl;
}
*/

void nldr_LLEmbed::matrixForEmbedding(tbox_DenseMatrix<double> &OmegaMx,
									  tbox_DenseMatrix<double> &EmbdCostMx)
{
	int i,j;
	int nRow = OmegaMx.getNrows();

	tbox_DenseMatrix<double> EmbdCostMxTemp;


	EmbdCostMxTemp.resize(nRow,nRow);
	EmbdCostMx.resize(nRow,nRow);
	for(i=0;i<nRow;i++)
	{	
		for(j=0;j<nRow;j++)
		{
			if(i==j)
			{
				EmbdCostMxTemp.setItem(i,j,1-OmegaMx.getItem(i,j));
			}
			else
			{
				EmbdCostMxTemp.setItem(i,j,-OmegaMx.getItem(i,j));
			}
		}
	}
	EmbdCostMxTemp.multMatrix('T',EmbdCostMxTemp,'N',EmbdCostMx);
}

/**
* Compute new coordinates 
**/
void nldr_LLEmbed::coordinates(int nDim,
							   tbox_DenseMatrix<double> &EigVec, 
							   tbox_DenseMatrix<double> &NewCoord)
{
	int i,j,nRow;

	nRow = EigVec.getNrows();
	
	NewCoord.resize(nRow,nDim);

	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nDim;j++)
		{
			NewCoord.setItem(i,j,EigVec(i,j+1));
		}	
	}
}

/**
* Compute residual variances for a list of possible dimensions
**/
void nldr_LLEmbed::intrinsicDim(tbox_DenseMatrix<double> &Input,
								tbox_DenseMatrix<double> &OmegaMx,
								tbox_DenseMatrix<double> &EigVec, 
						 	    vector<double> &residualVariances)
{
	int d,nCol,nFeatures;
	vector<double> OmegaMx_vector, NewOmegaMx_vector; // Distance matrix at projected space
	tbox_DenseMatrix<double> NewCoord,NewOmegaMx;

	nCol = Input.getNcols();

	if(nCol > d_nFeature4DimEstLimits)
	{
		nFeatures = d_nFeature4DimEstLimits;
	}
	else
	{
		nFeatures = nCol;
	}

	residualVariances.clear();
	for(d=1;d<=nFeatures;d++)
	{
		coordinates(d,EigVec,NewCoord);
		
		// Compute new weight matrix
		reconstrctWeights(NewCoord,NewOmegaMx);

		// Compare new weight matrix with the old one
		tboxDM2vector(OmegaMx,OmegaMx_vector);
		tboxDM2vector(NewOmegaMx,NewOmegaMx_vector);
		residualVariances.push_back(residualVar(OmegaMx_vector,NewOmegaMx_vector));
	}
}

/**
* Use reconstruction weight from Y on X
* Calculate the least square error
**/
void nldr_LLEmbed::reconstructionError(tbox_DenseMatrix<double> &Input,
									   tbox_DenseMatrix<double> &EigVec,
									   vector<double> &ReconError)
{
	int iter,i,j;
	int nRow,nCol;
	tbox_DenseMatrix<double> reconC;
	double LSerror,errorTemp;


	int d,nFeatures;
	vector<double> OmegaMx_vector, NewOmegaMx_vector; // Distance matrix at projected space
	tbox_DenseMatrix<double> NewCoord,NewOmegaMx;
	tbox_DenseMatrix<double> DNewCoord;

	nRow = Input.getNrows();
	nCol = Input.getNcols();

	// Compute Weight Metrix from Y
	if(nCol > d_nFeature4DimEstLimits)
	{
		nFeatures = d_nFeature4DimEstLimits;
	}
	else
	{
		nFeatures = nCol;
	}

	coordinates(nFeatures,EigVec,NewCoord);

	ReconError.clear();
	for(d=1;d<=nFeatures;d++)
	{

		// Resize to d-dimentional Y Coordinates
		
		DNewCoord.resize(nRow,d);
		for(j=0;j<d;j++)
		{
			for(i=0;i<nRow;i++)
			{
				DNewCoord(i,j) = NewCoord(i,j);
			}
		}
		
		// Compute new weight matrix
		//cout<<"Recon weights..."<<endl;
		reconstrctWeights(DNewCoord,NewOmegaMx);
		
/*		cout<<"d = "<<d<<endl;
		cout<<"New coordinates:"<<endl;
		for(i=0;i<nRow;i++)
		{
			for(j=0;j<d;j++)
			{
				cout<<DNewCoord(i,j)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;

		cout<<"New weight matrix:"<<endl;
		for(i=0;i<nRow;i++)
		{
			for(j=0;j<nRow;j++)
			{
				cout<<NewOmegaMx(i,j)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;
*/
		// Reconstructed Centers
		reconC.resize(nRow,nCol);
		for(iter=0;iter<nRow;iter++)
		{
			for(j=0;j<nCol;j++)
			{
				reconC(iter,j) = 0.0;
				for(i=0;i<nRow;i++)
				{
					reconC(iter,j) += Input(i,j)*NewOmegaMx(iter,i);
				}
			}
		}

/*		cout<<"Reconstructed Centers: "<<endl;
		for(i=0;i<nRow;i++)
		{
			for(j=0;j<nCol;j++)
			{
				cout<<reconC(i,j)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;
*/

		// Sum up reconstruction errors
		errorTemp = 0.0;
		
		for(iter=0;iter<nRow;iter++)
		{
			LSerror= 0.0;
			for(j=0;j<nCol;j++)
			{
				LSerror += pow(Input(iter,j) - reconC(iter,j),2);	
			}
			errorTemp += LSerror;
		}
		double error;
		//cout<<'d = '<<d<<', errorTemp= '<<errorTemp<<endl;
		error = sqrt(errorTemp);

		//cout<<"Reconstruction error: "<<error<<endl;
		//cout<<"Reconstruction error: "<<errorTemp<<endl;

		ReconError.push_back(error);
		//ReconError.push_back(errorTemp);
	}
}

#ifndef LACKS_NAMESPACE
}
#endif






