//
// File:        nldr_LTSA.cxx
// Package:     Sapphire Non-Linear Dimension Reduction
// Copyright:   
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Non-Linear Dimension Reduction Class for LTSA algorithm
// Authors:     Chandrika Kamath and Ya-Ju Fan, LLNL
//


#ifndef included_nldr_LTSA
#include "nldr_LTSA.hxx"
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
 * Class for the LTSA algorithm for nonlinear dimension reduction. 
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
nldr_LTSA::nldr_LTSA()
{
	d_nDim = 1;
}

/** 
* Copy Constructor
**/
nldr_LTSA::nldr_LTSA(const nldr_LTSA &rhs)
{
	*this = rhs;
}

/** 
* Constructor created with an appropriate nearest neighbor
* algorithm. The nearest neighbor object is created in the calling
* program and passed as a pointer to the base class.
**/
nldr_LTSA::nldr_LTSA(nns_NNSearch<double> *nns_algo)
{
	d_nDim = 1;
	d_nns_algo = nns_algo;
}


/** 
* Assignment operator
**/
nldr_LTSA &nldr_LTSA::operator=(const nldr_LTSA &rhs)
{
	*this = rhs;
}

/** 
* Destructor
**/
nldr_LTSA::~nldr_LTSA()
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


void nldr_LTSA::setNearestNbrAlgo(nns_NNSearch<double> *nns_algo)
{
	d_nns_algo = nns_algo;
}


void nldr_LTSA::setNumDimentions(int nDim)
{
	d_nDim = nDim;
}


int nldr_LTSA::getNumDimentions()
{
	return d_nDim;
}


/**
*  Applies the non-linear dimension reduction model to Input, filling
*  Output with result. 
*
*   @param Input Reference to Input instance array 
*   @param Output  Reference to Output instance array
**/
void nldr_LTSA::apply(const di_InstanceArray &Input_instArray,
                         di_InstanceArray &Output)
{
	int nRow,nCol;
	tbox_DenseMatrix<double> Input;
	tbox_DenseMatrix<double> NewCoord;

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

	LTSA_algorithm(Input,NewCoord);


	cout<<"Wraping results..."<<endl;
	wrapResults(Input_instArray,NewCoord,Output);
}

/**
* Compute the correlation matrix
* Perform singular value decomposition
* Obtain the eigenvalue matrix
**/
void nldr_LTSA::LTSA_algorithm(tbox_DenseMatrix<double> &Input,
								  tbox_DenseMatrix<double> &NewCoord)
{
	int i,j,k,kk,nRow,nCol;
	unsigned int numNbr;
	int searchType;
	vector<unsigned int> neighborIdx;
	vector<double> aveNbr;
	double tempSum;
	tbox_DenseMatrix<double> XI,CorMx_temp,CorMx,EmbedMx;

	
	nRow = Input.getNrows();
	nCol = Input.getNcols();
	numNbr = d_nns_algo->getNumNbrs();
	searchType = d_nns_algo->getSearchType();

	if(searchType != 1)
	{
		cout<<"Have to use KNN! Changing search type to KNN..."<<endl;
		d_nns_algo->setSearchType(1);
	}

	if(d_nDim > nCol)
	{
		d_nDim = nCol;
	}

	// Initialize
	XI.resize(numNbr,nCol);

	CorMx_temp.resize(numNbr,numNbr);
	CorMx.resize(numNbr,numNbr);
	for(j=0;j<numNbr;j++)
	{
		for(i=0;i<numNbr;i++)
		{
			CorMx_temp(i,j) = 0.0;
			CorMx(i,j) = 0.0;
		}
	}


	EmbedMx.resize(nRow,nRow);
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			EmbedMx(i,j) = 0.0;
		}
	}

	// to be finished ... //////////////////???

	// Constructing Alignment Matrix B
	tbox_DenseMatrix<double> AlignMxB; 

	AlignMxB.resize(nRow,nRow);
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			AlignMxB(i,j) = 0.0;
		}
	}

	// Weights on neighbors
	for(i=0;i<nRow;i++)
	{
		findNeighborIdx(i,neighborIdx);

		// Average neighbors & centered data
		// numNbr = neighborIdx.size();

/*		cout<<"Nbr indices:"<<endl;
		for(k=0;k<numNbr;k++)
		{
			cout<<neighborIdx[k]<<' ';
		}		
		cout<<endl;
*/
		aveNbr.clear();
		for(j=0;j<nCol;j++)
		{
			tempSum = 0.0;
			for(k=0;k<numNbr;k++)
			{
				tempSum += Input(neighborIdx[k],j);		
			}
			//cout<<tempSum/numNbr<<' ';
			aveNbr.push_back(tempSum/numNbr);
		}
		//cout<<endl;

/*		cout<<"Ave Nbr values:"<<endl;
		for(j=0;j<nCol;j++)
		{
			cout<<aveNbr[j]<<' ';
		}
		cout<<endl;
*/
		for(j=0;j<nCol;j++)
		{
			for(k=0;k<numNbr;k++)
			{
				XI(k,j) = Input(neighborIdx[k],j) - aveNbr[j];
			}
		}

/*		cout<<"XI:"<<endl;
		for(k=0;k<numNbr;k++)
		{
			for(j=0;j<nCol;j++)
			{
				cout<<XI(k,j)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;
*/
		// Correlation matrix
		double cr_temp;
		for(k=0;k<numNbr;k++)
		{
			for(kk=0;kk<numNbr;kk++)
			{
				cr_temp = 0.0;
				for(j=0;j<nCol;j++)
				{
					cr_temp += XI(k,j)*XI(kk,j);
				}
				CorMx_temp(k,kk) = cr_temp;
			}
		}

		// Symmetrize Correlation Matrix
		for(k=0;k<numNbr;k++)
		{
			for(kk=k;kk<numNbr;kk++)
			{
				CorMx(kk,k) = (CorMx_temp(kk,k) + CorMx_temp(k,kk)) / 2.0;
				CorMx(k,kk) = CorMx(kk,k);
			}
		}

/*		cout<<"W:"<<endl;
		for(k=0;k<numNbr;k++)
		{
			for(kk=0;kk<numNbr;kk++)
			{
				cout<<CorMx(k,kk)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;
*/

		// Singular value decomposition

		double *eigVal; // Eigenvalues
		tbox_DenseMatrix<double> EigVec;

		tbox_DenseMatrix<double> EigVec_temp;
		EigVec_temp.resize(numNbr,numNbr);
		vector<double> eigVal_temp,eigenValues;

		EigVec.resize(numNbr,numNbr);
		eigVal = new double[numNbr];

		//CorMx.eig(EigVec_temp,eigVal);

		int info;
	
		info = eigvalues(CorMx,EigVec_temp,eigVal);

		if(info != 0)
		{
			cout<<"info = "<<info<<" at iter "<<i<<endl;
		}

		array2vector(eigVal,eigVal_temp,numNbr);
		delete [] eigVal;

		orderEigVals(EigVec_temp,eigVal_temp,EigVec,eigenValues); // small to large
		
/*		for(k=0;k<numNbr;k++)
		{
			for(j=0;j<d_nDim;j++)
			{
				cout<<EigVec(k,numNbr-j-1)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;
*/

		// Make sure sufficient d
		if (numNbr < d_nDim)
		{
			cout<<"Warning: Target dimensionality reduced to "<<numNbr<<"."<<endl;
			d_nDim = numNbr;
		}

		// Take the d largest eigen vectors
		tbox_DenseMatrix<double> GiVectors;
		GiVectors.resize(d_nDim+1,numNbr);
		double gConst;

		gConst = (double)1.0/sqrt(numNbr);
		for(k=0;k<numNbr;k++)
		{
			GiVectors(0,k) = gConst;
			for(j=1;j<d_nDim+1;j++)
			{
				GiVectors(j,k) = EigVec(k,numNbr-j); // Need the largest d vectors
			}
		}
		
/*		if(i<5)
		{
			cout<<"Gi:"<<endl;
			for(j=0;j<3;j++)
			{
				for(k=0;k<5;k++)
				{
					cout<<GiVectors(j,k)<<' ';
				}
				cout<<endl;
			}
			cout<<endl;
		}
*/
		// Increment the alignment matrix B
		tbox_DenseMatrix<double> G_temp;
		double gMulti;

		G_temp.resize(numNbr,numNbr);
		for(k=0;k<numNbr;k++)
		{
			for(kk=k;kk<numNbr;kk++)
			{
				gMulti = 0.0;
				for(j=0;j<d_nDim+1;j++)
				{
					gMulti += GiVectors(j,k)*GiVectors(j,kk);
				}
				if(k == kk)
				{
					G_temp(k,kk) = 1.0 - gMulti;
				}
				else
				{ 
				    G_temp(k,kk) = -gMulti;		
					G_temp(kk,k) = -gMulti;		
				}

				//cout<<gMulti<<' ';
			}
			//cout<<endl;
		}
		//cout<<endl;

/*		cout<<"G:"<<endl;
		for(k=0;k<numNbr;k++)
		{
			for(kk=0;kk<numNbr;kk++)
			{
				cout<<G_temp(k,kk)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;
*/
		// Use the neighbor indices for alignment matrix B
		
		for(k=0;k<numNbr;k++)
		{
			for(kk=0;kk<numNbr;kk++)
			{
				AlignMxB(neighborIdx[k],neighborIdx[kk]) += G_temp(k,kk);	
			}
		}
	}

	// Symmetrize Alignment Matrix

	tbox_DenseMatrix<double> AlignMx; 

	AlignMx.resize(nRow,nRow);
	for(i=0;i<nRow;i++)
	{
		for(j=i;j<nRow;j++)
		{
			AlignMx(i,j) = (AlignMxB(i,j) + AlignMxB(j,i)) / 2.0;
			AlignMx(j,i) = AlignMx(i,j);
		}
	}

/*	// B(i,i) -= 1.0
	for(i=0;i<nRow;i++)
	{
		AlignMxB(i,i) -= 1.0;
	}
*/

/*	cout<<"B:"<<endl;
	for(i=0;i<5;i++)
	{
		for(j=0;j<5;j++)
		{
			cout<<AlignMx(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
	// Singular value decomposition

	double *eigValB; // Eigenvalues
	tbox_DenseMatrix<double> EigVecB;

	tbox_DenseMatrix<double> EigVecB_temp;
	EigVecB_temp.resize(nRow,nRow);
	vector<double> eigValB_temp,eigenValuesB;
	eigValB_temp.resize(nRow);
	eigenValuesB.resize(nRow);

	EigVecB.resize(nRow,nRow);
	eigValB = new double[nRow];


	cout<<"Final Eigen Decomp..."<<endl;
	//AlignMx.eig(EigVecB_temp,eigValB);


	int info;
	
	info = eigvalues(AlignMx,EigVecB_temp,eigValB);

	if(abs(info) < 1e-8)
	{
		cout<<"Eigen values calculated."<<endl;
	}
	else
	{
		cout<<"Warning:Failed to calculate eigen values."<<endl;
	}

	array2vector(eigValB,eigValB_temp,nRow);

	delete [] eigValB;

	orderEigVals(EigVecB_temp,eigValB_temp,EigVecB,eigenValuesB); // small to large
		
	// Take the 2 to d+1 smallest eigen vectors
	//tbox_DenseMatrix<double> NewCoord;
	cout<<"New Coordinates..."<<endl;
	NewCoord.resize(nRow,d_nDim);


/*	cout<<"Ordered:"<<endl;
	for(i=0;i<8;i++)
	{
		for(j=0;j<8;j++)
		{
			cout<<EigVecB(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
	// Ignore zero eigenvalues
	int start_nonzero;
	start_nonzero = 0;
	while(start_nonzero < nRow)
	{
		if(abs(eigenValuesB[start_nonzero]) < 1e-8)
		{
			start_nonzero += 1;
		}
		else
		{
			break;
		}
	}

	cout<<"nonzero eigen values start at "<<start_nonzero<<endl;
	cout<<"eigenvalue: "<<eigenValuesB[start_nonzero]<<endl;

	for(i=0;i<nRow;i++)
	{
		for(j=0;j<d_nDim;j++)
		{			
			NewCoord(i,j) = EigVecB(i,j+start_nonzero); // Ignore the first vector
		}
	}
}

int nldr_LTSA::eigvalues(tbox_DenseMatrix<double> &A,
						 tbox_DenseMatrix<double> &Evects,
                         double *Evals)
{
   int i,j,n,m;

   n=A.getNrows();
   m=A.getNcols();

 /*  cout<<"Target Matrix for Eigen Evaluation:"<<endl;
   for(i=0;i<5;i++)
   {
      for(j=0;j<5;j++)
      {
         cout<<A(i,j)<<' ';
      }
      cout<<endl;
   }
   cout<<endl;
*/
   if (n != m)
   {
      cerr <<"\nError in eig: matrix needs to be square." << endl;
   }

   int n1=Evects.getNrows();
   int m1=Evects.getNcols();

   if ( (n1!=n) || (m1!=m))
   {
      cerr <<"\nError in eig: dimension of the eigenvector matrix "
           << "has to match dimension of original matrix." << endl;
   }

   double *H = new double[n*n];
   double *E = new double[n];

   char jobz = 'V';         // compute the eigenvectors (V/N)
   char uplo = 'U';         // assume upper triangular is stored (U/L)

	// dgeev --------
	char jobvl = 'V';
	char jobvr = 'V';
	double *wr = new double[n];
	double *wi = new double[n];
	int ldvl = n;
	double *vl = new double[ldvl*n];
	int ldvr = n;
	double *vr = new double[ldvr*n];
	//---------------

	
   //int lwork = 3*n-1;       // size of work area
	//int lwork = 1 + 6*n + 2*n*n;
	int lwork = 4*n;

   double *work = new double[lwork];   // work area

	int liwork = 3 + 5*n;
	int *iwork = new int[liwork];	// work area

   int info = 0;

   for (j=0; j<n; j++)
   {
	  for (i=0;i<n;i++)
   	  {
   	     H[j*n+i] = A(i,j);
      }
   }

   //dsyev_(jobz, uplo, n, H, n, E, work, lwork, info);
	//dsyevd_(jobz, uplo, n, H, n, E, work, lwork, iwork, liwork, info);
	dgeev_(jobvl,jobvr,n,H,n,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info);

   for (j=0;j<n;j++)
   {
      for (i=0; i<n; i++)
      {
         //Evects(i,j)= H[j*n+i];
		 Evects(i,j) = vr[j*n+i];
      }
   }

/*   cout<<"Eigen vectors:"<<endl;
   for(i=0;i<5;i++)
   {
      for(j=0;j<5;j++)
      {
         cout<<Evects(i,j)<<' ';
      }
      cout<<endl;
   }
   cout<<endl;
*/



   for (i=0; i<n; i++)
   {
      //Evals[i] = E[i];
		Evals[i] = wr[i];
   }

   delete [] H;
   delete [] E;
   delete [] work;
 	delete [] iwork;

	delete [] wr;
	delete [] wi;
	delete [] vl;
	delete [] vr;

	return info;

}

/** 
* Find neighbor's indices
* using nearest neighbor search.
**/
void nldr_LTSA::findNeighborIdx(int centerIdx,
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





#ifndef LACKS_NAMESPACE
}
#endif






