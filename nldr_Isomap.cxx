//
// File:        nldr_Isomap.cxx
//


#ifndef included_nldr_Isomap
#include "nldr_Isomap.hxx"
#endif

#include <stdlib.h>
#include <cfloat>

#ifndef included_dgeev
#include <dgeev.h> 
#define included_dgeev
#endif



/** 
* Default Constructor
**/
nldr_Isomap::nldr_Isomap()
{
	d_nDim = 2;
}

/** 
* Copy Constructor
**/
nldr_Isomap::nldr_Isomap(const nldr_Isomap &rhs)
{
	*this = rhs;
}

/** 
* Constructor created with an appropriate nearest neighbor
* algorithm. The nearest neighbor object is created in the calling
* program and passed as a pointer to the base class.
**/
nldr_Isomap::nldr_Isomap(nns_NNSearch<double> *nns_algo)
{
	d_nDim = 2;
	d_nns_algo = nns_algo;
}


/** 
* Assignment operator
**/
nldr_Isomap& nldr_Isomap::operator=(const nldr_Isomap &rhs)
{
	*this = rhs;

	return *this;
} 

/** 
* Destructor
**/
nldr_Isomap::~nldr_Isomap()
{
;
}
 

/**
* Set the method to be used for finding the nearest neighbors. This
* should be created in the calling program and passed to this class as a
* pointer to the base nearest-neighbor class.
*
* If this object is not set by the user, the default is the 1-nearest
* neighbor implemented using a naive search. This is set in the constructor.
**/


void nldr_Isomap::setNearestNbrAlgo(nns_NNSearch<double> *nns_algo)
{
	d_nns_algo = nns_algo;
}


void nldr_Isomap::setNumDimentions(int nDim)
{
	d_nDim = nDim;
}


int nldr_Isomap::getNumDimentions()
{
	return d_nDim;
}


/**
*  Applies the non-linear dimension reduction model to input, filling
*  output with result. 
*
*   @param input Reference to input instance array 
*   @param output  Reference to output instance array
**/
void nldr_Isomap::apply(const di_InstanceArray &input_instArray,
                        di_InstanceArray &output)
{
	tbox_DenseMatrix<double> input;
	
	InstanceArray2DenseMatrix(input_instArray,input);

	const int nRow = input.getNrows();
	int nFeatures = input.getNcols();

	vector<double> Distances;

	vector<int> keptIdx;
	tbox_DenseMatrix<double> EigVec;

	cout<< "Number of samples: "<< nRow<<endl;
	cout<< "Number of features: "<< nFeatures<<endl;

	if(d_nDim > nFeatures)
	{
		d_nDim = nFeatures;
	}

	cout<< "Adjacency distances..."<<endl;
	adjacencyDist(Distances,nRow);

	//-----------------
	// Symmetrize
	//-----------------
	cout<< "Symmetrizing..."<<endl;
	symmetrize(Distances,nRow);

	//--------------------------
	// Geodesic Distance Matrix
	//--------------------------
	cout<< "Floyd's Geodesic..."<<endl;
	floyd(Distances,nRow);

	//-----------------------------
	// Remove outliers
	//-----------------------------
	vector<double> D_cleaned;
	int nRowC;
	cout<< "Remove outliers..."<<endl;
	removeOutlier(Distances,nRow,D_cleaned,nRowC,keptIdx);

	cout<<"nRowC = "<< nRowC<<endl;
	cout<<"size of kept indices = "<< keptIdx.size()<<endl;

/*	GeoDist.resize(nRowC,nRowC);
	for(i=0;i<nRowC;i++)
	{
		for(j=0;j<nRowC;j++)
		{
			GeoDist(i,j) = D_cleaned[i*nRowC+j];
		}
	}
*/
	//--------------------
	// Double centering
	//--------------------

	// Store values in column based matrix for LAPACK
	tbox_DenseMatrix<double> MxForEigDecom;
	
	cout<< "Double centering..."<<endl;
	doubleCentering(D_cleaned,nRowC,MxForEigDecom);

	//----------------------
	// Eigen-Decomposition
	//----------------------
	double *eigVal; // Eigenvalues

	cout<< "Eigenvalue decomposition..."<<endl;
 
	EigVec.resize(nRowC,nRowC);
	eigVal = new double[nRowC];
	MxForEigDecom.eig(EigVec,eigVal);


	//-----------------
	// New Coordinate
	//-----------------
	cout<<"New coordinates... "<<endl;
	vector<double> NewCoord;

	coordinates(d_nDim,nFeatures,eigVal,EigVec,NewCoord);

	//cout<<"output..."<<endl;
	tbox_DenseMatrix<double> tbNewCoord;

	vector2tboxDM(NewCoord,tbNewCoord,nRowC,getNumDimentions());
	wrapResults(input_instArray,tbNewCoord,output);


	//-------------------
	// Residual Variance
	//-------------------
/*	cout<<"Residual variances..."<<endl;
	vector<double> D_map; // Distance matrix at projected space
	double rv_temp;
	residualVariances.clear();

	for(d=1;d<=nFeatures;d++)
	{
		coordinates(d,nRowC,eigVal,EigVec,NewCoord);
		DistMatrix(NewCoord,nRowC,d,D_map);

		rv_temp = residualVar(D_map,D_cleaned);
		residualVariances.push_back(rv_temp);

		cout<<rv_temp<<' ';
	}
	cout<<endl;

	//array2tboxDM(LeftEigVec,EigenVectors,nRowC,nRowC);
	array2vector(eigVal,eigenValues,nRowC);
*/
	delete [] eigVal;

}

void nldr_Isomap::apply(const di_InstanceArray &input_instArray,
                        di_InstanceArray &output,
						vector<int> &keptIdx,
						vector<double> &eigenValues,						
						tbox_DenseMatrix<double> &EigVec,
						tbox_DenseMatrix<double> &GeoDist)
{
	int i,j;
	tbox_DenseMatrix<double> input;
	
	InstanceArray2DenseMatrix(input_instArray,input);

	const int nRow = input.getNrows();
	int nFeatures = input.getNcols();

	vector<double> Distances;

	cout<< "Number of samples: "<< nRow<<endl;
	cout<< "Number of features: "<< nFeatures<<endl;

	if(d_nDim > nFeatures)
	{
		d_nDim = nFeatures;
	}

	cout<< "Adjacency distances..."<<endl;
	adjacencyDist(Distances,nRow);


	//-----------------
	// Symmetrize
	//-----------------
	cout<< "Symmetrizing..."<<endl;
	symmetrize(Distances,nRow);

	//--------------------------
	// Geodesic Distance Matrix
	//--------------------------
	cout<< "Floyd's Geodesic..."<<endl;
	floyd(Distances,nRow);

	//-----------------------------
	// Remove outliers
	//-----------------------------
	vector<double> D_cleaned;
	int nRowC;
	cout<< "Remove outliers..."<<endl;
	removeOutlier(Distances,nRow,D_cleaned,nRowC,keptIdx);

	cout<<"nRowC = "<< nRowC<<endl;
	cout<<"size of kept indices = "<< keptIdx.size()<<endl;

	GeoDist.resize(nRowC,nRowC);
	for(i=0;i<nRowC;i++)
	{
		for(j=0;j<nRowC;j++)
		{
			GeoDist(i,j) = D_cleaned[i*nRowC+j];
		}
	}

	//--------------------
	// Double centering
	//--------------------

	// Store values in column based matrix for LAPACK
	tbox_DenseMatrix<double> MxForEigDecom;
	
	cout<< "Double centering..."<<endl;
	doubleCentering(D_cleaned,nRowC,MxForEigDecom);

	//----------------------
	// Eigen-Decomposition
	//----------------------
	double *eigVal; // Eigenvalues

	cout<< "Eigenvalue decomposition..."<<endl;

	EigVec.resize(nRowC,nRowC);
	eigVal = new double[nRowC];
	MxForEigDecom.eig(EigVec,eigVal);


	//-----------------
	// New Coordinate
	//-----------------
	cout<<"New coordinates... "<<endl;
	vector<double> NewCoord;

	coordinates(d_nDim,nFeatures,eigVal,EigVec,NewCoord);

	tbox_DenseMatrix<double> tbNewCoord;

	vector2tboxDM(NewCoord,tbNewCoord,nRowC,getNumDimentions());
	wrapResults(input_instArray,tbNewCoord,output);


	//-------------------
	// Residual Variance
	//-------------------
/*	cout<<"Residual variances..."<<endl;
	vector<double> D_map; // Distance matrix at projected space
	double rv_temp;
	residualVariances.clear();

	for(d=1;d<=nFeatures;d++)
	{
		coordinates(d,nRowC,eigVal,EigVec,NewCoord);
		DistMatrix(NewCoord,nRowC,d,D_map);

		rv_temp = residualVar(D_map,D_cleaned);
		residualVariances.push_back(rv_temp);

		cout<<rv_temp<<' ';
	}
	cout<<endl;

	array2vector(eigVal,eigenValues,nRowC);
*/
	delete [] eigVal;


}

void nldr_Isomap::apply(const di_InstanceArray &input_instArray,
                        di_InstanceArray &output,
						vector<int> &keptIdx,
						vector<double> &residualVariances,
						vector<double> &eigenValues,						
						tbox_DenseMatrix<double> &EigVec,
						tbox_DenseMatrix<double> &GeoDist)
{
	int i,j,d;
	tbox_DenseMatrix<double> input;
	
	InstanceArray2DenseMatrix(input_instArray,input);

	const int nRow = input.getNrows();
	int nFeatures = input.getNcols();

	vector<double> Distances;

	cout<< "Number of samples: "<< nRow<<endl;
	cout<< "Number of features: "<< nFeatures<<endl;

	if(d_nDim > nFeatures)
	{
		d_nDim = nFeatures;
	}

	cout<< "Adjacency distances..."<<endl;
	adjacencyDist(Distances,nRow);

	//-----------------
	// Symmetrize
	//-----------------
	cout<< "Symmetrizing..."<<endl;
	symmetrize(Distances,nRow);

	//--------------------------
	// Geodesic Distance Matrix
	//--------------------------
	cout<< "Floyd's Geodesic..."<<endl;
	floyd(Distances,nRow);

	//-----------------------------
	// Remove outliers
	//-----------------------------
	vector<double> D_cleaned;
	int nRowC;
	cout<< "Remove outliers..."<<endl;
	removeOutlier(Distances,nRow,D_cleaned,nRowC,keptIdx);

	cout<<"nRowC = "<< nRowC<<endl;
	cout<<"size of kept indices = "<< keptIdx.size()<<endl;

	GeoDist.resize(nRowC,nRowC);
	for(i=0;i<nRowC;i++)
	{
		for(j=0;j<nRowC;j++)
		{
			GeoDist(i,j) = D_cleaned[i*nRowC+j];
		}
	}

	//--------------------
	// Double centering
	//--------------------

	// Store values in column based matrix for LAPACK
	tbox_DenseMatrix<double> MxForEigDecom;
	
	cout<< "Double centering..."<<endl;
	doubleCentering(D_cleaned,nRowC,MxForEigDecom);

	//----------------------
	// Eigen-Decomposition
	//----------------------
	double *eigVal; // Eigenvalues

	cout<< "Eigenvalue decomposition..."<<endl;

	EigVec.resize(nRowC,nRowC);
	eigVal = new double[nRowC];
	MxForEigDecom.eig(EigVec,eigVal);


	//-----------------
	// New Coordinate
	//-----------------
	cout<<"New coordinates... "<<endl;
	vector<double> NewCoord;

	coordinates(d_nDim,nFeatures,eigVal,EigVec,NewCoord);

	tbox_DenseMatrix<double> tbNewCoord;

	vector2tboxDM(NewCoord,tbNewCoord,nRowC,getNumDimentions());
	wrapResults(input_instArray,tbNewCoord,output);


	//-------------------
	// Residual Variance
	//-------------------
	cout<<"Residual variances..."<<endl;
	vector<double> D_map; // Distance matrix at projected space
	double rv_temp;
	residualVariances.clear();

	for(d=1;d<=nFeatures;d++)
	{
		coordinates(d,nRowC,eigVal,EigVec,NewCoord);
		DistMatrix(NewCoord,nRowC,d,D_map);

		rv_temp = residualVar(D_map,D_cleaned);
		residualVariances.push_back(rv_temp);

		cout<<rv_temp<<' ';
	}
	cout<<endl;

	array2vector(eigVal,eigenValues,nRowC);

	delete [] eigVal;

}


/** 
* Compute Adjacency Distance Matrix
* Input: nRow - Number of rows (samples)
* Output: Distances - Adjacency distance matrix 
*         contains distances to neighbors
* The adjacency distance matrix is stored in a STL vector 
* for row-based matrix access.
**/
void nldr_Isomap::adjacencyDist(vector<double> &Distances,const int &nRow)
{
	int i,iter;
	int rSize;

	vector<nns_SearchResults<double> > search_results;
	
	Distances.resize(nRow*nRow);

	for(i=0;i<nRow*nRow;i++)
	{
		Distances[i] = DBL_MAX;
	}

	for(i=0;i<nRow;i++)
	{
		Distances[i*nRow + i] = 0.0;
	}

	//----------------------------------------
	// Adjacency Matrix
	//
	// Not efficient because: 
	// NNSearch calculated distances twice
	//----------------------------------------

	for(iter=0;iter<nRow;iter++)
	{
		search_results.clear();
		if (d_nns_algo->search(iter,search_results))
		{
			rSize = search_results.size();

			if (rSize > 0)
			{
				for(i=0;i<rSize;i++)
				{
					Distances[iter*nRow+search_results[i].index] = search_results[i].distance;
				}
			}
			else
			{
				cout<<"No neighbors exist for component "<<iter<<"."<<endl;
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

}

/** 
* Symmetrize Adjacency Matrix
**/
void nldr_Isomap::symmetrize(vector<double> &Distances,
							 const int &nRow)
{
	int i,j;

	for(i=0;i<nRow;i++)
	{
		for(j=i+1;j<nRow;j++)
		{
			Distances[i*nRow+j] = min(Distances[i*nRow+j],Distances[j*nRow+i]);
			Distances[j*nRow+i] = Distances[i*nRow+j];
		}
	}
}

/**
* Compute a matrix, each column contains nearest neighbor indices 
* ordered according to their closeness to that specific sample.
**/
void nldr_Isomap::neighborIdxMx(tbox_DenseMatrix<double> &NeighborhoodTable,
									   tbox_DenseMatrix<double> &AdjacencyMx)
{
	int i,iter;
	int nRow,tableWidth,rSize;
	vector<unsigned int> nnIdx;

	nRow = NeighborhoodTable.getNrows();
	tableWidth = NeighborhoodTable.getNcols();

	int d_nnSize;
	d_nnSize = d_nns_algo->getNumNbrs();
	

	if(nRow < d_nnSize)
	{
		cout<<"Number of samples is less than number of neighbors. Set k = "<<nRow<<" - 1."<<endl;
		d_nnSize = nRow - 1;
		NeighborhoodTable.resize(nRow,d_nnSize);
	}
	if(tableWidth != d_nnSize)
	{
		cout<<"Set Table Width = "<< d_nnSize<<endl;
		NeighborhoodTable.resize(nRow,d_nnSize);
	}

	// (+1) For using index to search 
	d_nns_algo->setNumNbrs(d_nnSize+1);


	for(iter=0;iter<nRow;iter++)
	{
		if (d_nns_algo->search(iter,nnIdx))
		{
			rSize = nnIdx.size();
		
			// First item is the sample itself
			if (rSize > d_nnSize + 1)
			{
				cout<<"Number of neighbors is larger than the table width (k). Set to be k."<<endl;
				rSize = d_nnSize + 1;
			}		
			if (rSize > 1)
			{
				for(i=0;i<d_nnSize;i++)
				{
					NeighborhoodTable(iter,i) = (int) nnIdx[i+1];
					// Column based Adjacencymx
					AdjacencyMx(iter,nnIdx[i+1]) = 1;
				}
			}
		}
		else
		{
			cout<<"No search of neighbors has been made."<<endl;
			cout<<"Invalid search type."<<endl;
			exit(-1);
		}
	}
} 

void nldr_Isomap::sharedNNstrength(tbox_DenseMatrix<double> &NeighborhoodTable,
										  tbox_DenseMatrix<double> &AdjacencyMx,
										  tbox_DenseMatrix<double> &SNNstrength)
{
	int i,j,s,t;
	int nRow,nnSize;
	int strength;

	nRow = SNNstrength.getNrows();	
	nnSize = NeighborhoodTable.getNcols();

	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			// If i is j's neighbor and j is i's neighbor
			if (AdjacencyMx(i,j) == 1 & AdjacencyMx(j,i) == 1) 
			{
				strength = 0;
				for(s=0;s<nnSize;s++)
				{
					for(t=0;t<nnSize;t++)
					{
						if(NeighborhoodTable(i,s) == NeighborhoodTable(j,t))
						{			
							strength += (nnSize-s+1)*(nnSize-t+1);
						}
					}
				}
				SNNstrength(i,j) = strength;
				SNNstrength(j,i) = strength;
			}
		}
	}

}


/**
* Floyd's Algorithm
**/
void nldr_Isomap::floyd(vector<double> &Distances,
						const int &nRow)
{
	double *D_new;
	D_new = new double [nRow*nRow];
	
	for(int k=0;k<nRow;k++){
		for(int i=0;i<nRow;i++){
			for(int j=0;j<nRow;j++){
				D_new[i*nRow+j] = min(Distances[i*nRow+j],Distances[i*nRow+k]+Distances[k*nRow+j]);
				Distances[i*nRow+j] = D_new[i*nRow+j];
			}
		}
	}
	delete [] D_new;
}

/**
* Remove outliers
**/
void nldr_Isomap::removeOutlier(vector<double> &inDistances,
								const int &nRow,
								vector<double> &outDistances,
								int &nRowC,
								vector<int> &keptIdx)
{	 
	int i,j,count,ctemp,largestG,largestGInx;
	int *firstInx, *CfirstInx, *yInx;	

	firstInx = new int[nRow];
	CfirstInx = new int[nRow];

	largestGInx = 0;

	// Index of first point each point connects to 
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			if(inDistances[i*nRow+j]!= DBL_MAX)
			{
					firstInx[i] = j;
					break;				
			}
		}
	}

	// Size of each connected component
	for(i=0;i<nRow;i++)
	{
		CfirstInx[i] = 0;
	}
	
	for(i=0;i<nRow;i++)
	{
		CfirstInx[firstInx[i]]++;
	}

	// Find the first index of the largest group
	// Take out zeros in CfirstInx? Faster?
	//-------------------------------------------
	largestG = 0;	//count number of element in the largest group
	ctemp = 0;
	for(i=0;i<nRow;i++)
	{
		if (CfirstInx[i] != 0)
			ctemp++;	
		if (CfirstInx[i] > largestG)
		{
			largestG = CfirstInx[i];
			largestGInx = i;
		}
	}
	
	keptIdx.clear();
	if (largestG < 1)
	{
		cout<< "No connected components."<<endl;
		exit(-1);
	}
	else
	{

		cout<<"	- Number of connected components:"<<ctemp<<'\n';
		cout<<"	- Embedding component "<<largestGInx<<" with "<<largestG<<" points.\n";

		// Find indices connected with the component "largestGInx"
		yInx = new int[largestG];
		count = 0;
		for(i=0;i<nRow;i++)
		{
			if(firstInx[i]==largestGInx)
			{
				keptIdx.push_back(i);
				yInx[count] = i;
				count++;
			}
		}

		// New distance matrix contained only the largestG group
		outDistances.clear();

		for(i=0;i<largestG;i++)
		{
			for (j=0;j<largestG;j++)
			{
				outDistances.push_back(inDistances[ yInx[i]*nRow + yInx[j] ]);
			}
		}			
	} // end of else

	delete firstInx;
	delete CfirstInx;
	delete yInx;
	
	nRowC = largestG;
}

/**
* Out put a double-centered tbox_DenseMatrix
**/
void nldr_Isomap::doubleCentering(vector<double> &D_cleaned,
							 	  int nRowC,
							 	  tbox_DenseMatrix<double> &m4eigDecom)
{
	int i,j;
	tbox_DenseMatrix<double> S,H;

	S.resize(nRowC,nRowC);
	H.resize(nRowC,nRowC);
	m4eigDecom.resize(nRowC,nRowC);
	double temp;
	double val;

	// n = input.getNrows();
	for (i=0;i<nRowC;i++){
		for (j=0;j<nRowC;j++){
			temp = D_cleaned[i*nRowC+j];
			S.setItem(i,j,temp*temp);			
			if (i==j){
				val = 1.0-1.0/(double)nRowC;
			}
			else{
				val = -1.0/(double)nRowC;
			}
			H.setItem(i,j,val);
		}
	}

	H.multMatrix('N',S,'N',m4eigDecom);
	m4eigDecom*= H;
	m4eigDecom*= -0.5;		
}


/**
* Compute new coordinates in the form of STL vector 
* 
* Note that "eigVal" contains eigenvalues in increasing order
**/
void nldr_Isomap::coordinates(int nDim,
							  const int nFeatures,
							  double *eigVal,
							  tbox_DenseMatrix<double> &EigVec, 
						 	  vector<double> &NewCoord)
{
	int i,j,nRow;
	nRow = EigVec.getNrows();

	if(nDim > nFeatures)
	{
		cout<<"Requested number of dimensions exceeds the maximum available."<<endl;
		cout<<"Set the number of dimensions to the total number of features."<<endl;
		nDim = nFeatures;
	}	

	NewCoord.resize(nRow*nDim);

	for (j=0;j<nDim;j++)
	{
		//cout<<"Eigen value at D = "<< j+1<<": "<< eigVal[nRow-j-1]<<endl;
		if(eigVal[nRow-j-1] < 0)
		{
			cerr<<"Unable to compute coordinates. "<<endl;
			cerr<<"Negative eigenvalue at dimension "<<nRow-j-1<<'.'<<endl;
		}
		for(i=0;i<nRow;i++)
		{
			// Note: NewCoord represents a row-based matrix
			// while LeftEigVec represents a column-based matrix
			NewCoord[i*nDim+j] = EigVec(i,nRow-j-1)*sqrt(eigVal[nRow-j-1]);
		}
	}

}

double nldr_Isomap::min(double value1,
						double value2)
{
	if(value1 < value2)
		return value1;
	else
		return value2;
}

void nldr_Isomap::DistMatrix(vector<double> &DataMatrix,
				  			   int nRow,
				  			   int nCol,
				  			   vector<double> &DistanceMatrix)
{
	int i,j,k;
	double temp;
	
	DistanceMatrix.resize(nRow*nRow);
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			temp = 0.0;
			for(k=0;k<nCol;k++)
			{
				temp += pow(DataMatrix[i*nCol+k]-DataMatrix[j*nCol+k],2);
			}
			DistanceMatrix[i*nRow+j] = pow(temp,0.5);
			DistanceMatrix[j*nRow+i] = DistanceMatrix[i*nRow+j];
		}
	}
}

double nldr_Isomap::residualVar(vector<double> &X,
				 			    vector<double> &Y)
{
	unsigned int i;
	double tempXY,tempXX,tempYY;
	double rVar;

	unsigned int n = X.size();
	if(n != Y.size())
	{
		cout<<"Lengths of the two vectors should be the same."<<endl;
		exit(-1);
	}
	double meanX; 
	double meanY;

	meanX = mean(X);
	meanY = mean(Y);

	// Sum xy
	tempXY = 0;
	tempXX = 0;
	tempYY = 0;
	for(i=0;i<n;i++)
	{
		tempXY += X[i]*Y[i];
		tempXX += pow(X[i],2);
		tempYY += pow(Y[i],2);
	}

	rVar = 1 - pow(tempXY-n*meanX*meanY,2)/
				((tempXX-n*meanX*meanX)*(tempYY-n*meanY*meanY));

	return rVar;
}

double nldr_Isomap::mean(vector<double> &values)
{
	double temp = 0.0;
	int n = values.size();
	
	for(int i=0;i<n;i++)
	{
		temp += values[i];
	}

	return temp/n;
}
				 	 




#ifndef LACKS_NAMESPACE
}
#endif






