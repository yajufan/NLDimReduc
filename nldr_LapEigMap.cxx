//
// File:        nldr_LapEigMap.cxx
//

#include <cfloat>

#ifndef included_nldr_LapEigMap
#include "nldr_LapEigMap.hxx"
#endif

#ifndef included_dggev
#include <dggev.h>
#define included_dggev
#endif


#define _DXX_



struct compareDouble
{
	bool operator () (const double &lhs, const double &rhs) const
	{
		return lhs < rhs;
	}
};



/** 
* Default Constructor
**/
nldr_LapEigMap::nldr_LapEigMap()
{
	d_nDim = 1;
	d_sigma = 5.0;
	d_weightType = 1;
}

/** 
* Copy Constructor
**/
nldr_LapEigMap::nldr_LapEigMap(const nldr_LapEigMap &rhs)
{
	*this = rhs;
}

/** 
* Constructor created with an appropriate nearest neighbor
* algorithm. The nearest neighbor object is created in the calling
* program and passed as a pointer to the base class.
**/
nldr_LapEigMap::nldr_LapEigMap(nns_NNSearch<double> *nns_algo)
{
	d_nDim = 1;
	d_sigma = 5.0;
	d_weightType = 1;
	d_nns_algo = nns_algo;
}


/** 
* Assignment operator
**/
nldr_LapEigMap &nldr_LapEigMap::operator=(const nldr_LapEigMap &rhs)
{
	*this = rhs;
}

/** 
* Destructor
**/
nldr_LapEigMap::~nldr_LapEigMap()
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


void nldr_LapEigMap::setNearestNbrAlgo(nns_NNSearch<double> *nns_algo)
{
	d_nns_algo = nns_algo;
}


void nldr_LapEigMap::setNumDimentions(int nDim)
{
	d_nDim = nDim;
}


int nldr_LapEigMap::getNumDimentions()
{
	return d_nDim;
}

void nldr_LapEigMap::setHeatKernelParam(double sigma)
{
	d_sigma = sigma;
}

double nldr_LapEigMap::getHeatKernelParam()
{
	return d_sigma;
}
/*
void nldr_LapEigMap::setWeightType(int weightType)
{
	d_weightType = weightType;
}

int nldr_LapEigMap::getWeightType()
{
	return d_weightType;
}
*/

/**
*  Applies the non-linear dimension reduction model to Input, filling
*  Output with result. 
*
*   @param Input Reference to Input instance array 
*   @param Output  Reference to Output instance array
**/
void nldr_LapEigMap::apply(const di_InstanceArray &Input,
                         di_InstanceArray &Output)
{

}

void nldr_LapEigMap::apply(const di_InstanceArray &Input_instArray,
                           di_InstanceArray &Output,
						   tbox_DenseMatrix<double> &EigVec,
						   vector<double> &eigVal,
						   vector<double> &objectiveValues,
						   int &estDim)
{
	tbox_DenseMatrix<double> Input;

	InstanceArray2DenseMatrix(Input_instArray,Input);

	unsigned int nCol = Input.getNcols();
	unsigned int nRow = Input.getNrows();	

	if(d_nDim >= nCol)
	{
		d_nDim = nCol;
	}

	if(d_nns_algo->getNumNbrs() > nRow)
	{
		cout<<"Number of neighbors required is larger than number of samples available."<<endl;
	}				

	//-------------------
	// Adjacency matrix
	//-------------------
	tbox_DenseMatrix<double> AdjacencyMx;

	cout<<"AdjacencyMx:"<<endl;
	AdjacencyMx.resize(nRow,nRow);
	adjacency(AdjacencyMx);
	
/*	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			cout<<AdjacencyMx(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;*/
	//------------------------------------------------------
	// Make sure matrix is symmetric
	// Need this especially when using k-nearest-neighbors
	//------------------------------------------------------
	symmetrize(AdjacencyMx);

	//-----------------------
	// Compute weights (W)
	//-----------------------
	cout<<"Weight Matrix:"<<endl;
	tbox_DenseMatrix<double> WeightMx;
	weights(Input,AdjacencyMx,WeightMx);

/*	cout<<"sigma = "<<d_sigma<<endl;
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			cout<<WeightMx(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
	//-------------------------
	// Diagonal weight matrix (D_{ii} = \sum_j W_{ji})
	//-------------------------
	cout<<"Diagnal Weight Matrix:"<<endl;
	tbox_DenseMatrix<double> DiagWMx;
	diagWeights(WeightMx,DiagWMx);


/*	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			cout<<DiagWMx(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
	//-------------------
	// Laplacian matrix (L = D - W)
	//-------------------
	cout<<"Laplacian Matrix:"<<endl;
	tbox_DenseMatrix<double> LapMx;
	laplacianMatrix(WeightMx,DiagWMx,LapMx);

/*	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			cout<<LapMx(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/

	//---------------------------------
	// Generalized Eigenvalue Problem
	//---------------------------------

	tbox_DenseMatrix<double> EigVec_temp;
	vector<double> eigVal_temp;
	genEigProb(LapMx,DiagWMx,EigVec_temp,eigVal_temp);

	cout<<"Ordering lambda:"<<endl;
	orderEigVals(EigVec_temp,eigVal_temp,EigVec,eigVal);


/*	cout<<"Eigen Vectors:"<<endl;
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			cout<<EigVec(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
/*	cout<<"Eigen values:"<<endl;
	for(j=0;j<nRow;j++)
	{
		cout<<eigVal[j]<<' ';
	}
	cout<<endl;
*/	

/*	cout<<"Intrinsic Dimensionality..."<<endl;
	intrinsicDim(nCol,WeightMx,EigVec,objectiveValues,estDim);

*/
	tbox_DenseMatrix<double> NewCoord;
	//------------------------------
	// New coordinates
	//------------------------------
	if(d_nDim > nCol)
	{
		// Give a warning.
		cout<<"The number of dimension required exceeds totol number of dimensions available."<<endl;
		d_nDim = nCol;
	}

	coordinates(d_nDim,eigVal,EigVec,NewCoord);
	NewCoord.resize(nRow,d_nDim);

/*	for(i=0;i<nRow;i++)
	{
		for(j=0;j<d_nDim;j++)
		{
			NewCoord.setItem(i,j,EigVec(i,j+1));
		}	
	}
*/
	wrapResults(Input_instArray,NewCoord,Output);
}

/** 
* Compute Adjacency Matrix
**/
void nldr_LapEigMap::adjacency(tbox_DenseMatrix<double> &AdjacencyMx)
{
	int i,j,iter,nRow;
	int rSize;

	vector<nns_SearchResults<double> > search_results;

	nRow = AdjacencyMx.getNrows();
	
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			AdjacencyMx(i,j) = 0.0;
		}
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
					AdjacencyMx(iter,search_results[i].index) = 1.0;

/*					if(iter != search_results[i].index)
					{
						AdjacencyMx(iter,search_results[i].index) = 1.0;
					}
*/
				}
			}
			else
			{
				cout<<"No neighbors exists for component "<<iter<<"."<<endl;
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
void nldr_LapEigMap::symmetrize(tbox_DenseMatrix<double> &AdjacencyMx)
{
	int i,j,nRow;
	
	nRow = AdjacencyMx.getNrows();

	for(i=0;i<nRow;i++){
		for(j=i+1;j<nRow;j++){
			AdjacencyMx(i,j)= max(AdjacencyMx(i,j),AdjacencyMx(j,i));
			AdjacencyMx(j,i) = AdjacencyMx(i,j);
		}
	}
}


/**
* Compute the matrix of optimal reconstruction weights 
**/
void nldr_LapEigMap::weights(tbox_DenseMatrix<double> &Input,
							 tbox_DenseMatrix<double> &AdjacencyMx,
							 tbox_DenseMatrix<double> &WeightMx)
{
	int i,j;
	int nRow = Input.getNrows();	
	vector<unsigned int> neighborIdx;
	vector<double> omega;

	WeightMx.resize(nRow,nRow);
	// Initialize
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			WeightMx(i,j) = 0.0;
		}
	}

	// Weights on neighbors
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			if( abs(AdjacencyMx(i,j)-1.0) <0.001)
			{
				WeightMx(i,j) = heatKernel(i,j,Input);
			}
		}
	}


}


/**
* Compute a matrix, each column contains nearest neighbor indices 
* ordered according to their closeness to that specific sample.
**/
void nldr_LapEigMap::neighborIdxMx(tbox_DenseMatrix<double> &NeighborhoodTable,
									   tbox_DenseMatrix<double> &AdjacencyMx)
{
	int i,iter;
	int nRow,tableWidth,rSize;
	vector<unsigned int> nnIdx;


	// Set the number of neighbors the same as the number of nearest neighbors 
	int d_nnSize;
	d_nnSize = d_nns_algo->getNumNbrs();

	nRow = NeighborhoodTable.getNrows();
	tableWidth = NeighborhoodTable.getNcols();

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

	// For using index to search 
	// (NNSearch has been modified. Now does not include search point itself.) 
	//d_nns_algo->setNumNbrs(d_nnSize+1);


	for(iter=0;iter<nRow;iter++)
	{
		if (d_nns_algo->search(iter,nnIdx))
		{			
			rSize = nnIdx.size();
		
			// First item is the sample itself (Not anymore)
			//if (rSize > d_nnSize + 1)
			if (rSize > d_nnSize)
			{
				cout<<"Number of neighbors is larger than the table width (k). Set to be k."<<endl;
				rSize = d_nnSize + 1;
			}		
			//if (rSize > 1)
			if (rSize > 0)
			{
				for(i=0;i<d_nnSize;i++)
				{
					//NeighborhoodTable(iter,i) = (int) nnIdx[i+1];
					// Column based Adjacencymx
					//AdjacencyMx(iter,nnIdx[i+1]) = 1;

					NeighborhoodTable(iter,i) = (int) nnIdx[i];
					// Column based Adjacencymx
					AdjacencyMx(iter,nnIdx[i]) = 1;
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

/*
void nldr_LapEigMap::sharedNNstrength(tbox_DenseMatrix<double> &NeighborhoodTable,
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
*/

double nldr_LapEigMap::heatKernel(int s, 
								  int t,
								  tbox_DenseMatrix<double> &Input)
{
	int j,nCol;
	double temp;
	
	nCol = Input.getNcols();


	temp = 0.0;
	for(j=0;j<nCol;j++)
	{
			temp += pow(Input(s,j)-Input(t,j),2);
	}
	return exp(-temp/d_sigma);

}

/**
* Compute the diagonal weight matrix
**/
void nldr_LapEigMap::diagWeights(tbox_DenseMatrix<double> &WeightMx,
							 	 tbox_DenseMatrix<double> &DiagWMx)
{
	int i,j,nRow;
	double temp;

	nRow = WeightMx.getNrows();
	DiagWMx.resize(nRow,nRow);

	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
				DiagWMx(i,j) = 0.0;
		}
	}

	for(i=0;i<nRow;i++)
	{
		temp = 0.0;
		for(j=0;j<nRow;j++)
		{
			temp += WeightMx(j,i);
		}
		DiagWMx(i,i) = temp;
	}
}

/**
* Compute the diagonal weight matrix
**/
void nldr_LapEigMap::laplacianMatrix(tbox_DenseMatrix<double> &WeightMx,
									 tbox_DenseMatrix<double> &DiagWMx,
								 	 tbox_DenseMatrix<double> &LapMx)
{	
	int i,j,nRow;
	
	nRow = WeightMx.getNrows();
	LapMx.resize(nRow,nRow);

	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			LapMx(i,j) = DiagWMx(i,j) - WeightMx(i,j);
			
		}
	}
}

/**
* Solve the generalized eigenvalue problem
* Obtaining the right eigenvector v(j):
* A*v(j) = lambada(j)*B*v(j)
* For Laplacian Eigenmaps,
* A = LapMx
* B = DiagWMx.
* eigVal_i = -1 if unable to calculate the i-th eigenvalue.
**/
int nldr_LapEigMap::genEigProb(tbox_DenseMatrix<double> &LapMx,
							   tbox_DenseMatrix<double> &DiagWMx,
							   tbox_DenseMatrix<double> &EigVec,
						 	   vector<double> &eigVal)
{
	int i,j;

	//---------------------------------------------
	// Set parameters for LAPACK function DGGEV()
	//---------------------------------------------
	// The order of the matrices
	int n;

	n = LapMx.getNrows();

	// 'N': Do not compute the left generalized eigenvectors
	// 'V': Compute the left generalized eigenvectors
	char jobvl = 'N'; 

	// 'N': Do not compute the right generalized eigenvectors
	// 'V': Compute the right generalized eigenvectors	
	char jobvr = 'V';

	double *A = new double[n*n];
	int lda = n;
	double *B = new double[n*n];
	int ldb = n;

	double *alphaR = new double[n];
	double *alphaI = new double[n];
	double *beta = new double[n];

	int ldvl;

	if(jobvl == 'V')
		ldvl = n;
	else
		ldvl = 1;
	// Left eigenvectors
	double *vl = new double[ldvl*n];

	int ldvr;

	if(jobvr == 'V')
		ldvr = n;
	else
		ldvr = 1;
	// Right eigenvectors
	double *vr = new double[ldvr*n];


	int lwork = 8*n;
	// Dimension of work should be max(1,lwork).
	double *work = new double[lwork];
	
	int info;

	for(j=0;j<n;j++)	
	{
		for(i=0;i<n;i++)
		{	
			A[j*n+i] = LapMx(i,j);
			B[j*n+i] = DiagWMx(i,j);
		}
	}


	dggev_(jobvl,jobvr,n,A,lda,B,ldb,alphaR,alphaI,beta,vl,ldvl,vr,ldvr,work,lwork,info);


	EigVec.resize(n,n);

	for(j=0;j<n;j++)
	{
		for(i=0;i<n;i++)
		{
			EigVec(i,j) = vr[j*n+i];
			//cout<<vr[j*n+i]<<' ';
		}
		//cout<<endl;
	}
	//cout<<endl;


	eigVal.clear();

	for(i=0;i<n;i++)
	{
		if(abs(beta[i]) < 1e-10)
		{
			cout<<"beta "<<i<<"is close to zero."<<endl;
			eigVal.push_back(-1.0);
		}
		else
		{
			eigVal.push_back(alphaR[i]/beta[i]);
		}
	}

/*
	for(i=0;i<n;i++)
	{
		cout<<beta[i]<<' ';
	}
	cout<<endl;

	cout<<info<<endl;
*/


	delete [] A;
	delete [] B;
	delete [] alphaR;
	delete [] alphaI;
	delete [] beta;
	delete [] vl;
	delete [] vr;
	delete [] work;

	if(info>0)
	{
		cout<<"Failed to compute generalized eigenvectors."<<endl;
		return 0;
	}


	return 1;
}


/* !!!!!!!!!!!Should delete this function!!!!!!!!!!! */
/**
* Order eigen values and their corresponding eigen vectors.
**/
/*void nldr_LapEigMap::orderEigVals(tbox_DenseMatrix<double> &EigVec_temp,
				  				  vector<double> &eigVal_temp,
				  				  tbox_DenseMatrix<double> &EigVec,
				  				  vector<double> &eigVal)
{
	int i,nRow,count;
	map<double,int> lambda;
	map<double,int>::iterator iter;
	

	nRow = eigVal_temp.size();

	for(i=0;i<nRow;i++)
	{		
		lambda[eigVal_temp[i]] = i;
	}
	cout<<endl;

	EigVec.resize(nRow,nRow);
	eigVal.clear();
	count = 0;
	for(iter=lambda.begin();iter != lambda.end();iter++)
	{
		eigVal.push_back((*iter).first);
		
		
		for(i=0;i<nRow;i++)
		{
			EigVec(i,count) = EigVec_temp(i,(*iter).second);
		}
		count++;
	}
}
*/

/**
* Compute new coordinates 
**/
void nldr_LapEigMap::coordinates(int nDim,
								 vector<double> &eigenValues,
							     tbox_DenseMatrix<double> &EigVec, 
							     tbox_DenseMatrix<double> &NewCoord)
{
	int i,j,nRow,nCol;

	nRow = EigVec.getNrows();	

	NewCoord.resize(nRow,nDim);
	
	nCol = eigenValues.size();


	/* Ignore the first eigenvalue entry. */
	for(j=0;j<nDim;j++)
	{
			for(i=0;i<nRow;i++)
			{
				NewCoord(i,j) = EigVec(i,j+1);
			}
	}

    /* Take out all zero eigenvalues */
/*
	count = 0;

	for(j=0;j<nCol;j++)
	{

		if(count < nDim & abs(eigenValues[j]) > 0.0001)
		{
			for(i=0;i<nRow;i++)
			{
				NewCoord(i,count) = EigVec(i,j);
			}
			count += 1;
		}
	}
*/

}

/**
* Compute residual variances for a list of possible dimensions
**/
/*void nldr_LapEigMap::intrinsicDim(tbox_DenseMatrix<double> &Input,
								tbox_DenseMatrix<double> &DiagWMx,
								tbox_DenseMatrix<double> &AdjacencyMx,
								tbox_DenseMatrix<double> &EigVec, 
						 	    vector<double> &residualVariances)
{
	int d,nFeatures;
	vector<double> DiagWMx_vector, NewDiagWMx_vector; // Distance matrix at projected space
	tbox_DenseMatrix<double> NewCoord,NewWeightMx;
	tbox_DenseMatrix<double> NewDiagWMx;

	nFeatures = Input.getNcols();
	residualVariances.clear();

	for(d=1;d<=nFeatures;d++)
	{
		coordinates(d,EigVec,NewCoord);
		
		// Compute new weight matrix
		weights(NewCoord,AdjacencyMx,NewWeightMx);
		diagWeights(NewWeightMx,NewDiagWMx);


		// Compare new weight matrix with the old one
		tboxDM2vector(NewDiagWMx,NewDiagWMx_vector);
		tboxDM2vector(DiagWMx,DiagWMx_vector);
		residualVariances.push_back(residualVar(DiagWMx_vector,NewDiagWMx_vector));
	}
}
*/
/*
void nldr_LapEigMap::intrinsicDim(int &nFeatures,
								  tbox_DenseMatrix<double> &WeightMx,
								  tbox_DenseMatrix<double> &EigVec, 
						 	      vector<double> &objectiveValues,
								  int &estDim)
{
	int d;
	tbox_DenseMatrix<double> NewCoord;

	objectiveValues.clear();

	for(d=1;d<=nFeatures;d++)
	{
		coordinates(d,EigVec,NewCoord);
		
		objectiveValues.push_back(objective(WeightMx,NewCoord));
	}
	estDim = minIndex(objectiveValues) + 1;
	
}
*/
double nldr_LapEigMap::objective(tbox_DenseMatrix<double> &WeightMx,
								 tbox_DenseMatrix<double> &NewCoord)
{
	unsigned int i,j,t,nRow,nCol;
	double obj,tempObj,temp,tempSum;
	
	
	nRow = WeightMx.getNrows();
	nCol = NewCoord.getNcols();

	if(nRow != NewCoord.getNrows())
	{
		cout<<"The Number of samples in weight matrix is not the same as the one in Y!"<<endl;
		return 0.0;
	}
	// Second version
/*	
	obj = 0.0;
	for(i=0;i<nRow;i++)
	{
		tempSum = 0.0;
		tempObj = 0.0;
		for(j=0;j<nRow;j++)
		{
			temp = 0.0;
			for(t=0;t<nCol;t++)
			{
				temp += pow((NewCoord(i,t)-NewCoord(j,t)),2);
			}
			tempSum += temp;
			tempObj += temp*WeightMx(i,j);
		}
		obj += (tempObj / tempSum);
	}
*/
	// Third version
	tempSum = 0.0;
	tempObj = 0.0;
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			temp = 0.0;
			for(t=0;t<nCol;t++)
			{
				temp += pow((NewCoord(i,t)-NewCoord(j,t)),2);
			}
			tempObj += temp*WeightMx(i,j);
			tempSum += temp;
		}
	}
	obj = tempObj/tempSum;

	return obj;
}

int nldr_LapEigMap::minIndex(vector<double> &obj)
{
	int i,n,minIdx;
	double min;

	n = obj.size();
	minIdx = 0;

	min = DBL_MAX;
	for(i=0;i<n;i++)
	{
		if(obj[i] < min)
		{
			minIdx = i;
			min = obj[i];
		}
	}
	
	return minIdx;
}

#ifndef LACKS_NAMESPACE
}
#endif






