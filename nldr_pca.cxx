//
// File:        nldr_pca.cxx
//


#ifndef included_nldr_pca
#include "nldr_pca.hxx"
#endif

#ifndef included_nldr_NLDimRed
#include "nldr_NLDimRed.hxx"
#define included_nldr_NLDimRed
#endif

#ifndef included_cfloat
#include <cfloat>
#define included_cfloat
#endif


#define _DXX_


/** 
* Default Constructor
**/
nldr_pca::nldr_pca()
{
	d_nDim = 2;
	d_optDim_flag = 0;
}

/** 
* Copy Constructor
**/
nldr_pca::nldr_pca(const nldr_pca &rhs)
{
	*this = rhs;
}

/** 
* Assignment operator
**/
nldr_pca &nldr_pca::operator=(const nldr_pca &rhs)
{
	*this = rhs;
}

/** 
* Destructor
**/
nldr_pca::~nldr_pca()
{

}

void nldr_pca::setNumDimentions(int nDim)
{
	d_nDim = nDim;
}


int nldr_pca::getNumDimentions()
{
	return d_nDim;
}

void nldr_pca::setOptDimFlag(int optDim_flag)
{
	d_optDim_flag = optDim_flag;
}

int nldr_pca::getOptDimFlag()
{
	return d_optDim_flag;
}

void nldr_pca::apply(const di_InstanceArray &Input_InstArray,
                     di_InstanceArray &Output)
{
	int i,j,nRow,nCol;
	tbox_DenseMatrix<double> Input;
	tbox_DenseMatrix<double> PCAoutput;
	tbox_DenseMatrix<double> Y;
	int pca_nDim;
	vector<double> eigenValues;
	tbox_DenseMatrix<double> EigVec;

	InstanceArray2DenseMatrix(Input_InstArray,Input);

	// Function PCA requires tbox_DenseMatrix
	//di_InstanceArrayToMatrix(Input,PCAinput);

	nCol = Input.getNcols();
	d_nDim = nCol;
	//pca_nDim = min(d_nDim,nCol);
	pca_nDim = nCol;

	cout<<"Computing PCA..."<<endl;
	pca(Input,pca_nDim,eigenValues,EigVec,Y);

	cout<<"Number of features after PCA: "<<Y.getNcols()<<endl;

	wrapResults(Input_InstArray,Y,Output);

}

void nldr_pca::apply(const di_InstanceArray &Input_InstArray,
					 vector<double> &eigenValues,
					 tbox_DenseMatrix<double> &EigVec,
                     di_InstanceArray &Output)
{
	int i,j,nRow,nCol;
	tbox_DenseMatrix<double> Input;
	tbox_DenseMatrix<double> PCAoutput;
	tbox_DenseMatrix<double> Y;
	int pca_nDim;

	InstanceArray2DenseMatrix(Input_InstArray,Input);

	// Function PCA requires tbox_DenseMatrix
	//di_InstanceArrayToMatrix(Input,PCAinput);

	nCol = Input.getNcols();
	d_nDim = nCol;
	//pca_nDim = min(d_nDim,nCol);
	pca_nDim = nCol;

	cout<<"Computing PCA..."<<endl;
	pca(Input,pca_nDim,eigenValues,EigVec,Y);

	cout<<"Number of features after PCA: "<<Y.getNcols()<<endl;

	wrapResults(Input_InstArray,Y,Output);

}

/**
* PCA
**/
void nldr_pca::pca(tbox_DenseMatrix<double> &X,
					int &initDim,
					vector<double> &eigVal,
					tbox_DenseMatrix<double> &EigVec,
					tbox_DenseMatrix<double> &Y)
{
	int i,j,l,nRow,nCol;

	tbox_DenseMatrix<double> CovarM;	// Covariance matrix
	tbox_DenseMatrix<double> EigVec_temp;
	vector<double> eigVal_temp;
	double *eigVal_dbl;
	double tempSum;

	nRow = X.getNrows();
	nCol = X.getNcols();

	cout<<"nRow = "<<nRow<<endl;
	cout<<"nCol = "<<nCol<<endl;

	//----------------------------------------------------------
	// Subtract off the mean for each dimension
	//----------------------------------------------------------
	zeroMean(X);

	//--------------------------------
	// Compute the covariance matrix
	//--------------------------------
	covarMx(X,CovarM);

	EigVec_temp.resize(nCol,nCol);
	eigVal_dbl = new double[nCol];

	EigVec.resize(nCol,nCol);
	eigVal_temp.resize(nCol);

	CovarM.eig(EigVec_temp,eigVal_dbl);

/*	for(i=0;i<nCol;i++)
	{
		cout<<eigVal_dbl[i]<<endl;
	}
*/
	array2vector(eigVal_dbl,eigVal_temp,nCol);

	orderEigVals(EigVec_temp,eigVal_temp,EigVec,eigVal);

/*
	cout<<"Eigen values: "<<endl;
	for(i=0;i<nCol;i++)
	{
		cout<<eigVal[i]<<endl;
	}
	cout<<endl;

	cout<<"Eigen vectors: "<<endl;
	for(i=0;i<nCol;i++)
	{
		for(j=0;j<nCol;j++)
		{
			cout<<EigVec(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
	// Eigen values are in increasing order
	// Want the largest ones


	Y.resize(nRow,initDim);
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<initDim;j++)
		{
			tempSum = 0.0; 
			for(l=0;l<nCol;l++)
			{
				tempSum += EigVec(l,nCol-j-1)*X(i,l);
			}
			Y(i,j) = tempSum;
		}
	}
		
	cout<<"Completed PCA..."<<endl;

	delete [] eigVal_dbl;
}

/**
* Covariance Matrix
**/
void nldr_pca::covarMx(tbox_DenseMatrix<double> &X,
						tbox_DenseMatrix<double> &CovarM)
{
	int i,j,l,nRow,nCol;
	double tempSum;

	nRow = X.getNrows();
	nCol = X.getNcols();

	CovarM.resize(nCol,nCol);

	for(j=0;j<nCol;j++)
	{
		for(i=j;i<nCol;i++)
		{
			tempSum = 0.0;
			for(l=0;l<nRow;l++)
			{
				tempSum += X(l,i)*X(l,j);
			}	
			CovarM(i,j) = tempSum / (nRow-1);
			//CovarM(i,j) = tempSum;
			CovarM(j,i) = CovarM(i,j);
		}
	}
}


/**
* Subtract off the mean
**/
void nldr_pca::zeroMean(tbox_DenseMatrix<double> &Y)
{
	int i,j,nRow,nCol;
	double sTemp;

	nRow = Y.getNrows();
	nCol = Y.getNcols();

	for(j=0;j<nCol;j++)
	{
		sTemp = 0.0;
		for(i=0;i<nRow;i++)
		{
			sTemp += Y(i,j);
		}
		sTemp /= nRow;
		
		for(i=0;i<nRow;i++)
		{
			Y(i,j) -= sTemp; 
		}
	}	
}


#ifndef LACKS_NAMESPACE
}
#endif






