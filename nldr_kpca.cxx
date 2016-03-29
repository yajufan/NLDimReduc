//
// File:        nldr_kpca.cxx
//


#ifndef included_nldr_kpca
#include "nldr_kpca.hxx"
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
nldr_kpca::nldr_kpca()
{
	d_nDim = 2;
	d_optDim_flag = 0;
	d_polynomialD = 2;
	d_sigma = 10;
}

/** 
* Copy Constructor
**/
nldr_kpca::nldr_kpca(const nldr_kpca &rhs)
{
	*this = rhs;
}

/** 
* Assignment operator
**/
nldr_kpca &nldr_kpca::operator=(const nldr_kpca &rhs)
{
	*this = rhs;
}

/** 
* Destructor
**/
nldr_kpca::~nldr_kpca()
{

}

void nldr_kpca::setNumDimentions(int nDim)
{
	d_nDim = nDim;
}


int nldr_kpca::getNumDimentions()
{
	return d_nDim;
}

void nldr_kpca::setOptDimFlag(int optDim_flag)
{
	d_optDim_flag = optDim_flag;
}

int nldr_kpca::getOptDimFlag()
{
	return d_optDim_flag;
}

void nldr_kpca::setPolynomialD(int polynomialD)
{
	d_polynomialD = polynomialD;
}

int nldr_kpca::getPolynomialD()
{
	return d_polynomialD;
}

void nldr_kpca::setGaussianSigma(double sigma)
{
	d_sigma = sigma;
}

double nldr_kpca::getGaussianSigma()
{
	return d_sigma;
}


void nldr_kpca::apply(const di_InstanceArray &Input,
                     di_InstanceArray &Output)
{

}

void nldr_kpca::apply(const di_InstanceArray &Input_InstArray,
					 vector<double> &eigenValues,
                     di_InstanceArray &Output)
{
	int i,j,nRow,nCol;
	tbox_DenseMatrix<double> Input;
	tbox_DenseMatrix<double> kpcaoutput;
	tbox_DenseMatrix<double> Y;
	int kpca_nDim;

	InstanceArray2DenseMatrix(Input_InstArray,Input);

	nCol = Input.getNcols();
	// Function kpca requires tbox_DenseMatrix
	//di_InstanceArrayToMatrix(Input,kpcainput);
	if(d_nDim > nCol)
	{
		cout<<"Number of desired dimensions exceeds total available. Set to total available."<<endl;
		
	}

	kpca_nDim = min(d_nDim,nCol);

	cout<<"Computing kpca..."<<endl;
	kpca(Input,kpca_nDim,eigenValues,Y);

		// Apply kpca Transform
/*		dr_kpcaTransform<double> kpcat;
		kpcat.setNumDim(nCol-1); // One of the eigenvalue is zero
		kpcat.buildAndApply(kpcainput,X);	
*/
	cout<<"Number of features after kpca: "<<Y.getNcols()<<endl;


	wrapResults(Input_InstArray,Y,Output);


}

/**
* kpca
**/
void nldr_kpca::kpca(tbox_DenseMatrix<double> &X,
					int &initDim,
					vector<double> &eigVal,
					tbox_DenseMatrix<double> &Y)
{
	int i,j,l,nRow,nCol;

	tbox_DenseMatrix<double> KernelM,CenteredKernelM;	
	tbox_DenseMatrix<double> EigVec_temp,EigVec;
	vector<double> eigVal_temp;
	double *eigVal_dbl;
	double tempSum;

	nRow = X.getNrows();
	nCol = X.getNcols();

	cout<<"nRow = "<<nRow<<endl;
	cout<<"nCol = "<<nCol<<endl;


	//------------------------
	// Compute kernel matrix
	//------------------------
	kernelMx(X,KernelM);

	//---------------------------------
	// Compute centered kernel matrix
	//---------------------------------
	centerKernelMx(KernelM,CenteredKernelM);

	//-----------------------
	// Compute Eigen Values
	//-----------------------
	EigVec_temp.resize(nRow,nRow);
	eigVal_dbl = new double[nRow];

	EigVec.resize(nRow,nRow);
	eigVal_temp.resize(nRow);

	CenteredKernelM.eig(EigVec_temp,eigVal_dbl);

/*	for(i=0;i<nRow;i++)
	{
		cout<<eigVal_dbl[i]<<endl;
	}
*/
	array2vector(eigVal_dbl,eigVal_temp,nRow);

	orderEigVals(EigVec_temp,eigVal_temp,EigVec,eigVal);


	cout<<"Eigen values: "<<endl;
	for(i=0;i<nRow;i++)
	{
		cout<<eigVal[i]<<endl;
	}
	cout<<endl;

/*	cout<<"Eigen vectors: "<<endl;
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
	// Eigen values are in increasing order
	// Want the largest ones


	Y.resize(nRow,initDim);
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<initDim;j++)
		{
			tempSum = 0.0; 
			for(l=0;l<nRow;l++)
			{
				tempSum += EigVec(l,nRow-j-1)*CenteredKernelM(i,l);
			}
			Y(i,j) = tempSum;
		}
	}
		
	cout<<"Completed kpca..."<<endl;

	delete [] eigVal_dbl;
}

/**
* Covariance Matrix
**/
void nldr_kpca::kernelMx(tbox_DenseMatrix<double> &X,
						 tbox_DenseMatrix<double> &KernelM)
{
	int i,j,l,nRow,nCol;
	double tempSum,tempSigma2;

	nRow = X.getNrows();
	nCol = X.getNcols();

	KernelM.resize(nRow,nRow);

	//------------------
	// Gaussian Kernel
	//------------------
	cout<<"Gaussian Kernel with sigma = "<<d_sigma<<endl;
	tempSigma2 = d_sigma*d_sigma;

	for(j=0;j<nRow;j++)
	{
		for(i=j;i<nRow;i++)
		{
			// norm
			tempSum = 0.0;
			for(l=0;l<nCol;l++)
			{
				tempSum += pow(X(i,l) - X(j,l),2.0);
			}
			KernelM(i,j) = exp(-tempSum/tempSigma2);
			KernelM(j,i) = KernelM(i,j);
		}
	}
	//---------------------
	// Polynomial Kernels
	//---------------------
/*	cout<<"Polynomial Kernel with D = "<< d_polynomialD <<endl;
	for(j=0;j<nRow;j++)
	{
		for(i=j;i<nRow;i++)
		{
			tempSum = 0.0;
			for(l=0;l<nCol;l++)
			{
				tempSum += X(i,l)*X(j,l);
			}	
			KernelM(i,j) = pow(tempSum,d_polynomialD);
			KernelM(j,i) = KernelM(i,j);
		}
	}
*/
}

/**
* Center the kernel matrix so the projected data has zero mean
**/
void nldr_kpca::centerKernelMx(tbox_DenseMatrix<double> &KernelM,
							   tbox_DenseMatrix<double> &CenteredKernelM)
{
	int i,j,nRow;
	double tempSum;
	vector<double> columnSum;
	double matrixSum;
	

	nRow = KernelM.getNrows();

	// Column sum
	columnSum.clear();
	matrixSum = 0.0;
	for(j=0;j<nRow;j++)
	{
		tempSum = 0.0;
		for(i=0;i<nRow;i++)
		{
			tempSum += KernelM(i,j);
		}
		columnSum.push_back(tempSum/nRow);
		matrixSum += tempSum;
	}

	matrixSum = (double) matrixSum/(nRow*nRow);

	// Centering
	CenteredKernelM.resize(nRow,nRow);
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			CenteredKernelM(i,j) = KernelM(i,j) - columnSum[j] - columnSum[i] + matrixSum;
		}
	}
}

/**
* Subtract off the mean
**/
void nldr_kpca::zeroMean(tbox_DenseMatrix<double> &Y)
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






