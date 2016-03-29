//
// File:        nldr_tSNE.cxx
//


#ifndef included_nldr_tSNE
#include "nldr_tSNE.hxx"
#endif

#ifndef included_dr_PCATransform
#include "dr_PCATransform.hxx"
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
nldr_tSNE::nldr_tSNE()
{
	d_nDim = 2;
	d_perp = 15.0;
	d_multiP_flag = 1;
	d_multiP = 4.0;
	d_multiP_stop_iter = 100;
	d_tolerance = 0.00001;
	d_optTolerance = 0.00001;
	d_numIter = 1000;
	d_alpha_init = 0.5;
	d_alpha = 0.8;
	d_alpha_switch_iter = 250;
	d_eta_init = 500.0;
	d_adapt_multi = 0.8;
	d_adapt_plus = 0.2;
	d_adapt_min = 0.01;
	d_adapt_flag = 1;
	d_initSol_flag = 1;
	d_pca_flag = 1;
	d_srand_flag = 1;
	d_pca_nDim = 100;
}

/** 
* Copy Constructor
**/
nldr_tSNE::nldr_tSNE(const nldr_tSNE &rhs)
{
	*this = rhs;
}

/** 
* Constructor created with an appropriate nearest neighbor
* algorithm. The nearest neighbor object is created in the calling
* program and passed as a pointer to the base class.
**/
/*nldr_tSNE::nldr_tSNE()
{
	d_nDim = 2;
	d_perp = 15.0;
	d_numIter = 1000;
	d_alpha_init = 0.5;
	d_alpha = 0.8;
	d_alpha_switch_iter = 250;
	d_eta_init = 100.0;

}
*/

/** 
* Assignment operator
**/
nldr_tSNE &nldr_tSNE::operator=(const nldr_tSNE &rhs)
{
	*this = rhs;
}

/** 
* Destructor
**/
nldr_tSNE::~nldr_tSNE()
{

}

void nldr_tSNE::setNumDimentions(int nDim)
{
	d_nDim = nDim;
}


int nldr_tSNE::getNumDimentions()
{
	return d_nDim;
}

void nldr_tSNE::setPerplexity(double perp)
{
	d_perp = perp;
}

double nldr_tSNE::getPerplexity()
{
	return d_perp;
}

void nldr_tSNE::setMultiPflag(int multiP_flag)
{
	d_multiP_flag = multiP_flag;
}

int nldr_tSNE::getMultiPflag()
{
	return d_multiP_flag;
}

void nldr_tSNE::setMultiP(double multiP)
{
	d_multiP = multiP;
}

double nldr_tSNE::getMultiP()
{
	return d_multiP;
}

void nldr_tSNE::setMultiPstopIter(int multiP_stop_iter)
{
	d_multiP_stop_iter = multiP_stop_iter;
}

int nldr_tSNE::getMultiPstopIter()
{
	return d_multiP_stop_iter;
}

void nldr_tSNE::setTolerance(double tolerance)
{
	d_tolerance = tolerance;
}

double nldr_tSNE::getTolerance()
{
	return d_tolerance;
}

void nldr_tSNE::setOptTolerance(double optTolerance)
{
	d_optTolerance = optTolerance;
}

double nldr_tSNE::getOptTolerance()
{
	return d_optTolerance;
}

void nldr_tSNE::setNumIter(int numIter)
{
	d_numIter = numIter;
}

int nldr_tSNE::getNumIter()
{
   return d_numIter;			
} 

/** Initial Momentum **/
void nldr_tSNE::setAlphaInit(double alpha_init)
{
	d_alpha_init = alpha_init;
}

double nldr_tSNE::getAlphaInit()
{
   return d_alpha_init;     
}

/** Momentum **/
void nldr_tSNE::setAlpha(double alpha)
{
	d_alpha = alpha;
}

double nldr_tSNE::getAlpha()
{   
   return d_alpha;          
}

/** The number of iterations that switches the alpha value **/
void nldr_tSNE::setAlphaSwitchIter(int alpha_switch_iter)
{
	d_alpha_switch_iter = alpha_switch_iter;
}

int nldr_tSNE::getAlphaSwitchIter()
{
   return d_alpha_switch_iter;	
} 

/** Inital learning rate **/
void nldr_tSNE::setEtaInit(double eta_init)
{
	d_eta_init = eta_init;
}

double nldr_tSNE::getEtaInit()
{
	return d_eta_init;		
}

/** Constant Adaptive learning rate, the multiplication **/
void nldr_tSNE::setAdaptMulti(double adapt_multi)
{
	d_adapt_multi = adapt_multi;
}

double nldr_tSNE::getAdaptMulti()
{
	return d_adapt_multi;		
}

/** Constant adaptive learning rate, the plus **/
void nldr_tSNE::setAdaptPlus(double adapt_plus)
{
	d_adapt_plus = adapt_plus;
}

double nldr_tSNE::getAdaptPlus()
{
	return d_adapt_plus;		
}

/** Constant adaptive learning rate, the min **/
void nldr_tSNE::setAdaptMin(double adapt_min)
{
	d_adapt_min = adapt_min;
}

double nldr_tSNE::getAdaptMin()
{
	return d_adapt_min;		
}

/** Constant adaptive learning rate, the flag **/
void nldr_tSNE::setAdaptFlag(int adapt_flag)
{
	d_adapt_flag = adapt_flag;
}

int nldr_tSNE::getAdaptFlag()
{
	return d_adapt_flag;		
}

/** Initial solution flag **/
void nldr_tSNE::setInitSolFlag(int initSol_flag)
{
	d_initSol_flag = initSol_flag;
}

int nldr_tSNE::getInitSolFlag()
{
	return d_initSol_flag;		
}

/** PCA flag **/
void nldr_tSNE::setPCAFlag(int pca_flag)
{
	d_pca_flag = pca_flag;
}

int nldr_tSNE::getPCAFlag()
{
	return d_pca_flag;		
}

/** random seed flag **/
void nldr_tSNE::setSrandFlag(int srand_flag)
{
	 d_srand_flag = srand_flag;
}

int nldr_tSNE::getSrandFlag()
{
	return d_srand_flag;		
}

/** Number of dimensions output from PCA **/
void nldr_tSNE::setPCAnDim(int pca_nDim)
{
	d_pca_nDim = pca_nDim;
}


int nldr_tSNE::getPCAnDim()
{
	return d_pca_nDim;
}


/**
*  Applies the non-linear dimension reduction model to Input, filling
*  Output with result. 
*
*   @param Input Reference to Input instance array 
*   @param Output  Reference to Output instance array
**/
/**
*  Applies the non-linear dimension reduction model to Input, filling
*  Output with result. 
*
*   @param Input Reference to Input instance array 
*   @param Output  Reference to Output instance array
**/
void nldr_tSNE::apply(const di_InstanceArray &Input_InstArray,
                         di_InstanceArray &Output)
{
	int i,j,nRow,nCol;
	tbox_DenseMatrix<double> Input,Xnormalized,X,Dsquared;
	tbox_DenseMatrix<double> PCAoutput;
	tbox_DenseMatrix<double> Y;
	tbox_DenseMatrix<double> InitialSolution;
	tbox_DenseMatrix<double> Pcond;
	vector<double> beta;
	int pca_nDim;

	d_initSol_flag = 0;

	InstanceArray2DenseMatrix(Input_InstArray,Input);

	nRow = Input.getNrows();
	nCol = Input.getNcols();

	if(d_nDim > nCol)
	{
		d_nDim = nCol;
		cout<<"Number of desired dimensions exceeds the maximum available."<<endl;
	}
	cout<<"Dimention = " << d_nDim <<endl;

	// Normalize Input Data
	
	normalize(Input,Xnormalized);

	if(d_pca_flag == 1)
	{
		// Function PCA requires tbox_DenseMatrix
		//di_InstanceArrayToMatrix(Input,PCAinput);

		pca_nDim = min(d_nDim,nCol);

		cout<<"Computing PCA..."<<endl;
		pca(Xnormalized,pca_nDim,X);

		cout<<"Number of features after PCA: "<<X.getNcols()<<endl;

	}
	else
	{
		// Convert Input to X
		//di_InstanceArrayToMatrix(Input,X);
		X.resize(nRow,nCol);
		for(j=0;j<nCol;j++)
		{
			for(i=0;i<nRow;i++)
			{
				X(i,j) = Xnormalized(i,j);
			}
		}
	}


	// Compute the distance matrix
	squaredDistMx(X,Dsquared);
	// Compute pairwise affinities p with perplexity Prep
	d2p(Dsquared,Pcond,beta);


/*	cout<<"Initial solution:"<<InitSol.getNumFeatures()<<endl;

	if(InitSol.getNumFeatures() < d_nDim)
	{
		cout<<"Dimension of the inital solution is less the desired dimension, "<<d_nDim<<"."<<endl;
		cout<<"Use random inital solution."<<endl;
		d_initSol_flag = 0;
	}
	else if(InitSol.getNumInstances() != nRow)
	{
		cout<<"The number of instances of initial solution is not sufficient."<<endl;
		cout<<"Use random inital solution."<<endl;
		d_initSol_flag = 0;		
	}
	else
	{
		cout<<"Use the given inital solution."<<endl;
		// Instance Array to Dense Matrix
		di_InstanceArrayToMatrix(InitSol,InitialSolution);
		d_initSol_flag == 1;
	}
*/
	cout<<"t-SNE algorithm"<<endl;
	// t-SNE algorithm
	tSNEp(Pcond,Y,InitialSolution);

	cout<<"Number of featurs before wraping results: "<< Y.getNcols()<<endl;
	
	wrapResults(Input_InstArray,Y,Output);
}
void nldr_tSNE::apply(const di_InstanceArray &Input_InstArray,
                      di_InstanceArray &Output,
					  di_InstanceArray &InitSol)
{
	int i,j,nRow,nCol;
	tbox_DenseMatrix<double> Input,Xnormalized,X,Dsquared;
	tbox_DenseMatrix<double> PCAoutput;
	tbox_DenseMatrix<double> Y;
	tbox_DenseMatrix<double> InitialSolution;
	tbox_DenseMatrix<double> Pcond;
	vector<double> beta;
	double sigma;
	int pca_nDim;


	InstanceArray2DenseMatrix(Input_InstArray,Input);

	nRow = Input.getNrows();
	nCol = Input.getNcols();

	if(d_nDim > nCol)
	{
		d_nDim = nCol;
		cout<<"Number of desired dimensions exceeds the maximum available."<<endl;
	}
	cout<<"Dimention = " << d_nDim <<endl;

	// Normalize Input Data
	
	normalize(Input,Xnormalized);

/*	Xnormalized.resize(nRow,nCol);
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nCol;j++)
		{
			Xnormalized(i,j) = Input(i,j);
		}	
	}
*/

/*	for(i=0;i<nRow;i++)
	{
		for(j=0;j<Xnormalized.getNcols();j++)
		{

			cout<<Xnormalized(i,j) <<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
	
	if(d_pca_flag == 1)
	{
		// Function PCA requires tbox_DenseMatrix
		//di_InstanceArrayToMatrix(Input,PCAinput);

		pca_nDim = min(d_nDim,nCol);

		cout<<"Computing PCA..."<<endl;
		pca(Xnormalized,pca_nDim,X);

/*		cout<<"Data after PCA:"<<endl;
		for(i=0;i<nRow;i++)
		{
			for(j=0;j<X.getNcols();j++)
			{
				cout<<X(i,j)<<' ';
			}
			cout<<endl;
		}
*/
		// Apply PCA Transform
/*		dr_PCATransform<double> pcat;
		pcat.setNumDim(nCol-1); // One of the eigenvalue is zero
		pcat.buildAndApply(PCAinput,X);	
*/
		cout<<"Number of features after PCA: "<<X.getNcols()<<endl;

	}
	else
	{
		// Convert Input to X
		//di_InstanceArrayToMatrix(Input,X);
		X.resize(nRow,nCol);
		for(j=0;j<nCol;j++)
		{
			for(i=0;i<nRow;i++)
			{
				X(i,j) = Xnormalized(i,j);
			}
		}
	}

/*	for(i=0;i<nRow;i++)
	{
		for(j=0;j<X.getNcols();j++)
		{

			cout<<X(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/


	// Compute the distance matrix
	squaredDistMx(X,Dsquared);
/*
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{

			cout<<Dsquared(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
	// Compute pairwise affinities p with perplexity Prep
	d2p(Dsquared,Pcond,beta);
/*
	sigma = 0.0;
	for(j=0;j<nRow;j++)
	{
		sigma += sqrt(1/beta[j]);
	}
	cout<<"Mean value of sigma: "<<sigma/nRow<<endl;
*/

/*	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			cout<<Pcond(i,j)<<' ';
		}
		cout<<endl;
	}
*/
	//----------------------------
	// Print out Pcond for MATLAB
	//----------------------------
/*	string location("");
	location += "../../tools/tSNE/Simple_tSNE/matrix.dat";
	ofstream testD(location.data());

	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			testD<<Pcond(i,j)<<'\t';
		}
		testD<<endl;
	}
	testD.close();
	*/


	//---------------------------
	// Testing Pcond from MATLAB
	//---------------------------
/*	Pcond.resize(nRow,nRow);
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			Pcond(i,j) = Input[i][j];
		}
	}*/
	//---------------------------

	cout<<"Initial solution:"<<InitSol.getNumFeatures()<<endl;

	if(InitSol.getNumFeatures() < d_nDim)
	{
		cout<<"Dimension of the inital solution is less the desired dimension, "<<d_nDim<<"."<<endl;
		cout<<"Use random inital solution."<<endl;
		d_initSol_flag = 0;
	}
	else if(InitSol.getNumInstances() != nRow)
	{
		cout<<"The number of instances of initial solution is not sufficient."<<endl;
		cout<<"Use random inital solution."<<endl;
		d_initSol_flag = 0;		
	}
	else
	{
		cout<<"Use the given inital solution."<<endl;
		// Instance Array to Dense Matrix
		di_InstanceArrayToMatrix(InitSol,InitialSolution);
		d_initSol_flag == 1;
	}

	cout<<"t-SNE algorithm"<<endl;
	// t-SNE algorithm
	tSNEp(Pcond,Y,InitialSolution);

	cout<<"Number of featurs before wraping results: "<< Y.getNcols()<<endl;
	
	wrapResults(Input_InstArray,Y,Output);
}

/**
* Normalize Input Data
* X = (X - min(X)) / max(X)
**/
void nldr_tSNE::normalize(tbox_DenseMatrix<double> &Input,
						  tbox_DenseMatrix<double> &X)
{
	int i,j;
	int nRow,nCol;
	double minX,maxX,diff;

	nRow = Input.getNrows();
	nCol = Input.getNcols();
	
	// Find minimum and maximum of X
	minX = DBL_MAX;
	maxX = -DBL_MAX;
	for(j=0;j<nCol;j++)	
	{
		for(i=0;i<nRow;i++)
		{
			if(Input(i,j) < minX)
			{
				minX = Input(i,j);
			}
			if(Input(i,j) > maxX)
			{
				maxX = Input(i,j);
			}
		}
	}

	diff = maxX - minX;
	X.resize(nRow,nCol);
	for(j=0;j<nCol;j++)
	{
		for(i=0;i<nRow;i++)
		{
			X(i,j) = (Input(i,j) - minX)/diff;
		}
	}
}


/**
* PCA
**/
void nldr_tSNE::pca(tbox_DenseMatrix<double> &X,
					int &initDim,
					tbox_DenseMatrix<double> &Y)
{
	int i,j,l,nRow,nCol;

	tbox_DenseMatrix<double> CovarM;	// Covariance matrix
	tbox_DenseMatrix<double> EigVec_temp,EigVec;
	vector<double> eigVal, eigVal_temp;
	double *eigVal_dbl;
	double tempSum;

	nRow = X.getNrows();
	nCol = X.getNcols();

	cout<<"nRow = "<<nRow<<endl;
	cout<<"nCol = "<<nCol<<endl;


	//----------------------------------------------------------
	// Subtract off the mean for each dimension
	// ?? tSNE Matlab code has zero mean here... maybe wrong!!
	//----------------------------------------------------------
	zeroMean(X);

	// Compute the covariance matrix
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


/*	cout<<"Eigen values: "<<endl;
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
}

/**
* Covariance Matrix
**/
void nldr_tSNE::covarMx(tbox_DenseMatrix<double> &X,
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
* Compute squared distance matrix
* for Gaussian
**/
void nldr_tSNE::squaredDistMx(tbox_DenseMatrix<double> &X,
							  tbox_DenseMatrix<double> &Dsquared)
{
	int i,j,l,nRow,nCol;
	double dTemp;
	
	nRow = X.getNrows();
	nCol = X.getNcols();

	Dsquared.resize(nRow,nRow);

	for(i=0;i<nRow;i++)
	{
		for(j=i+1;j<nRow;j++)
		{
			// Compute squared Euclidean 
			dTemp = 0.0;
			for(l=0;l<nCol;l++)
			{
				dTemp += pow(X(i,l) - X(j,l),2);
			}
			Dsquared(i,j) = dTemp;
			Dsquared(j,i) = dTemp;
		}
	}
}

/**
* Compute conditional probabilities
* Identifies appropriate sigma's that satisfy the given perplexity
* up to some tolerance.
* Input: D.^2, perplexity and tolerance
* Output: Conditional probabilities, beta=1/(2*(sigma^2))
**/
void nldr_tSNE::d2p(tbox_DenseMatrix<double> &Dsquared,
					tbox_DenseMatrix<double> &Pcond,
					vector<double> &beta)
{
	int i,ii,iter,j,nRow;
	double logPerp; // log of perplexity = entropy (H) 
	double beta_s,beta_max,beta_min;
	double h,h_diff;
	int count;
	vector<double> PcondV;

	int max_iter = 50;

	Timer timer=Timer();


	nRow = Dsquared.getNrows();

	// Initialize variables
	logPerp = log2(d_perp);
	//logPerp = log(d_perp);
	Pcond.resize(nRow,nRow);

		
	beta.clear();

	// Run over all data points
	for(i=0;i<nRow;i++)
	{	
		beta_s = 1.0; 
		beta_max = DBL_MAX;
		beta_min = -DBL_MAX;

		// Conditional probabilities
		condProb(Dsquared,i,beta_s,PcondV);
/*
		double sumPtemp;
		sumPtemp = 0.0;
		for (ii=0;ii<nRow;ii++)
		{
			sumPtemp += Pcond(ii,i);
		}
		cout<<"sum of Pcond = "<< sumPtemp<<endl;
*/

		// Shannon entropy
		shannonEntropy(PcondV,i,h);

		//cout<<"h = "<<h<<endl;
		//cout<<"logPerp = "<<logPerp<<endl;

		// Evaluate whether the perplexity is within tolerance
		h_diff = h - logPerp;
		count = 0;

		while(abs(h_diff) > d_tolerance && count < max_iter)
		{
			// If h is higher than desired, increase beta.
			if(h_diff > 0.0)
			{
				//cout<<"-> ";
				beta_min = beta_s;
				if(beta_max == DBL_MAX)
				{
					beta_s = beta_s*2.0;
				}
				else
				{
					beta_s = (beta_s + beta_max)/2.0;
				}
			}
			else // If h is lower than desired, decrease beta
			{
				//cout<<"<- ";
				beta_max = beta_s;
				if(beta_min == -DBL_MAX)
				{
					beta_s = beta_s/2.0;
				}
				else
				{
					beta_s = (beta_s + beta_min)/2.0;
				}
			}
			// Recompute conditional probilities and Shannon entropy
			// Conditional probabilities
			condProb(Dsquared,i,beta_s,PcondV);
			// Shannon entropy
			shannonEntropy(PcondV,i,h);

			h_diff = h - logPerp;
			count += 1;
		}

		beta.push_back(beta_s);

		if(i%500 < 0.1)
		{		
			cout<<"Computed P-values "<<i<< " of "<<nRow<<" components..."<<endl;
		}

		for(j=0;j<nRow;j++)
		{
			Pcond(j,i) = PcondV[j];
		}
/*
		if(count < max_iter)
		{
			cout<<"h_diff"<< h_diff<<endl;
		}		
*/
	}
}


/**
* For a fixed point, idx,
* Compute temporary conditional probabilities given beta=1/(2*(sigma^2))
**/
// The for loops adapt the column based matrix.
void nldr_tSNE::condProb(tbox_DenseMatrix<double> &Dsquared,
						 int &idx,
						 double &beta_s,
						 vector<double> &PcondV)
{
	int i,j,nRow;
	double numerator, denominator;
	vector<double> DsTimesBeta;

	nRow = Dsquared.getNrows();
	DsTimesBeta.clear();
	PcondV.clear();

	for(j=0;j<nRow;j++)
	{
		DsTimesBeta.push_back(exp(-Dsquared(j,idx)*beta_s));
	}

	denominator = 0.0;
	for(j=0;j<nRow;j++)
	{
		if(j != idx)
		{
			denominator += DsTimesBeta[j];
		}
	}	

	for(j=0;j<nRow;j++)
	{	
		PcondV.push_back( DsTimesBeta[j]/denominator );
	}
	PcondV[idx] = 0.0;
}


/**
* For a fixed point, idx,
* Compute Shannon Entropy, H(P_i), using joint probabilities
**/
// The for loops adapt the column based matrix.
void nldr_tSNE::shannonEntropy(vector<double> &PcondV,
							   int &idx,
							   double &h)
{
	int i,nRow;
	double h_temp;
	
	nRow = PcondV.size();
	
	h_temp = 0.0;
	for(i=0;i<nRow;i++)
	{
		if(i != idx)
		{
			h_temp += PcondV[i]*log2(PcondV[i]);
			//h_temp += PcondV[i]*log(PcondV[i]);
		}
	}
	h = -h_temp;
}

/**
* Given conditional probabilities, apply tSNE and 
* output new coordinates
* If there is no initial solution, set initial solution flag = 0.
**/
void nldr_tSNE::tSNEp(tbox_DenseMatrix<double> &Pcond,
					  tbox_DenseMatrix<double> &Y,
					  tbox_DenseMatrix<double> &InitialSolution)
{
	int i,j,iter,iCount,nRow;

	nRow = Pcond.getNrows();

	tbox_DenseMatrix<double> Ptemp,Ptemp2,P,Ptrue,Plied;
	// tbox_DenseMatrix<double> Y;
	tbox_DenseMatrix<double> numeratorM,Qtemp,Q;
	tbox_DenseMatrix<double> Y_grads,Y_incs;
	tbox_DenseMatrix<double> Adapt;
	double eta,alpha;
	double constCost,obj;
	double grads_norm;

	//Timer timer=Timer();

	cout<<"Symmetrize P"<<endl;
	// Symmetrize p
	symmetrizeP(Pcond,Ptemp);
					
	// Make sure P sum to one and 
	// prevent having P too close to zero
	PsumOne(Ptemp,Ptrue);

	// Prevent having P too close to zero
	//preventZero(Ptemp2,Ptrue);

	// Compute constant part of objective function before exaggerating P
	constCostFun(Ptrue,constCost);
	cout<<"Constant term in cost funtion = "<<constCost<<endl;

	// Exaggerate P
	P.resize(nRow,nRow);
	if(d_multiP_flag == 1)
	{
		exaggerateP(Ptrue,Plied);
		
		for(j=0;j<nRow;j++)
		{
			for(i=0;i<nRow;i++)
			{	
				P(i,j) = Plied(i,j);
			}
		}
	}
	else
	{
		for(j=0;j<nRow;j++)
		{
			for(i=0;i<nRow;i++)
			{	
				P(i,j) = Ptrue(i,j);
			}
		}
	}

	// Compute initial solution
	Y.resize(nRow,d_nDim);
	if(d_initSol_flag == 1)
	{
		for(j=0;j<d_nDim;j++)
		{
			for(i=0;i<nRow;i++)
			{
				Y(i,j)= InitialSolution(i,j);
			}
		}
	}
	else
	{
		initialY(Y);
	}


//-------------------------------
// MATLAB values
/*	double Ytemp[18] = { 
				 -0.0001616,-0.0000492,-0.0000969,-0.0000474,
   				 -0.0000666,0.0000478,-0.0001041,-0.0000821,
				  0.0000101,0.0001084,-0.0001855, 0.0001800,
				 -0.0000878,0.0000648, 0.0001050, 0.0001323,-0.0000517,0.0000024};
	
	for(j=0;j<d_nDim;j++)
	{
		for(i=0;i<nRow;i++)
		{
			Y(i,j)= Ytemp[j*nRow+i];
			cout<<Y(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
//----------------------------------
*/

	cout<<"First element in Y: "<<Y(0,0)<<endl;

	Y_incs.resize(nRow,d_nDim);
	for(j=0;j<d_nDim;j++)
	{
		for(i=0;i<nRow;i++)
		{
			Y_incs(i,j) = 0.0;
		}
	}

	alpha = d_alpha_init;
	eta = d_eta_init;	

	Adapt.resize(nRow,d_nDim);
	for(j=0;j<d_nDim;j++)
	{
		for(i=0;i<nRow;i++)
		{
			Adapt(i,j) = 1.0;
		}
	}


	// Optimization iterations
	iCount = 0;
	//for(iter=0;iter<d_numIter;iter++)
	//while(iCount < 3)
	while(iCount < d_numIter)
	{		
		
		iCount += 1;
		iter = iCount;

		// Compute Q
//		cout<<"Compute Q"<<endl;
		mappedQ(Y,numeratorM,Qtemp);
		PsumOne(Qtemp,Q);

/*		if(iCount == 1)
		{
		cout<<"Q:"<<endl;
		for(i=0;i<nRow;i++)
		{
			for(j=0;j<nRow;j++)
			{
				cout<<Q(i,j)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;

		}
*/
		// Compute the gradients
//		cout<<"Gradients"<<endl;		
		gradients(P,Q,numeratorM,Y,Y_grads);

/*		if(iCount == 1)
		{
		cout<<"Gradients:"<<endl;
		for(i=0;i<nRow;i++)
		{
			for(j=0;j<d_nDim;j++)
			{
				cout<<Y_grads(i,j)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;

		}
*/
		// Compute norm of the gradient
//		cout<<"Norm"<<endl;
		grads_norm = Gnorm(Y_grads);

		// Update the solution
//		cout<<"Update solution"<<endl;
		updateSolution(eta,alpha,Adapt,Y_grads,Y_incs,Y);

/*		cout<<"Y incs:"<<endl;
		for(i=0;i<nRow;i++)
		{
			for(j=0;j<d_nDim;j++)
			{
				cout<<Y_incs(i,j)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;
*/
/*		if(iCount == 1)
		{
		cout<<"Updated Y:"<<endl;
		for(i=0;i<nRow;i++)
		{
			for(j=0;j<d_nDim;j++)
			{
				cout<<setprecision(15)<<Y(i,j)<<' ';
			}
			cout<<endl;
		}
		cout<<endl;
		}
*/
		// Update the momentum if necessary

		if(iter == d_alpha_switch_iter)
		{
			alpha = d_alpha;
		}
	
		// Stop exaggerating
		if(d_multiP_flag == 1)
		{
			if(iter == d_multiP_stop_iter)
			{
				for(j=0;j<nRow;j++)
				{
					for(i=0;i<nRow;i++)
					{
						P(i,j) = Ptrue(i,j);
					}
				}
				d_multiP_flag = 0;

/*				for(i=0;i<nRow;i++)
				{
					for(j=0;j<nRow;j++)
					{
						cout<<P(i,j)<<' ';
					}
					cout<<endl;
				}
				cout<<endl;
*/
			}
		}	

		if(iter%10 < 0.1)
		{
			// Print out progress
			costFunction(P,Q,constCost,obj);
		
			cout<<"At iter = "<<iter<< ",\tcost function = "<<setprecision(15)<<obj<<endl;
		}
		
		if(grads_norm < d_optTolerance & d_multiP_flag == 0)
		{
			// cout<<"At iter = "<<iter<<"Norm of the Gradient is close to zero."<<endl;
			
			// break;
		}

	}

	cout<<"Total number of iterations: "<<iCount<<endl;
	// Print out optimal solution
/*	cout<<"Y:"<<endl;
	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nRow;j++)
		{
			cout<<Y(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;*/
	costFunction(P,Q,constCost,obj);
	cout<< "Cost function = "<<obj<<endl;
}

/**
* Symmetrize probabilities
**/
void nldr_tSNE::symmetrizeP(tbox_DenseMatrix<double> &Pcond,
							tbox_DenseMatrix<double> &P)
{
	int i,j,nRow;

	nRow = Pcond.getNrows();

	P.resize(nRow,nRow);

	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			//P(i,j) = 0.5*(Pcond(i,j)+Pcond(j,i))/nRow;
			P(i,j) = 0.5*(Pcond(i,j)+Pcond(j,i));
		}
	}
}

/**
* Make sure P sum to one
**/
void nldr_tSNE::PsumOne(tbox_DenseMatrix<double> &Ptemp,
						tbox_DenseMatrix<double> &P)
{
	int i,j,nRow;
	double sum;	

	nRow = Ptemp.getNrows();
	
	P.resize(nRow,nRow);

	sum = 0.0;
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			sum += Ptemp(i,j);
		}
	}
		
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			if(i!=j)
			{
				P(i,j) = max(Ptemp(i,j)/sum,DBL_MIN);
				
			}
		}
		P(j,j) = 0.0;
/*		for(i=j+1;i<nRow;i++)
		{
			P(i,j) = max(Ptemp(i,j)/sum,DBL_MIN);
			P(j,i) = P(i,j);
		}
		P(j,j) = 0.0;*/
	}
}

/**
* Replace the close to zero P with the DBL_MIN
**/
void nldr_tSNE::preventZero(tbox_DenseMatrix<double> &Ptemp,
						tbox_DenseMatrix<double> &P)
{
	int i,j,nRow;
	nRow = Ptemp.getNrows();
	
	P.resize(nRow,nRow);
	
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			if(i!=j)
			{
				P(i,j) = max(Ptemp(i,j),DBL_MIN);
				
			}
		}
	}

}

/**
* Exaggerate P
**/
void nldr_tSNE::exaggerateP(tbox_DenseMatrix<double> &P,
							tbox_DenseMatrix<double> &Plied)
{
	int i,j,nRow;
	
	nRow = P.getNrows();
	
	Plied.resize(nRow,nRow);
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			Plied(i,j) = d_multiP * P(i,j);
		}
	}
}

/**
* Initial solution
**/
void nldr_tSNE::initialY(tbox_DenseMatrix<double> &Y)
{
	int i,j,nRow;

	double std = 0.0001;
	
	nRow = Y.getNrows();

	// Normal distribution

	if(d_srand_flag == 1)
	{
		srand(time(NULL));
	}
	else
	{
		srand(0);
		//srand(2);
	}
	for(j=0;j<d_nDim;j++)
	{
		for(i=0;i<nRow;i++)
		{
			Y(i,j) = std*randnorm();
		}
	}
}

/**
* Normally distributed random numbers
* with zero mean and unit variance
**/
double nldr_tSNE::randnorm()
{


	double v1,v2,r_squared,factor;
	
	r_squared = 0.0;
	while(r_squared >= 1.0 || r_squared == 0.0)
	{
		// uniform numbers from -1 to 1
		v1 = (rand()%200)/100.0 -1.0;	
		v2 = (rand()%200)/100.0 -1.0;

		r_squared = v1*v1 + v2*v2;
	} 
	
	factor = sqrt(-2.0*log(r_squared)/r_squared);
	
	return v1*factor; 
}

/**
* Compute Q that has Student t-distribution
**/
// Column based
void nldr_tSNE::mappedQ(tbox_DenseMatrix<double> &Y,
						tbox_DenseMatrix<double> &numeratorM,
						tbox_DenseMatrix<double> &Q)
{
	int i,j,k,l,nRow,nCol;
	tbox_DenseMatrix<double> DsquaredY;
	double denominator;	

	nRow = Y.getNrows();	
	nCol = Y.getNcols();

	Q.resize(nRow,nRow);
	numeratorM.resize(nRow,nRow);
	
	if(nCol != d_nDim)
	{
		cout<<"Error in computing Q."<<endl;
		cout<<"Number of columns in Y is not the same as the number of required dimention."<<endl;
	}

	squaredDistMx(Y,DsquaredY);


	for(j=0;j<nRow;j++)
	{
		for(i=j+1;i<nRow;i++)
		{
			numeratorM(i,j) = 1.0/(1.0 + DsquaredY(i,j));
			numeratorM(j,i) = numeratorM(i,j);
		}
		numeratorM(j,j) = 0.0;
	}

	// Denominator
	denominator = 0.0;
	for(k=0;k<nRow;k++)
	{
		for(l=0;l<nRow;l++)
		{
			if(l!=k)
			{
				denominator += numeratorM(l,k);
			}
		}
	}	

	for(j=0;j<nRow;j++)
	{
		for(i=j+1;i<nRow;i++)
		{
			Q(i,j) = numeratorM(i,j)/denominator;
			Q(j,i) = Q(i,j);
		}
		Q(j,j) = 0.0;
	}
}

/**
* Compute the gradients
* Output: Y_grads
**/
void nldr_tSNE::gradients(tbox_DenseMatrix<double> &P,
						  tbox_DenseMatrix<double> &Q,
						  tbox_DenseMatrix<double> &numeratorM,
						  tbox_DenseMatrix<double> &Y,
						  tbox_DenseMatrix<double> &Y_grads)
{
	int i,j,t,nRow,nCol;
	tbox_DenseMatrix<double> L;
	vector<double> L_colSum;
	double Ltemp,gTemp;

	nRow = Y.getNrows();
	nCol = Y.getNcols();

	L.resize(nRow,nRow);
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			L(i,j) = (P(i,j) - Q(i,j))*numeratorM(i,j);
		}
	}

	// Compute sums along the column of L 
	L_colSum.clear();
	for(i=0;i<nRow;i++)
	{
		Ltemp = 0.0;
		for(j=0;j<nRow;j++)
		{
			Ltemp += L(i,j);
		}
		L_colSum.push_back(Ltemp);
	}

	Y_grads.resize(nRow,nCol);	
	
	for(t=0;t<nCol;t++)
	{
		for(i=0;i<nRow;i++)
		{
			gTemp = 0.0;
			for(j=0;j<nRow;j++)
			{
					gTemp += L(i,j)*Y(j,t);
			}
			Y_grads(i,t) = 4*(L_colSum[i]*Y(i,t) - gTemp);
		}
	}
}

/**
* Update the solution Y and Y_incs
**/
void nldr_tSNE::updateSolution(double eta,
							   double alpha,
							   tbox_DenseMatrix<double> &Adapt,
							   tbox_DenseMatrix<double> &Y_grads,
							   tbox_DenseMatrix<double> &Y_incs,
							   tbox_DenseMatrix<double> &Y)
{
	int i,j,nRow,nCol;
	double sTemp;

	nRow = Y.getNrows();
	nCol = Y.getNcols();

	if(d_adapt_flag == 1)
	{
		// Adaptive learning rate
		adaptiveLearning(Y_grads, Y_incs, Adapt);
	}

/*	for(i=0;i<nRow;i++)
	{
		for(j=0;j<nCol;j++)
		{
			cout<<Adapt(i,j)<<' ';
		}
		cout<<endl;
	}
	cout<<endl;
*/
	// Update Y
	for(j=0;j<nCol;j++)
	{
		for(i=0;i<nRow;i++)
		{
			Y_incs(i,j) = - eta*Adapt(i,j)*Y_grads(i,j) + alpha*Y_incs(i,j); //??????????
			
			Y(i,j) += Y_incs(i,j);
		}
	}

	// Adjust Y to have zero mean
	zeroMean(Y);
/*
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
	}	*/
}

/**
* Subtract off the mean
**/
void nldr_tSNE::zeroMean(tbox_DenseMatrix<double> &Y)
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

/**
* Adaptive learning rate
**/
void nldr_tSNE::adaptiveLearning(tbox_DenseMatrix<double> &Y_grads,
								 tbox_DenseMatrix<double> &Y_incs,
				 				 tbox_DenseMatrix<double> &Adapt)
{
	int i,j,nRow,nCol;
	
	nRow = Y_grads.getNrows();
	nCol = Y_grads.getNcols();

	for(j=0;j<nCol;j++)
	{
		for(i=0;i<nRow;i++)
		{
			// If same sign
			if( sameSign(Y_grads(i,j), Y_incs(i,j)) )
			{
				Adapt(i,j) = max(d_adapt_multi * Adapt(i,j), d_adapt_min);
			}
			// If not same sign
			else
			{
				Adapt(i,j) = max(Adapt(i,j) + d_adapt_plus, d_adapt_min);
			}
		}
	}
}


/**
* Compute constant part of cost function
**/
void nldr_tSNE::constCostFun(tbox_DenseMatrix<double> &P,
							 double &constCost)
{
	int i,j,nRow;

	nRow = P.getNrows();

	constCost = 0.0;
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			if(i != j)
			{
				constCost += P(i,j)*log(P(i,j));			
			}
		}
	}
}

/**
* Cost function
**/
void nldr_tSNE::costFunction(tbox_DenseMatrix<double> &P,
							 tbox_DenseMatrix<double> &Q,
							 double &constCost,
							 double &obj)
{
	int i,j,nRow;
	double sTemp;
	
	nRow = P.getNrows();
	
	sTemp = 0.0;
	for(j=0;j<nRow;j++)
	{
		for(i=0;i<nRow;i++)
		{
			if(i != j)
			{
				sTemp += P(i,j)*log(Q(i,j));
			}
		}
	}
	obj = constCost - sTemp;
}


/**
* Compute the norm of gradients
* Why is Y_grads changed after appying svd()??
**/
double nldr_tSNE::Gnorm(tbox_DenseMatrix<double> &Y_grads)
{
	int i,j,nRow,nCol;
	tbox_DenseMatrix<double> G,U,S,V;

	nRow = Y_grads.getNrows();
	nCol = Y_grads.getNcols();

	G.resize(nRow,nCol);
	U.resize(nRow,nRow);
	S.resize(nRow,nCol);
	V.resize(nCol,nCol);

	G = Y_grads;

	G.svd(U,S,V);

	return S(0,0);
}

bool nldr_tSNE::sameSign(double value1,
			  			 double value2)
{
	return value1*value2 > 0.0;
}


#ifndef LACKS_NAMESPACE
}
#endif






