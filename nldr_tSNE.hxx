//
// File:        nldr_tSNE.hxx
//

#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>

#ifndef included_nldr_tSNE
#define included_nldr_tSNE

#ifndef included_Sapphire_config
#include "Sapphire_config.hxx"
#endif

#ifndef included_nldr_NLDimRed
#include "nldr_NLDimRed.hxx"
#define included_nldr_NLDimRed
#endif

#ifndef included_di_InstanceArray
#include "di_InstanceArray.hxx"
#endif

#ifndef included_tbox_DenseMatrix
#include "tbox_DenseMatrix.hxx"
#endif

#ifndef included_nns_NNSearch
#include "nns_NNSearch.hxx"
#define included_nns_NNSearch
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
 * Class for the tSNE algorithm for nonlinear dimension reduction. 
 *
 * TODO: Yaru- pls add any references used in the implementation. You will
 * probably need to include other header files, such as the one for
 * eigenvalue decomp. You will also need to add functions to get and set
 * different variables used in the code or to provide the user relevant
 * information e.g. the eigenvalues.
 *
 */

class nldr_tSNE: public nldr_NLDimRed
{
public:

   /** 
   * Default Constructor
   **/
   nldr_tSNE();

   /** 
   * Copy Constructor
   **/
   nldr_tSNE(const nldr_tSNE &rhs);

   /** 
   * Constructor created with an appropriate nearest neighbor
   * algorithm. The nearest neighbor object is created in the calling
   * program and passed as a pointer to the base class.
   **/
   /*nldr_tSNE();*/


  /** 
   * Assignment operator
   **/
   nldr_tSNE &operator=(const nldr_tSNE &rhs);

  /** 
   * Destructor
   **/
   virtual ~nldr_tSNE();

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


   int getNumDimentions();


   void setPerplexity(double perp);

   double getPerplexity();

   void setMultiPflag(int multiP_flag);

   int getMultiPflag();

   void setMultiP(double multiP);

   double getMultiP();

   void setMultiPstopIter(int multiP_stop_iter);

   int getMultiPstopIter();

   void setTolerance(double tolerance);

   double getTolerance();

   void setOptTolerance(double tolerance);

   double getOptTolerance();

   void setNumIter(int numIter);

   int getNumIter();

   /** Initial Momentum **/
   void setAlphaInit(double alpha_init);

   double getAlphaInit();

   /** Momentum **/
   void setAlpha(double alpha);

   double getAlpha();

   /** The number of iterations that switches the alpha value **/
   void setAlphaSwitchIter(int alpha_switch_iter);

   int getAlphaSwitchIter();

   /** Inital learning rate **/
   void setEtaInit(double eta_init);

   double getEtaInit();

   /** Constant adaptive learning rate, the multiplication **/
   void setAdaptMulti(double adapt_multi);

   double getAdaptMulti();

   /** Constant adaptive learning rate, the plus **/
   void setAdaptPlus(double adapt_plus);

   double getAdaptPlus();

   /** Constant adaptive learning rate, the min **/
   void setAdaptMin(double adapt_min);

   double getAdaptMin();

   /** Constant adaptive learning rate, the flag **/
   void setAdaptFlag(int adapt_flag);

   int getAdaptFlag();

   /** Initial solution flag **/
   void setInitSolFlag(int initSol_flag);

   int getInitSolFlag();

   /** PCA flag **/
   void setPCAFlag(int pca_flag);

   int getPCAFlag();

   /** Random seed flag **/
   void setSrandFlag(int srand_flag);

   int getSrandFlag();

   /** Number of dimensions output from PCA **/
   void setPCAnDim(int pca_nDim);

   int getPCAnDim();

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
			  di_InstanceArray &InitSol);

private:

   /**
   * Normalize Input Data
   * X = (X - min(X)) / max(X)
   **/
   void normalize(tbox_DenseMatrix<double> &Input,
				  tbox_DenseMatrix<double> &X);


   /**
   * PCA
   **/
   void pca(tbox_DenseMatrix<double> &X,
			int &initDim,
			tbox_DenseMatrix<double> &Y);

   /**
   * Covariance Matrix
   **/
   void covarMx(tbox_DenseMatrix<double> &X,
				tbox_DenseMatrix<double> &CovarM);

   /**
   * Compute squared distance matrix in Euclidean
   * for Gaussian
   **/
   void squaredDistMx(tbox_DenseMatrix<double> &X,
							  tbox_DenseMatrix<double> &Dsquared);

   /**
   * Compute conditional probabilities
   * Identifies appropriate sigma's that satisfy the given perplexity
   * up to some tolerance.
   * Input: D.^2, perplexity and tolerance
   * Output: Conditional probabilities, beta=1/(2*(sigma^2))
   **/
   void d2p(tbox_DenseMatrix<double> &Dsquared,
			tbox_DenseMatrix<double> &Pcond,
			vector<double> &beta);

   /**
   * For a fixed point, idx,
   * Compute temporary conditional probabilities given beta=1/(2*(sigma^2))
   **/
   void condProb(tbox_DenseMatrix<double> &Dsquared,
				 int &idx,
				 double &beta_s,
				 vector<double> &PcondV);

   /**
   * For a fixed point, idx,
   * Compute Shannon Entropy, H(P_i), using joint probabilities
   **/
   void shannonEntropy(vector<double> &PcondV,
					   int &idx,
					   double &h);

   /**
   * Given conditional probabilities, apply tSNE and 
   * output new coordinates
   * If there is no initial solution, set initial solution flag = 0.
   **/
   void tSNEp(tbox_DenseMatrix<double> &Pcond,
			  tbox_DenseMatrix<double> &NewCoord,
			  tbox_DenseMatrix<double> &InitialSolution);

   /**
   * Symmetrize probabilities
   **/
   void symmetrizeP(tbox_DenseMatrix<double> &Pcond,
					tbox_DenseMatrix<double> &P);


   /**
   * Make sure P sum to one
   **/
   void PsumOne(tbox_DenseMatrix<double> &Ptemp,
		  	    tbox_DenseMatrix<double> &P);

   /**
   * Replace the close to zero P with the DBL_MIN
   **/
   void preventZero(tbox_DenseMatrix<double> &Ptemp,
					tbox_DenseMatrix<double> &P);

   /**
   * Exaggerate P
   **/
   void exaggerateP(tbox_DenseMatrix<double> &P,
					tbox_DenseMatrix<double> &Plied);

   /**
   * Initial solution
   **/
   void initialY(tbox_DenseMatrix<double> &Y);


   /**
   * Normally distributed random numbers
   * with zero mean and unit variance
   **/
   double randnorm();


   /**
   * Compute Q that has Student t-distribution
   **/
   void mappedQ(tbox_DenseMatrix<double> &Y,
				tbox_DenseMatrix<double> &numeratorM,
				tbox_DenseMatrix<double> &Q);

   /**
   * Compute the gradients
   **/
   void gradients(tbox_DenseMatrix<double> &P,
				  tbox_DenseMatrix<double> &Q,
				  tbox_DenseMatrix<double> &numeratorM,
				  tbox_DenseMatrix<double> &Y,
				  tbox_DenseMatrix<double> &Y_grads);

   /**
   * Update the solution Y and Y_incs
   **/
   void updateSolution(double eta,
					   double alpha,
					   tbox_DenseMatrix<double> &Adapt,
					   tbox_DenseMatrix<double> &Y_grads,
					   tbox_DenseMatrix<double> &Y_incs,
					   tbox_DenseMatrix<double> &Y);

   /**
   * Subtract off the mean
   **/
   void zeroMean(tbox_DenseMatrix<double> &Y);


   /**
   * Adaptive learning rate
   **/
   void adaptiveLearning(tbox_DenseMatrix<double> &Y_grads,
			 			 tbox_DenseMatrix<double> &Y_incs,
			 			 tbox_DenseMatrix<double> &Adapt);

   /**
   * Compute constant part of cost function
   **/
   void constCostFun(tbox_DenseMatrix<double> &P,
					 double &constCost);

   /**
   * Cost function
   **/
   void costFunction(tbox_DenseMatrix<double> &P,
					 tbox_DenseMatrix<double> &Q,
					 double &constCost,
					 double &obj);

   /**
   * Compute the norm of gradients
   **/
   double Gnorm(tbox_DenseMatrix<double> &Y_grads);

   bool sameSign(double value1,
			  	 double value2);

   /** Number of dimension desired at output **/
   int d_nDim;
   /** Perplexity **/
   double d_perp;
   /** Exaggeration flag, default is 0: not using it **/
   int d_multiP_flag;
   /** The constant used for exaggerating the probability, P **/
   double d_multiP;
   /** The number of terations to stop exaggerating **/
   int d_multiP_stop_iter;
   /** Tolerance on desired perplexity **/
   double d_tolerance;
   /** Tolerance in optimization **/
   double d_optTolerance;
   /** Number of iterations intended to use in the optimization algorithm **/
   int d_numIter;			
   /** Initial Momentum **/
   double d_alpha_init;     
   /** Momentum **/
   double d_alpha;          
   /** The number of iterations that switches the alpha value **/
   int d_alpha_switch_iter;	
   /** Inital learning rate **/
   double d_eta_init;		
   /** Constant for adaptive learning rate **/
   double d_adapt_multi;
   double d_adapt_plus;
   double d_adapt_min;
   int d_adapt_flag;
   int d_initSol_flag;
   int d_pca_flag;
   int d_srand_flag;
   int d_pca_nDim;
};


#ifndef LACKS_NAMESPACE
}
#endif

#endif






