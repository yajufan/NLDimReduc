//
// File:        nldr_NLDimRed.cxx
//


#ifndef included_nldr_NLDimRed
#include "nldr_NLDimRed.hxx"
#endif

#ifndef included_cfloat
#include <cfloat>
#endif

#ifndef included_map
#include <map>
#endif



/**
* Transform InstanceArray to DenseMatrix
**/
void nldr_NLDimRed::InstanceArray2DenseMatrix(const di_InstanceArray &Input,
											  tbox_DenseMatrix<double> &InputDMx)
{
	int i,j,nRow,nCol;
	int count,count2;
	nRow = Input.getNumInstances();
	nCol = Input.getNumFeatures();

	// loop for the columns
	count = 0;
	for(j=0;j<nCol;j++)
	{
		// Check selected features
		if((Input.getFeatureType(j) == di_kNumeric) &&
		   (Input.getIgnoreFlag(j) == false) )
		{
			count++;	
		}
	}

	//---------------------------------
	// Resize if data set is too large
	//---------------------------------
/*	if(nRow > 500)
	{
		cout<<"resize..."<<endl;
		nRow = 600;
	}	
*/
	// Resize DenseMatrix
	InputDMx.resize(nRow,count);
	count2 = 0;
	for(j=0;j<nCol;j++)
	{
		// Check selected features
		if((Input.getFeatureType(j) == di_kNumeric) &&
		   (Input.getIgnoreFlag(j) == false) )
		{
			for(i=0;i<nRow;i++)
			{
				InputDMx(i,count2) = Input[i][j];
			}
			count2 += 1;
		}
	}

}

double nldr_NLDimRed::residualVar(vector<double> &X,
				 			    vector<double> &Y)
{
	int i,j;
	double tempXY,tempXX,tempYY;
	double rVar;

	int n = X.size();
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

/**
* Wrap results will keep the ignored features.
* Then append mat at the end.
**/
void nldr_NLDimRed::wrapResults(const di_InstanceArray &input,
                                tbox_DenseMatrix<double> &mat,
                                di_InstanceArray &output)
{

   int sj, m;
   unsigned i, j, index;
   unsigned numInstances, numOrigFeatures, numOrigIgnored, numNewFeatures;
   di_InstanceInfo newinfo;
   di_InstanceInfo temp_info;
   tbox_Array1D<bool> selfeatures;
   char *stringX;


   di_InstanceInfo originfo = *input.getInstanceInfo();
   numOrigFeatures = input.getNumFeatures();
   numInstances = input.getNumInstances();
   m = mat.getNcols();

   selfeatures.resize(numOrigFeatures);


   // count and identify which features were ignored in the input InstanceArray
   numOrigIgnored = 0;
   for(j=0; j < numOrigFeatures; j++)
   {
      selfeatures[j] = false;

      if(originfo.getIgnoreFlag(j) == true ||
         originfo.getFeatureType(j) == di_kNominal)
      {
         numOrigIgnored++;
         selfeatures[j] = true;
      }
   }

   // reduced dimensions + presumably housekeeping features
   numNewFeatures = mat.getNcols() + numOrigIgnored;

   // create instance info with names and types for new (reduced)
   // variables
   temp_info.setNumFeatures(m);
   stringX = new char[10];
   for (j=0; j < (unsigned)m; j++)
   {
      temp_info.setFeatureType(j, di_kNumeric);
      temp_info.setIgnoreFlag(j, false);
      sprintf(stringX, "X%d", j);
      temp_info.setFeatureName(j, stringX);
   }
   delete [] stringX;

   newinfo.select(originfo, selfeatures);
   newinfo.appendInfo(temp_info);
   output.setInstanceInfo(newinfo);
   output.resize(numInstances);

   // copy features ignored from input array and features selected
   for(i=0; i < numInstances; i++)
   {
      index = 0;
      for(j=0; j < numOrigFeatures; j++)
      {
         if (originfo.getIgnoreFlag(j) == true || originfo.getFeatureType(j) == di_kNominal)
         {
            output[i][index] = input[i][j];
            index++;
         }
      }
      for(sj=0; sj < m; sj++)
      {
         output[i][index] = mat.getItem(i, sj);
         index++;
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(index == numNewFeatures);
#endif
      output[i].setClassLabel(input[i].getClassLabel());
      output[i].setID(input[i].getID());
      output[i].setWeight(input[i].getWeight());
   }
}

void nldr_NLDimRed::orderEigVals(tbox_DenseMatrix<double> &EigVec_temp,
				  				  vector<double> &eigVal_temp,
				  				  tbox_DenseMatrix<double> &EigVec,
				  				  vector<double> &eigVal)
{
	int i,nRow,count;
	multimap<double,int> lambda;
	multimap<double,int>::iterator iter;
	
//	map<double,int> lambda;
//	map<double,int>::iterator iter;

	nRow = eigVal_temp.size();

	lambda.clear();
	for(i=0;i<nRow;i++)
	{		
		lambda.insert( pair<double,int>(eigVal_temp[i],i) );
	}

	EigVec.resize(nRow,nRow);

	eigVal.clear();
	count = 0;
	for(iter=lambda.begin();iter != lambda.end();iter++)
	{
		eigVal.push_back(eigVal_temp[(*iter).second]);
		
		
		for(i=0;i<nRow;i++)
		{
			EigVec(i,count) = EigVec_temp(i,(*iter).second);
		}
		count++;
	}
}



void nldr_NLDimRed::vector2tboxDM(vector<double> &rowBasedMatrix,
								tbox_DenseMatrix<double> &columnBasedMatrix,
								int nRow,
								int nCol)
{
	columnBasedMatrix.resize(nRow,nCol);
	
	for(int i=0;i<nRow;i++)
	{
		for(int j=0;j<nCol;j++)
		{
			columnBasedMatrix.setItem(i,j,rowBasedMatrix[i*nCol+j]);
		}
	}
}

void nldr_NLDimRed::array2tboxDM(double* rowBasedMatrix,
								  tbox_DenseMatrix<double> &columnBasedMatrix,
								  int nRow,
								  int nCol)
{
	columnBasedMatrix.resize(nRow,nCol);
	
	for(int i=0;i<nRow;i++)
	{
		for(int j=0;j<nCol;j++)
		{
			columnBasedMatrix.setItem(i,j,rowBasedMatrix[i*nCol+j]);
		}
	}
}

void nldr_NLDimRed::tboxDM2vector(tbox_DenseMatrix<double> &columnBasedMatrix,
								  vector<double> &output_vector)
{
	int nRow,nCol;

	nRow = columnBasedMatrix.getNrows();
	nCol = columnBasedMatrix.getNcols();

	output_vector.clear();
	for(int i=0;i<nRow;i++)
	{
		for(int j=0;j<nCol;j++)
		{
			output_vector.push_back(columnBasedMatrix(i,j));
		}
	}


}

void nldr_NLDimRed::array2vector(double *input_array,
	     			 			 vector<double> &output_vector,
					 			 int nRow)
{
	int i;

	output_vector.clear();
	
	for(i=0;i<nRow;i++)
	{
		output_vector.push_back(input_array[i]);
	}
	
}

double nldr_NLDimRed::min(double value1,
						  double value2)
{
	if(value1 < value2)
		return value1;
	else
		return value2;
}

int nldr_NLDimRed::min(int value1,
					   int value2)
{
	if(value1 < value2)
		return value1;
	else
		return value2;
}


double nldr_NLDimRed::max(double value1,
						double value2)
{
	if(value1 > value2)
		return value1;
	else
		return value2;
}

double nldr_NLDimRed::mean(vector<double> &values)
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







