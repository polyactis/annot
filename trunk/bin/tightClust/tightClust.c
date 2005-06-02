#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "cluster.h"

const double NA=1E30;
const int maxStringLength=500;
void tab_fscanf(FILE* stream, char* format, char s[]);
void write_data(int nrow, int ncol, double** data, char** sampleName, char** geneID, char** annotation, char path[]);
void write_double_matrix(int nrow, int ncol, double** data, char path[]);
void write_int_matrix(int nrow, int ncol, int** data, char path[]);
void write_cluster_treeview(int nrow, int ncol, double** data, char** sampleName, char** geneID, char** annotation, int targetClustNum, int clusterid[], char path[], int emptyrows_betwn_cluster);
int which_minDist(int ncol, double** data, double** cdata, int** mask, int** cmask, 
	double weight[], int index1, int nclusters);
void calcAveComemMatrix(int nclusters, int nrow, int ncol, double** data, int** mask, 
	double weight[], int npass, int resampNum, double subSampPercnt, short** aveComemMatrix);
void select_candidateClust_One(short** aveComemMatrix, int nrow, int** clustIndex, int* nclusters, 
	int clustSize[], int resampNum);
void select_candidateClust_alpha(int nclusters, int nrow, int ncol, double** data, int** mask, 
	double weight[], int npass, int resampNum, double subSampPercnt, int topNum, 
	double alpha, int** clustIndex, int clustSize[]);
int extend_oneCandidateCluster(int** clustIndex, int clustSize[], int index, short** aveComemMatrix, int nrow, double alpha, int excludedObject[], int excludedObjectNum, int resampNum);
void extend_candidateClusters(int** clustIndex, int* pnclusters, int clustSize[], short** aveComemMatrix, int nrow, double alpha, int resampNum);
void get_topClusters(int* pnclusters, int** clustIndex, int clustSize[], int topNum, int nrow);
void tightClust(int targetClustNum, int nrow, int ncol, double** data, char** sampleName,
	char** geneID, char** annotation, int min_k, int max_k, int resampNum, int topNum, 
	int seqNum, double alpha, double beta, int* seed, double subSampPercnt, 
	double remainGenePercnt, int clusterid[], double**cdata, int** mask, 
	double weight[], int npass, char logFileName[], char tempClustFileName[]);
void find_single_tightClust(int nrow, int ncol, double** data, int min_k, int max_k, 
		int resampNum, int topNum, int seqNum, double alpha, double beta, int* seed, 
		double subSampPercnt, int clustIndex[], int* clustSize, double**cdata, int** mask, 
		double weight[], int npass, FILE* logStream);
double calc_IOU(int clustIndex1[], int clustSize1, int clustIndex2[], int clustSize2);
void combineChars(char first[], char second[], char combined[]);
void assignChars(char oldChars[], char newChars[]);
/*
//data  (input) double array
//The data array containing the second vector.

//mask  (input) int array
//This array which elements in data1 are missing. If mask1[i][j]==0, then data[i][j] is missing.
*/

//modify cluster.c: remark initran in kcluster so that same clustering result is obtained in different run 5/25/2004
//modify cluster.c:  #include<ranlib.h> is changed to #include "ranlib.h"  5/25/2004
//modify cluster.c kcluster: 1. Add return parameter "cmask". 
//					2. Remove memory allocation for "cmask"  
//                  3. copy "tcmask" to "cmask" when "cdata" replaced by "tcdata"  5/25/2004 
// modify cluster.c and cluster.h: 1. change "static double euclid (int n,..." to "double CALL euclid (int n,...".
//					2. add "euclid" function to "cluster.h"

int main(int argc, char* argv[])
{	
	FILE *stream, *logStream;
	char s[500];
	int nrow=10, ncol=10;
	int i,j;
	double** data;
	double* weight;
	int** mask;
	char** sampleName;
	char** geneID;
	char** annotation;
	int* clusterid;
	double** cdata;
	int** cmask;
	double error;
	//int ifound;
	int min_k=120, max_k=125;
	int resampNum=5, topNum=7, seqNum=2, npass=1, seed=12345, targetClustNum=3;
	double alpha=0.1, beta=0.6, subSampPercnt=0.7, remainGenePercnt=0.1;
	char workingDir[200];
	char dataFile[200];
	char fullDataFileName[400];
	char outputFileName[400];
	char buffer[10];
	char fullLogFileName[400];
	char fullTempClustFileName[400];
	char fullTempDataFileName[400];  //this is used for users to confirm that their data are read in correctly
	char tempString[200];
	//temp variables
	//short** aveComemMatrix;
	int nclusters;
	int** clustIndex;
	int* clustSize;
	
	assignChars(workingDir, argv[1]);
	assignChars(dataFile, argv[2]);
	targetClustNum=atoi(argv[3]);
	min_k=atoi(argv[4]);
	max_k=atoi(argv[5]);
	alpha=atof(argv[6]);
	beta=atof(argv[7]);
	topNum=atoi(argv[8]);
	seqNum=atoi(argv[9]);
	resampNum=atoi(argv[10]);
	subSampPercnt=atof(argv[11]);
	npass=atoi(argv[12]);
		
	combineChars(workingDir, dataFile, fullDataFileName);
	combineChars(workingDir, "temp_confirmData.txt", fullTempDataFileName);
	combineChars(workingDir, "temp_clust.txt", fullTempClustFileName);
	//outputFileName
	combineChars(workingDir, "tightClust_", outputFileName);
	//Note: _itoa and _gcvt only works in Visual C++ in Windows
	combineChars(outputFileName, argv[3], outputFileName);
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[4], outputFileName);
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[5], outputFileName);	
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[6], outputFileName);
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[7], outputFileName);
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[8], outputFileName);	
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[9], outputFileName);	
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[10], outputFileName);	
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[11], outputFileName);
	combineChars(outputFileName, "_", outputFileName);
	combineChars(outputFileName, argv[12], outputFileName);	
	combineChars(outputFileName, ".txt", outputFileName);

	combineChars(workingDir, "temp_log_", fullLogFileName);
	combineChars(fullLogFileName, argv[3], fullLogFileName);
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[4], fullLogFileName);
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[5], fullLogFileName);	
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[6], fullLogFileName);
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[7], fullLogFileName);
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[8], fullLogFileName);	
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[9], fullLogFileName);	
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[10], fullLogFileName);	
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[11], fullLogFileName);
	combineChars(fullLogFileName, "_", fullLogFileName);
	combineChars(fullLogFileName, argv[12], fullLogFileName);	
	combineChars(fullLogFileName, ".txt", fullLogFileName);


	//fseek( stream, 0L, SEEK_SET );
	if( (stream  = fopen( fullDataFileName, "r" )) == NULL )
      	printf( "The data file was not opened\n" );
	tab_fscanf( stream, "%s", s );
	nrow=atoi(s);
	tab_fscanf( stream, "%s", s );
	ncol=atoi(s);

	if( (logStream  = fopen( fullLogFileName, "w" )) == NULL )
		printf( "The log file was not opened\n" );
	printf("read in %d genes and %d samples\n\n", nrow, ncol);
	fprintf(logStream, "read in %d genes and %d samples\n\n", nrow, ncol);
	fclose(logStream);

	//temp
	//aveComemMatrix = (short**)malloc((size_t)nrow*sizeof(short*));
	//for(i=0; i<nrow; i++)
	//	aveComemMatrix[i] = (short*)malloc((size_t)nrow*sizeof(short));	
	
	clustIndex = (int**)malloc((size_t)nrow*sizeof(int*));
	for(i=0; i<nrow; i++)
		clustIndex[i] = (int*)malloc((size_t)nrow*sizeof(int));
	clustSize= (int*)malloc((size_t)nrow*sizeof(int));


	//arrange memory for sampleName, geneID and annotation
	data = (double**)malloc((size_t)nrow*sizeof(double*));
	for(i=0; i<nrow; i++)
		data[i] = (double*)malloc((size_t)ncol*sizeof(double));	
	mask = (int**)malloc((size_t)nrow*sizeof(int*));  
	for(i=0; i<nrow; i++)
		mask[i] = (int*)malloc((size_t)ncol*sizeof(int));	
	weight = (double*)malloc((size_t)ncol*sizeof(double)); 
	sampleName = (char**)malloc((size_t)ncol*sizeof(char*));
	for(i=0; i<ncol; i++)
		sampleName[i] = (char*)malloc((size_t)maxStringLength*sizeof(char));
	geneID = (char**)malloc((size_t)nrow*sizeof(char*));
	for(i=0; i<nrow; i++)
		geneID[i] = (char*)malloc((size_t)maxStringLength*sizeof(char));
	annotation = (char**)malloc((size_t)nrow*sizeof(char*));
	for(i=0; i<nrow; i++)
		annotation[i] = (char*)malloc((size_t)maxStringLength*sizeof(char));

	clusterid = (int*)malloc((size_t)nrow*sizeof(int));	
	cdata = (double**)malloc((size_t)400*sizeof(double*));
	for(i=0; i<400; i++)
		cdata[i] = (double*)malloc((size_t)ncol*sizeof(double));
	cmask = (int**)malloc((size_t)400*sizeof(int*));
	for(i=0; i<400; i++)
		cmask[i] = (int*)malloc((size_t)ncol*sizeof(int));
	
	/////////////////////////////////////////////////////
	for(i=0; i<ncol; i++)
		weight[i]=1;
	
	tab_fscanf( stream, "%s", s ); //read in string "geneID"
	tab_fscanf( stream, "%s", s ); //read in string "annotation"

	for(i=0; i<ncol; i++)
		tab_fscanf(stream, "%s", sampleName[i]);  //read in sampleName

	for(i=0; i<nrow; i++)
	{
		tab_fscanf(stream, "%s", geneID[i]);
		tab_fscanf(stream, "%s", annotation[i]);
		for(j=0; j<ncol; j++)
		{
			tab_fscanf(stream, "%s", s);
			if(s[0]=='\0')
				mask[i][j]=0;   //missing value
			else
			{
				data[i][j]=atof(s);	
				mask[i][j]=1;
			}
		}
	}
	fclose(stream);
	
	write_double_matrix(nrow, ncol, data, fullTempDataFileName);

	//select_candidateClust_alpha(3, nrow, ncol, data, mask, weight, 2, 10, 
	//	0.7, 7, 0, clustIndex, clustSize);

	tightClust(targetClustNum, nrow, ncol, data, sampleName, geneID, annotation, 
		min_k, max_k, resampNum, topNum, seqNum, alpha, beta, &seed, subSampPercnt, 
		remainGenePercnt, clusterid, cdata, mask, weight, npass, 
		fullLogFileName, fullTempClustFileName);

	//kcluster(25, nrow, ncol, data, mask, weight, 0, 1, 'a', 'e', clusterid, cdata, cmask, &error, &ifound);
	//printf("help\n");
	write_cluster_treeview(nrow, ncol, data, sampleName, geneID, annotation, targetClustNum+1, clusterid, outputFileName, 3);
	/*calcAveComemMatrix(25, nrow, ncol, data, mask, weight, 3, 5, 0.7, aveComemMatrix);
	printf("help\n");
	write_double_matrix(nrow, nrow, aveComemMatrix, "c:\\temp_aveComemMatrix.txt");
	selectClust_One(aveComemMatrix, nrow, clustIndex, &nclusters, clustSize);
	write_int_matrix(nclusters, nrow, clustIndex, "c:\\temp_clustIndex.txt");
	printf("nclusters: %d\n", nclusters);
	for(i=0; i<nclusters; i++)
		printf("%d\n", clustSize[i]);
	for(i=0; i<nrow; i++) clusterid[i]=100;
	for(i=0; i<100; i++)
	{
		for(j=0; j<clustSize[i]; j++)
		{
			clusterid[clustIndex[i][j]]=i;
		}
	}
	write_cluster_treeview(nrow, ncol, data, sampleName, geneID, annotation, 101, clusterid, "c:\\temp_clusteredData.txt", 3);*/
	//scanf("%d", &i);
	return 0;
}
// extend a cluster (the "index"th cluster) from all ACM (average comembership matrix) index equal 1 to ACM larger than or equal to 1-alpha
// clustIndex: cluster indexes. clustIndex[2][3] represents the 4th object in the 3rd cluster
// temp_nclusters: total number of candidate clusters
// clustSize: cluster size. clustSize[2] represents the size of the 3rd cluster
// aveComemMatrix: average comembership matrix  (nrow by nrow)
int extend_oneCandidateCluster(int** clustIndex, int clustSize[], int index, short** aveComemMatrix, int nrow, double alpha, int excludedObject[], int excludedObjectNum, int resampNum)
{
	int i,j,k,l, isFound;
	int temp_badPairNum, max_badPairNum=-1, max_badClustIndex;
	int isClusterExcluded=0;
	int* isObjectExcluded = (int*)malloc((size_t)nrow*sizeof(int));	

	for(i=0; i<nrow; i++)
		isObjectExcluded[i]=0;
	for(i=0; i<excludedObjectNum; i++)
		isObjectExcluded[excludedObject[i]]=1;

	if(index>0)
	{
		for(i=0; i<excludedObjectNum; i++)
		{
			if(clustIndex[index][0]==excludedObject[i])
				isClusterExcluded=1;
		}
	}

	if(isClusterExcluded==0)
	{
		// add to the cluster all the objects that have average comembership index larger than 1-alpha with the cluster members
		for(i=0; i<nrow; i++)
		{
			if(isObjectExcluded[i]==0)
			{
				if((aveComemMatrix[i][clustIndex[index][0]]+0.0)/resampNum>=1-alpha-0.0001 && (aveComemMatrix[i][clustIndex[index][0]]+0.0)/resampNum<1-0.0001)
				{
					clustIndex[index][clustSize[index]]=i;					
					clustSize[index]++;
				}
			}
		}
		
		// filter objects (one by one) that have the largest number of bad pairs (pairs with average comembership index less that 1-alpha within the cluster)
		do
		{
			//identify the object in the cluster ("max_badClustIndex") that has the largest number of pairs 
			//  that have average comembership index less than 1-alpha
			max_badPairNum=-1;
			for(i=0; i<clustSize[index]; i++)
			{
				temp_badPairNum=0;
				for(j=0; j<clustSize[index]; j++)
				{
					if((aveComemMatrix[clustIndex[index][i]][clustIndex[index][j]]+0.0)/resampNum<1-alpha-0.0001)
						temp_badPairNum++;
				}
				if(temp_badPairNum>max_badPairNum)
				{
					max_badClustIndex=i;
					max_badPairNum=temp_badPairNum;
				}
			}
			
			if(max_badPairNum>0)
			{
				for(i=max_badClustIndex; i<clustSize[index]-1; i++)
				{
					clustIndex[index][i]=clustIndex[index][i+1];
				}
				clustSize[index]=clustSize[index]-1;
			}
		}while(max_badPairNum>0);
	}

	free(isObjectExcluded);	
	
	if(isClusterExcluded==1)
		return(0);
	else return(1);
}


// extend all clusters from all ACM (average comembership matrix) index equal 1 to ACM larger than or equal to 1-alpha
// clustIndex: cluster indexes. clustIndex[2][3] represents the 4th object in the 3rd cluster
// nclusters: total number of candidate clusters
// clustSize: cluster size. clustSize[2] represents the size of the 3rd cluster
// aveComemMatrix: average comembership matrix  (nrow by nrow)
void extend_candidateClusters(int** clustIndex, int* pnclusters, int clustSize[], short** aveComemMatrix, int nrow, double alpha, int resampNum)
{
	int i,j, a, nclusters=*(pnclusters);
	int* isEffectiveClust = (int*)malloc((size_t)nclusters*sizeof(int));	
	int* excludedObject = (int*)malloc((size_t)nrow*sizeof(int));	
	int excludedObjectNum=0;
	int temp_nclusters=nclusters;
	int** temp_clustIndex;
	int* temp_clustSize = (int*)malloc((size_t)temp_nclusters*sizeof(int));	
	
	temp_clustIndex = (int**)malloc((size_t)temp_nclusters*sizeof(int*));
	for(i=0; i<temp_nclusters; i++)
	{
		temp_clustSize[i]=clustSize[i];
		temp_clustIndex[i]=(int*)malloc((size_t)nrow*sizeof(int));	
		for(j=0; j<clustSize[i]; j++)
			temp_clustIndex[i][j]=clustIndex[i][j];
	}

	for(i=0; i<temp_nclusters; i++)
		isEffectiveClust[i]=0;

	for(i=0; i<temp_nclusters; i++)
	{
		isEffectiveClust[i]=extend_oneCandidateCluster(temp_clustIndex, temp_clustSize, i, aveComemMatrix, nrow, alpha, excludedObject, excludedObjectNum, resampNum);
		if(isEffectiveClust[i]==1)
		{
			for(j=0; j<temp_clustSize[i]; j++)
				excludedObject[excludedObjectNum+j]=temp_clustIndex[i][j];
			excludedObjectNum += temp_clustSize[i];
		}
	}

	nclusters=0;
	for(i=0; i<temp_nclusters; i++)
	{
		if(isEffectiveClust[i]==1)
		{
			clustSize[nclusters]=temp_clustSize[i];
			for(j=0; j<temp_clustSize[i]; j++)
			{
				clustIndex[nclusters][j]=temp_clustIndex[i][j];
			}
			nclusters++;
		}
	}

	*(pnclusters)=nclusters;

	free(isEffectiveClust);	
	free(excludedObject);	
	for(i=0; i<temp_nclusters; i++)
		free(temp_clustIndex[i]);
	free(temp_clustIndex);
	free(temp_clustSize);	
}

void select_candidateClust_alpha(int nclusters, int nrow, int ncol, double** data, int** mask, 
	double weight[], int npass, int resampNum, double subSampPercnt, int topNum, 
	double alpha, int** clustIndex, int clustSize[])
{
	short** aveComemMatrix;
	int i, j;
	int temp_nclusters;  //number of clusters selected by select_candidateClust_One
		
	aveComemMatrix = (short**)malloc((size_t)nrow*sizeof(short*));
	for(i=0; i<nrow; i++)
		aveComemMatrix[i] = (short*)malloc((size_t)nrow*sizeof(short));	

	calcAveComemMatrix(nclusters, nrow, ncol, data, mask, weight, npass, resampNum, 
		subSampPercnt, aveComemMatrix);
	//
	//write_double_matrix(nrow, nrow, aveComemMatrix, "c:\\temp_aveComemMatrix.txt");
//printf("111\n");
	select_candidateClust_One(aveComemMatrix, nrow, clustIndex, &temp_nclusters, clustSize, resampNum);

	//
	//write_int_matrix(temp_nclusters, nrow, clustIndex, "c:\\temp_clustOneIndex.txt");

	//if(temp_nclusters<topNum)
	//	printf("Warning: less than %d clusters selected by select_candidateClust_One\n", topNum);
//printf("222\n");
	extend_candidateClusters(clustIndex, &temp_nclusters, clustSize, aveComemMatrix, nrow, alpha, resampNum);   
	//
	//write_int_matrix(temp_nclusters, nrow, clustIndex, "c:\\temp_clustAlphaIndex.txt");
//printf("333\n");
	get_topClusters(&temp_nclusters, clustIndex, clustSize, topNum, nrow);
//printf("444\n");
	for(i=0; i<nrow; i++)
		free(aveComemMatrix[i]);
	free(aveComemMatrix);
}

// get the top "topNum" clusters with largest size
// clustIndex: an nrow by nrow matrix. 
// nclusters: an integer, the number of clusters with all pairwise stability strength larger than 1-alpha
// clustSize: an array of size nrow. The first "nclusters" numbers are effective. 
// topNum: 
// clustIndex, nclusters and clustSize will be updated. nclusters becomes topNum
// if nclusters<topNum, empty clusters will be assigned to clusters nclusters to topNum
void get_topClusters(int* pnclusters, int** clustIndex, int clustSize[], int topNum, int nrow)
{
	int i, j, nclusters=*pnclusters;
	int* remainIndex = (int*)malloc((size_t)nclusters*sizeof(int));	
	int* topIndex = (int*)malloc((size_t)topNum*sizeof(int));	
	int remainSize=nclusters;
	int max_clustSize, max_clust;
	int** temp_clustIndex;
	int* temp_clustSize=(int*)malloc((size_t)nclusters*sizeof(int));	

	temp_clustIndex=(int**)malloc((size_t)nclusters*sizeof(int*));
	for(i=0; i<nclusters; i++)
		temp_clustIndex[i]=(int*)malloc((size_t)nrow*sizeof(int));


	if(nclusters>topNum)
	{
		//get indexes of the "topNum" clusters with largest size
		for(i=0; i<nclusters; i++)
			remainIndex[i]=i;

		for(i=0; i<topNum; i++)
		{
			max_clustSize=0;
			for(j=0; j<remainSize; j++)
			{
				if(clustSize[remainIndex[j]]>max_clustSize)
				{
					max_clustSize=clustSize[remainIndex[j]];
					max_clust=j;
					topIndex[i]=remainIndex[j];
				}
			}
			for(j=max_clust; j<remainSize-1; j++)
			{
				remainIndex[j]=remainIndex[j+1];
			}
			remainSize--;			
		}		

		for(i=0; i<nclusters; i++)
		{
			temp_clustSize[i]=clustSize[i];
			for(j=0; j<clustSize[i]; j++)
				temp_clustIndex[i][j]=clustIndex[i][j];
		}

		for(i=0; i<topNum; i++)
		{
			clustSize[i]=temp_clustSize[topIndex[i]];
			for(j=0; j<temp_clustSize[topIndex[i]]; j++)
				clustIndex[i][j]=temp_clustIndex[topIndex[i]][j];
		}

		*pnclusters=topNum;
		//for(i=0; i<topNum; i++)
		//	printf("top Index: %d  %d\n", topIndex[i], clustSize[topIndex[i]]);	
	}
	else
	{
		for(i=nclusters; i<topNum; i++)
			clustSize[i]=0;
		*pnclusters=topNum;
	}
	free(remainIndex);	
	free(topIndex);	
	for(i=0; i<nclusters; i++)
		free(temp_clustIndex[i]);
	free(temp_clustIndex);
	free(temp_clustSize);	
}



// identify clusters such that all pairwise average comembership indexes within cluster equal one
// aveComemMatrix: nrow*nrow matrix The average comembership index of pair i and j
// 
// return parameters:
// clustIndex: nrow*nrow matrix. 
// nclusters: an integer, the number of clusters with all pairwise average comembership indexes equal one
// clustSize: an array of size nrow. The first "nclusters" numbers are effective. e.g. clustSize[2] representing length of cluster 3
void select_candidateClust_One(short** aveComemMatrix, int nrow, int** clustIndex, int* pnclusters, 
					 int clustSize[], int resampNum)
{
	int i,j,k, isFound;

	int temp_nclusters;
	int** temp_clustIndex;
	int* temp_clustSize= (int*)malloc((size_t)nrow*sizeof(int));
	temp_clustIndex = (int**)malloc((size_t)nrow*sizeof(int*));
	for(i=0; i<nrow; i++)
		temp_clustIndex[i] = (int*)malloc((size_t)nrow*sizeof(int));
	
	for(i=0; i<nrow; i++)
	{
		temp_clustIndex[i][0]=i;
		temp_clustSize[i]=1;
	}
	temp_nclusters=nrow;

	*(pnclusters)=0;
	
	for(i=nrow-1; i>=0; i--)
	{
		if(temp_clustSize[i]==1)
		{
			isFound=0;
			j=0;
			while(isFound==0 && j<i)
			{
				if((aveComemMatrix[i][j]+0.0)/resampNum>=0.9999)
				{
					temp_clustIndex[j][temp_clustSize[j]]=i;
					temp_clustSize[j]++;
					isFound=1;
				}
				j++;
			}
		}
		else
		{
			clustSize[*(pnclusters)]=temp_clustSize[i];
			for(j=0; j<clustSize[*(pnclusters)]; j++)
				clustIndex[*(pnclusters)][j]=temp_clustIndex[i][j];
			*pnclusters=*pnclusters+1;
		}
	}

	for(i=0; i<nrow; i++) free(temp_clustIndex[i]);
	free(temp_clustIndex);
	free(temp_clustSize);
}

// modified fscanf: read until \t or \n
void tab_fscanf(FILE* stream, char *format, char s[])
{
	char temp;
	int sLength=0;
	fscanf(stream, "%c", &temp);
	while(temp!='\t' && temp!='\n')
	{
		s[sLength]=temp;
		sLength++;
		fscanf(stream, "%c", &temp);
	}
	s[sLength]='\0';
} 

// write out matrix (majorly for debug)
void write_double_matrix(int nrow, int ncol, double** data, char path[])
{
	FILE* tempStream;
	int i,j;
	if( (tempStream  = fopen( path, "w+" )) == NULL )
      printf( "The file matrix data was not opened\n" );
	fprintf(tempStream, "%d\t%d\n", nrow, ncol);
	for(i=0; i<nrow; i++)
	{
		for(j=0; j<ncol; j++)
			fprintf(tempStream, "%f\t", data[i][j]);
		fprintf(tempStream, "\n");
	}
	fclose(tempStream);
}

void write_int_matrix(int nrow, int ncol, int** data, char path[])
{
	FILE* tempStream;
	int i,j;
	if( (tempStream  = fopen( path, "w+" )) == NULL )
      printf( "The file matrix data was not opened\n" );
	fprintf(tempStream, "%d\t%d\n", nrow, ncol);
	for(i=0; i<nrow; i++)
	{
		for(j=0; j<ncol; j++)
			fprintf(tempStream, "%d\t", data[i][j]);
		fprintf(tempStream, "\n");
	}
	fclose(tempStream);
}

//write data matrix
void write_data(int nrow, int ncol, double** data, char** sampleName, char** geneID, char** annotation, char path[])
{
	FILE* tempStream;
	int k, i,j;
	if( (tempStream  = fopen( path, "w+" )) == NULL )
      printf( "The file 'data' was not opened\n" );
	fprintf(tempStream, "%d\t%d\n", nrow, ncol);

	fprintf(tempStream, "%s\t%s\t", "geneID", "annotation");
	for(i=0; i<ncol; i++)
	{
		if(i!=ncol-1)
			fprintf(tempStream, "%s\t", sampleName[i]);
		else fprintf(tempStream, "%s\n", sampleName[i]);
	}

	for(i=0; i<nrow; i++)
	{
			fprintf(tempStream, "%s\t%s\t", geneID[i], annotation[i]);
			for(j=0; j<ncol; j++)
			{
				if(j!=ncol-1)
					fprintf(tempStream, "%f\t", data[i][j]);
				else fprintf(tempStream, "%f\n", data[i][j]);
			}		
	}
	fclose(tempStream);
}

//write cluster result in treeview format
void write_cluster_treeview(int nrow, int ncol, double** data, char** sampleName, char** geneID, char** annotation, int targetClustNum, int clusterid[], char path[], int emptyrows_betwn_cluster)
{
	FILE* tempStream;
	int k, i,j;
	if( (tempStream  = fopen( path, "w+" )) == NULL )
      printf( "The file 'data' was not opened\n" );
	fprintf(tempStream, "%d\t%d\n", nrow, ncol);

	fprintf(tempStream, "%s\t%s\t", "geneID", "annotation");
	for(i=0; i<ncol; i++)
	{
		if(i!=ncol-1)
			fprintf(tempStream, "%s\t", sampleName[i]);
		else fprintf(tempStream, "%s\n", sampleName[i]);
	}

	for(k=0; k<targetClustNum; k++)
	{
		for(i=0; i<nrow; i++)
		{
			if(clusterid[i]==k)
			{
				fprintf(tempStream, "%s\t%s\t", geneID[i], annotation[i]);
				for(j=0; j<ncol; j++)
				{
					if(j!=ncol-1)
						fprintf(tempStream, "%f\t", data[i][j]);
					else fprintf(tempStream, "%f\n", data[i][j]);
				}
			}
		}
		for(i=0; i<emptyrows_betwn_cluster; i++)
			fprintf(tempStream, "NONE\n");
	}
	fclose(tempStream);
}


// calculate which cluster is closest to the chosen data row
// data: nrow*ncol
// cdata: nclusters(ncol
// index1: index to the chosen data row in "data"
int which_minDist(int ncol, double** data, double** cdata, int** mask, int** cmask, 
	double weight[], int index1, int nclusters)
{
	int i;
	int minIndex;
	double minDist=1E30, dist;

	for(i=0; i<nclusters; i++)
	{
		dist=euclid(ncol, data, cdata, mask, cmask, weight, index1, i, 0);
		if(dist<minDist)
		{
			minIndex=i;
			minDist=dist;
		}
	}
	return(minIndex);
}


// calculate average comembership matrix in "resampNum" times of resampling
void calcAveComemMatrix(int nclusters, int nrow, int ncol, double** data, int** mask, 
	double weight[], int npass, int resampNum, double subSampPercnt, short** aveComemMatrix)
{
	int subSampNum=nrow*subSampPercnt;
	int i,j,k;
	double** subData;
	double** subCdata;
	int** subCmask;
	int** subMask;
	int* subClusterid=(int*)malloc((size_t)subSampNum*sizeof(int)); 
	int* allClusterid=(int*)malloc((size_t)nrow*sizeof(int)); 
	int* originalIndex=(int*)malloc((size_t)nrow*sizeof(int)); 
	double error; int ifound;

	subData = (double**)malloc((size_t)subSampNum*sizeof(double*)); 
	for(i=0; i<subSampNum; i++)
		subData[i] = (double*)malloc((size_t)ncol*sizeof(double)); 
	subCdata = (double**)malloc((size_t)nclusters*sizeof(double*)); 
	for(i=0; i<nclusters; i++)
		subCdata[i] = (double*)malloc((size_t)ncol*sizeof(double)); 
	subMask= (int**)malloc((size_t)subSampNum*sizeof(int*)); 
	for(i=0; i<subSampNum; i++)
		subMask[i] = (int*)malloc((size_t)ncol*sizeof(int)); 
	subCmask= (int**)malloc((size_t)nclusters*sizeof(int*)); 
	for(i=0; i<nclusters; i++)
		subCmask[i] = (int*)malloc((size_t)ncol*sizeof(int)); 

	for(i=0; i<nrow; i++)
		for(j=0; j<nrow; j++)
			aveComemMatrix[i][j]=0;
	for(i=0; i<nrow; i++)
		originalIndex[i]=i;

	for(k=0; k<resampNum; k++)
	{
		//printf("the %dth resampling!\n", k);
		genprm(originalIndex, nrow);
		for(i=0; i<subSampNum; i++)
		{
			for(j=0; j<ncol; j++)
			{
				subData[i][j]=data[originalIndex[i]][j];
				subMask[i][j]=mask[originalIndex[i]][j];
			}
		}

		if(nclusters<subSampNum)
			kcluster(nclusters, subSampNum, ncol, subData, subMask, weight, 0, npass, 'a', 'e', subClusterid, subCdata, subCmask, &error, &ifound);
		else {printf("Error: nclusters larger than subSampNum!!"); scanf("%d", &i);}

		for(i=0; i<nrow; i++)
			allClusterid[i]=which_minDist(ncol, data, subCdata, mask, subCmask, weight, i, nclusters);

		for(i=0; i<nrow; i++)
		{
			for(j=0; j<nrow; j++)
			{
				if(allClusterid[i]==allClusterid[j])
					aveComemMatrix[i][j] += 1;
			}
		}
	}

	//free memory
	for(i=0; i<subSampNum; i++) free(subData[i]);
	for(i=0; i<subSampNum; i++) free(subMask[i]);
	for(i=0; i<nclusters; i++) free(subCdata[i]);
	for(i=0; i<nclusters; i++) free(subCmask[i]);
	free(subData);
	free(subCdata);
	free(subMask);
	free(subCmask);
	free(subClusterid);
	free(allClusterid);
	free(originalIndex);
}

// tight cluster
// Following Wing's suggestion to start from seqNum 3 and automatically switch to 2 if cluster not found
// clusterNum: total number of clusters to be identified
// topNum: the number of top clusters (in size) to be examined in a sequence of excessClustNum
// seqNum: The number that clusters are stablized (with intersectOverUnion larger than beta) in a sequence of excessClustNum
// beta: The threshold for intersectOverUnion (takes care of stability)
// logFileName: the file name to save screen outputs
// tempClustFileName: the file name to save currently obtained tight clusters
// return parameter:  clusterid
void tightClust(int targetClustNum, int nrow, int ncol, double** data, char** sampleName,
	char** geneID, char** annotation, int min_k, int max_k, int resampNum, int topNum, 
	int seqNum, double alpha, double beta, int* seed, double subSampPercnt, 
	double remainGenePercnt, int clusterid[], double**cdata, int** mask, 
	double weight[], int npass, char logFileName[], char tempClustFileName[])
{
	FILE* logStream;
	FILE* tempStream;
	int k, i, j, remainNrow, oldRemainNrow;
	double** remainData;
	double** tempData;
	//char** remainGeneID;
	//char** remainAnnotation;
	int** remainMask;
	int** tempMask;
	int* remainOriginalIndex;
	int* tempOriginalIndex;
	int* tempIntRow;
	int* singleClustIndex;
	int singleClustSize;
	int notFoundTimes;
	int temp_remainNrow, isInCluster;

	for(i=0; i<nrow; i++)
		clusterid[i]=targetClustNum;



	tempIntRow = (int*)malloc((size_t)nrow*sizeof(int)); 
	singleClustIndex = (int*)malloc((size_t)nrow*sizeof(int)); 
	remainData = (double**)malloc((size_t)nrow*sizeof(double*));
	tempData = (double**)malloc((size_t)nrow*sizeof(double*));
	remainMask = (int**)malloc((size_t)nrow*sizeof(int*));
	tempMask = (int**)malloc((size_t)nrow*sizeof(int*));
	//remainGeneID = (char**)malloc((size_t)nrow*sizeof(char*));
	//remainAnnotation = (char**)malloc((size_t)nrow*sizeof(char*));
	for(i=0; i<nrow; i++)
	{
		remainMask[i] = (int*)malloc((size_t)ncol*sizeof(int));
		tempMask[i] = (int*)malloc((size_t)ncol*sizeof(int));
		remainData[i] = (double*)malloc((size_t)ncol*sizeof(double));
		tempData[i] = (double*)malloc((size_t)ncol*sizeof(double));
		//remainGeneID[i] = (char*)malloc((size_t)maxStringLength*sizeof(char));
		//remainAnnotation[i] = (char*)malloc((size_t)maxStringLength*sizeof(char));
	}
	remainOriginalIndex = (int*)malloc((size_t)nrow*sizeof(int)); 
	tempOriginalIndex = (int*)malloc((size_t)nrow*sizeof(int)); 

	for(i=0; i<nrow; i++)
	{
		remainOriginalIndex[i]=i;
		for(j=0; j<ncol; j++)
		{
			remainData[i][j]=data[i][j];
			remainMask[i][j]=mask[i][j];
		}
		for(j=0; j<maxStringLength; j++)
		{
			//remainGeneID[i][j]=geneID[i][j];				
			//remainAnnotation[i][j]=annotation[i][j];
		}
	}
	remainNrow=nrow;


	if( (logStream  = fopen( logFileName, "a" )) == NULL )
		printf( "The log file was not opened\n" );

	if( (tempStream  = fopen( tempClustFileName, "w+" )) == NULL )
		printf( "The temp clust file was not opened\n" );

	fprintf(tempStream, "%d\t%d\n", nrow, ncol);
	fprintf(tempStream, "%s\t%s\t", "geneID", "annotation");
	for(i=0; i<ncol; i++)
	{
		if(i!=ncol-1)
			fprintf(tempStream, "%s\t", sampleName[i]);
		else fprintf(tempStream, "%s\n", sampleName[i]);
	}
	fclose(tempStream);

	for(i=0; i<targetClustNum; i++)
	{
		printf("Looking for the %dth tight cluster:\n", i+1);
		fprintf(logStream, "Looking for the %dth tight cluster:\n", i+1);
		notFoundTimes=0;
		do{
			find_single_tightClust(remainNrow, ncol, remainData, min_k, max_k, 
				resampNum, topNum, seqNum, alpha, beta, seed, subSampPercnt, 
				singleClustIndex, &singleClustSize, cdata, remainMask, weight, npass, logStream);
			notFoundTimes++;
		}while(notFoundTimes<3 && singleClustSize==0);

		if(seqNum>2 && singleClustSize==0) 
		{
			seqNum--;
			notFoundTimes=0;
			do{
				find_single_tightClust(remainNrow, ncol, remainData, min_k, max_k, 
					resampNum, topNum, seqNum, alpha, beta, seed, subSampPercnt, 
					singleClustIndex, &singleClustSize, cdata, remainMask, weight, npass, logStream);
				notFoundTimes++;
			}while(notFoundTimes<3 && singleClustSize==0);
		}

		//for(j=0; j<singleClustSize; j++)
		//	printf("%d ", remainOriginalIndex[singleClustIndex[j]]);
		//printf("\n\n");

		if(singleClustSize>0)
		{
			// output logFile and tempClustFile
			printf("The %dth cluster found!!\nClustSize: %d. Remaining number of genes: %d\n\n", i+1, singleClustSize, remainNrow-singleClustSize);
			fprintf(logStream, "The %dth cluster found!!\nClustSize: %d. Remaining number of genes: %d\n\n", i+1, singleClustSize, remainNrow-singleClustSize);
			if( (tempStream  = fopen( tempClustFileName, "a" )) == NULL )		
				printf( "The temp clust file was not opened\n" );
			for(j=0; j<singleClustSize; j++)
			{
				fprintf(tempStream, "%s\t", geneID[remainOriginalIndex[singleClustIndex[j]]]);
				fprintf(tempStream, "%s\t", annotation[remainOriginalIndex[singleClustIndex[j]]]);					
				for(k=0; k<ncol; k++)
				{
					if(k!=ncol-1)
						fprintf(tempStream, "%f\t", data[remainOriginalIndex[singleClustIndex[j]]][k]);
					else fprintf(tempStream, "%f\n", data[remainOriginalIndex[singleClustIndex[j]]][k]);
				}
			}
			fprintf(tempStream, "NONE\nNONE\nNONE\n");
			fclose(tempStream);

		//update remainNrow, remainData, remainMask, remainOriginalIndex
			//copy remainData,  remainMask, remainOriginalIndex to tempData, tempMask, tempOriginalIndex
			for(j=0; j<remainNrow; j++)
			{
				for(k=0; k<ncol; k++)
				{
					tempData[j][k]=remainData[j][k];
					tempMask[j][k]=remainMask[j][k];
				}
				tempOriginalIndex[j]=remainOriginalIndex[j];
			}

			//assign returned parameter, clusterid
			for(j=0; j<singleClustSize; j++)
				clusterid[remainOriginalIndex[singleClustIndex[j]]]=i;

			//take out the identified tight cluster
			temp_remainNrow=0;
			for(j=0; j<remainNrow; j++)
			{
				isInCluster=0;
				for(k=0; k<singleClustSize; k++)
				{
					if(singleClustIndex[k]==j)
						isInCluster=1;
				}
				
				// add to remainData, remainMask and remainOriginalIndex if not in the cluster
				if(isInCluster==0)
				{
					for(k=0; k<ncol; k++)
					{
						remainData[temp_remainNrow][k]=tempData[j][k];
						remainMask[temp_remainNrow][k]=tempMask[j][k];						
					}
					remainOriginalIndex[temp_remainNrow]=tempOriginalIndex[j];
					temp_remainNrow++;
				}
			}
			if(remainNrow!=temp_remainNrow+singleClustSize) printf("Error here!");
			remainNrow=temp_remainNrow;			
		}
		else
		{
			printf("Warning: the %dth cluster not found!\n", i);
			fprintf(logStream, "Warning: the %dth cluster not found!\n", i);
		}

		/*for(j=0; j<remainNrow; j++)
			printf("%d ", remainOriginalIndex[j]);
		printf("\n");
		scanf("%d", &j);*/
		min_k--;
		max_k--;	
	}
	
	for(i=0; i<nrow; i++)
	{
		free(remainData[i]);
		free(tempData[i]);
		free(remainMask[i]);
		free(tempMask[i]);
	}
	free(remainData);
	free(tempData);
	free(remainMask);
	free(tempMask);
	free(remainOriginalIndex);
	free(tempOriginalIndex);
	free(tempIntRow);
	free(singleClustIndex);
	
}

// find a single tight cluster
// return parameters:
// clustIndex and clustSize
void find_single_tightClust(int nrow, int ncol, double** data, int min_k, int max_k, 
		int resampNum, int topNum, int seqNum, double alpha, double beta, int* seed, 
		double subSampPercnt, int clustIndex[], int* clustSize, double**cdata, int** mask, 
		double weight[], int npass, FILE* logStream)
{
	int isTightClustFound, i,j,k, nclusters;
	double max_IOU, temp_IOU;
	int max_IOUIndex;
	int** old_clustIndex;
	int* old_clustSize=(int*)malloc((size_t)nrow*sizeof(int));
	int** new_clustIndex;
	int* new_clustSize=(int*)malloc((size_t)nrow*sizeof(int));
	int* old_cumu_seqNum=(int*)malloc((size_t)topNum*sizeof(int));
	int* new_cumu_seqNum=(int*)malloc((size_t)topNum*sizeof(int));


	old_clustIndex = (int**)malloc((size_t)nrow*sizeof(int*));
	new_clustIndex = (int**)malloc((size_t)nrow*sizeof(int*));
	for(i=0; i<nrow; i++)
	{
		old_clustIndex[i] = (int*)malloc((size_t)nrow*sizeof(int));
		new_clustIndex[i] = (int*)malloc((size_t)nrow*sizeof(int));
	}

	isTightClustFound=0; i=0; *clustSize=0;
	do
	{
		nclusters=min_k+i;
		printf("nclusters: %d\n", nclusters);
		fprintf(logStream, "nclusters: %d\n", nclusters);
		//printf("nclusters: %d\n", nclusters);
		if(i==0) 
		{
			select_candidateClust_alpha(nclusters, nrow, ncol, data, mask, 
				weight, npass, resampNum, subSampPercnt, topNum, alpha, 
				old_clustIndex, old_clustSize);
			for(j=0; j<topNum; j++)
				old_cumu_seqNum[j]=0;

			/*for(j=0; j<topNum; j++)
			{
				printf("round %d: ", i);
				for(k=0; k<old_clustSize[j]; k++)
					printf("%d ", old_clustIndex[j][k]);
				printf("\n");
			}*/
		}
		else
		{
			select_candidateClust_alpha(nclusters, nrow, ncol, data, mask, 
				weight, npass, resampNum, subSampPercnt, topNum, alpha, 
				new_clustIndex, new_clustSize);

			/*for(j=0; j<topNum; j++)
			{
				printf("round %d: ", i);
				for(k=0; k<new_clustSize[j]; k++)
					printf("%d ", new_clustIndex[j][k]);
				printf("\n");
			}*/

			for(j=0; j<topNum; j++)
			{
				max_IOU=0;
				for(k=0; k<topNum; k++)
				{
					if(new_clustSize[j]>0 && old_clustSize[k]>0)
						temp_IOU=calc_IOU(new_clustIndex[j], new_clustSize[j], old_clustIndex[k], old_clustSize[k]);
					else temp_IOU=-1;
					if(temp_IOU>max_IOU)
					{
						max_IOU=temp_IOU;
						max_IOUIndex=k;
					}
				}
				

				if(max_IOU>beta)
					new_cumu_seqNum[j]=old_cumu_seqNum[max_IOUIndex]+1;
				else new_cumu_seqNum[j]=0;

				//printf("max_IOU: %f  max_IOUIndex: %d  old_cumu: %d  new_cumu: %d\n", max_IOU, max_IOUIndex, old_cumu_seqNum[j], new_cumu_seqNum[j]);

				if(new_cumu_seqNum[j]==seqNum-1)
				{
					isTightClustFound=1;
					if(old_clustSize[max_IOUIndex]>*clustSize)
					{
						*clustSize=old_clustSize[max_IOUIndex];
						for(k=0; k<*clustSize; k++)
							clustIndex[k]=old_clustIndex[max_IOUIndex][k];
					}
				}
			}

			//copy new_clustIndex, new_clustSize, new_cumu_seqNum to old
			for(j=0; j<topNum; j++)
			{
				old_clustSize[j]=new_clustSize[j];
				old_cumu_seqNum[j]=new_cumu_seqNum[j];
				for(k=0; k<nrow; k++)
					old_clustIndex[j][k]=new_clustIndex[j][k];
			}
		}
		i++;
		//scanf("%d", &j);
	}while(i<=max_k-min_k && isTightClustFound==0);

	/*printf("isFound: %d\n", isTightClustFound);
	printf("tight Clust: ");
	for(k=0; k<*clustSize; k++)
		printf("%d ", clustIndex[k]);
	printf("\n");
	printf("clustSize: %d\n", *clustSize);
	scanf("%d", &j);*/

	for(i=0; i<nrow; i++)
	{
		free(old_clustIndex[i]);
		free(new_clustIndex[i]);
	}
	free(old_clustIndex);
	free(new_clustIndex);
	free(old_clustSize);	
	free(new_clustSize);
	free(old_cumu_seqNum);
	free(new_cumu_seqNum);

}

// calculate similarity index of two clusters: (size of intersection)/(size of union)
double calc_IOU(int clustIndex1[], int clustSize1, int clustIndex2[], int clustSize2)
{
	int intersectSize, unionSize, isFound,i,j;

	//calculate unionSize
	unionSize=clustSize1;
	for(i=0; i<clustSize2; i++)
	{
		isFound=0;
		for(j=0; j<clustSize1; j++)
		{
			if(clustIndex1[j]==clustIndex2[i])
				isFound=1;
		}
		if(isFound==0) unionSize++;
	}
	
	//calculate intersection
	intersectSize=0;
	for(i=0; i<clustSize2; i++)
	{
		isFound=0;
		for(j=0; j<clustSize1; j++)
		{
			if(clustIndex1[j]==clustIndex2[i])
				isFound=1;
		}
		if(isFound==1) intersectSize++;
	}
	return((intersectSize+0.0)/(unionSize+0.0));
}
void combineChars(char first[], char second[], char combined[])
{
	int i=0, j=0;
	while(first[i]!='\0')
	{
		combined[j]=first[i];
		i++;
		j++;
	}

	i=0;
	while(second[i]!='\0')
	{
		combined[j]=second[i];
		i++;
		j++;
	}
	combined[j]='\0';
}

void assignChars(char oldChars[], char newChars[])
{
	int i=0;
	while(newChars[i]!='\0')
	{
		oldChars[i]=newChars[i];
		i++;
	}
	oldChars[i]='\0';
}


