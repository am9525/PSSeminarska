#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdint>
#include <pthread.h>
#include <assert.h>

#include "hough.h"

int** generateImg();
double ** int2double(int** a, int rows, int cols);
FuncArgs createGenericArgs();
//GLOBAL VARS
int numThetas = 240;
int orgRows = 70;
int orgCols = 70;
double diag = 0.0;
int minNumOfPixelsOnLine = 4;
int numThreads=1;
int tfRows;
int tfCols;


int main() {
	int** img = generateImg();

	FuncArgs threadArgs[numThreads];
	pthread_t threads[numThreads];

	diag = sqrt(orgRows * orgRows + orgCols * orgCols);
	tfRows=2 + ((int) (2 * diag));
	tfCols=2 + numThetas;
	double** orgImg=int2double(img, orgRows, orgCols);

	double step=1e-10+(double)orgRows/(double)numThreads;
	int threadCompleted;

	double** htr=makeDoublearray(tfRows,  tfCols);
	int i=0;
	for(;i<numThreads;i++){
		threadArgs[i]=createGenericArgs();
		threadArgs[i].arrIn=orgImg;
		threadArgs[i].arrOut=htr;

		threadArgs[i].threadStartRow=(int)((double)i*step);
		threadArgs[i].threadEndRow=(int)((double)(i+1)*step);


		//calculating hough field
		//houghField(&hough[i]);
		threadCompleted=pthread_create(threads+i, NULL, houghField, &threadArgs[i]);
		assert(!threadCompleted);
	}
	i=0;
	for(;i<numThreads;i++){
		threadCompleted = pthread_join( threads[i], NULL );
		assert(!threadCompleted);
	}

	//changing step size. we have been bound to the
	step=1e-10+(double)tfRows/(double)numThreads;

	FuncArgs linArgs=createGenericArgs();
	linArgs.arrIn=htr;
	linArgs.threadStartRow=0;
	linArgs.threadEndRow=linArgs.tfRows;
	//fixing edges of created transformation
	fixEdges(&linArgs);
	double ** fixedEdges=linArgs.arrOut;
	double** blurHtr=makeDoublearray(tfRows, tfCols);

	i=0;
	for(;i<numThreads;i++){
		threadArgs[i]=createGenericArgs();
		threadArgs[i].arrIn=fixedEdges;
		threadArgs[i].arrOut=blurHtr;

		threadArgs[i].threadStartRow=(int)((double)i*step);
		threadArgs[i].threadEndRow=(int)((double)(i+1)*step);

		//blurring the transformation to dampen the smaller peaks that may have occurred
		//blur(&hough[i]);
		threadCompleted=pthread_create(threads+i, NULL, blur, &threadArgs[i]);
		assert(!threadCompleted);
	}
	i=0;
	for(;i<numThreads;i++){
		threadCompleted = pthread_join( threads[i], NULL );
		assert(!threadCompleted);
	}
	//outImgD(blurHtr, tfRows, tfCols);




	linArgs.arrIn=blurHtr;
	linArgs.arrOut=htr;
	//outImgD(fa.arrIn, fa.tfRows, fa.tfCols);

	//getting desired lines
	getDesiredLines(&linArgs);
	//merging any simmilar lines
	mergeLines(&linArgs);
	//rotating lines that are close to the upper border of theta down (newTheta=theta-1; newR=-r;)
	rotateLines(&linArgs);
	//merging any simmilar lines
	mergeLines(&linArgs);
	linArgs.lines=removeUndesiredLines(&linArgs);
	//outputting gotten lines
	sortLines(&linArgs);
	//outLines(&fa);

	polar2cartesian(&linArgs);
	outLines(&linArgs);


	/*double** intersectTable = generateIntersectTable(&linArgs);
	outIntersectTable(intersectTable);

	double** distances = calculateNormalizedDistances(intersectTable);
	outDistances(distances);*/
}

FuncArgs createGenericArgs(){
	FuncArgs fa;
	fa.tfCols=tfCols;
	fa.tfRows=tfRows;
	fa.diag=diag;
	fa.numThetas=numThetas;
	fa.minNumPixelsOnLine=minNumOfPixelsOnLine;
	fa.orgRows=orgRows;
	fa.orgCols=orgCols;

	return fa;
}

int** generateImg() {
	/*int **img = makeIntArr(X, Y);
	int sX=10,sY=10;
	int edge=20;
	for (int i = 0; i <= edge; i++) {
		/*img[5][i] = 1;
		img[20][i] = 1;
		img[i][5] = 1;
		img[i][20] = 1;* /
		img[sX+i][sY+edge+i]=1;
		img[sX+i][sY+edge-i]=1;
		img[sX+2*edge-i][sY+edge-i]=1;
		img[sX+2*edge-i][sY+edge+i]=1;

	}
	//outArr(img, X,orgCols);
	return img;*/
	int sizeX, sizeY;
	scanf("%d", &sizeX);
	scanf("%d", &sizeY);
	orgRows=sizeX;
	orgCols=sizeY;
	int** img= (int**)malloc(sizeX*sizeof(int*));
	int i =0;
	for (; i < sizeX; i++) {
		img[i]= (int*)malloc(sizeY*sizeof(int));
		int j = 0;
		for (; j < sizeY; j++) {
			scanf("%d", &img[i][j]);
			if(img[i][j]<200)img[i][j]=0;
		}
	}
	return img;
}

double ** int2double(int** a, int rows, int cols) {
	double ** ret = (double**) malloc(sizeof(double*)*rows);
	int i = 0;
	ret[0] = (double*) malloc(sizeof(double) *cols);
	for (; i < rows; i++) {
		ret[i] = (double*) malloc(sizeof(double) *cols);
		int j = 0;
		for (; j < cols; j++) {
			ret[i][j] = (double) (a[i][j]);
		}
	}

	return ret;
}

