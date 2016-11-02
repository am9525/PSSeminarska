#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int numThetas=40;
int X= 20;
int Y=20;
double PI = 3.14159265;



double** blur(int **arr, int x, int y);
int** houghField(int** img, int x, int y);
int** makeIntArr(int x, int y);
void outImg(int** field, int x, int y);
void outImgD(double** field, int x, int y);

int main(){
	int **img = makeIntArr(X, Y);
	int i=0;
	for(;i<X;i++)img[X-i-1][i]=1;
	outImg(img, X,Y);
	houghField(img, X, Y);
}

int** houghField(int** img, int x, int y){
	double diag = sqrt(x*x + y*y);
	// adding 2 because we want an edge arround it, for the blurring.
	int** htr = makeIntArr(2+((int)(2*diag)), numThetas+2); 
	int i=0;
	for(;i<x;i++){
		int j=0;
		for(;j<y;j++){
			if(img[i][j]==0)
				continue;
			int theta=0;
			for(;theta<numThetas;theta++){
				double angle = PI*(((double)theta)/((double)numThetas));
				double r = (i+1)*cos(angle)+(j+1)*sin(angle);
				r+=diag;
				htr[1+(int)r][1+theta]++;
			}
		}
	}
	outImg(htr, 2+((int)(2*diag)), 2+numThetas);
	
	
	double** blurred = blur(htr,2+((int)(2*diag)), 2+numThetas);
	int maxR=0,maxTheta=0;
	i=0;
	for(;i<((int)(2*diag));i++){
		int j=0;
		for(;j<numThetas;j++){
			if(blurred[i][j]>blurred[maxR][maxTheta]){
				maxR=i;
				maxTheta=j;
			}
		}
	}
	maxR-=diag;
	printf("MAXES:\n R:%d, Th:%f\n", maxR, ((double)maxTheta/(double)numThetas));
	double a=tan(PI*((double)maxTheta/(double)numThetas));
	double b=(double)maxR/cos(PI*((double)maxTheta/(double)numThetas));
	printf("\nA:%f,	B:%f\n",a,b);
	
	return htr;
}

int** makeIntArr(int x, int y){
	int **ret = (int**)malloc(x*sizeof(int*));
	int i=0;
	for(;i<x;i++){
		ret[i]=(int*)malloc(y*sizeof(int));
		int j=0;
		for(;j<y;j++){
			ret[i][j] = 0;
		}
	}
	return ret;
}

double** blur(int **arr, int x, int y){
	double box[3][3]={
	{1.0/18.0, 1.0/9.0, 1.0/18.0},
	{1.0/9.0, 1.0/3.0, 1.0/9.0},
	{1.0/18.0, 1.0/9.0, 1.0/18.0}};
	double ** ret= malloc(sizeof(double*)*(x-2));
	int i=1;
	for(;i<x-1;i++){
		ret[i-1]=malloc(sizeof(double)*(y-2));
		int j=1;
		for(;j<y-1;j++){
			ret[i-1][j-1]= ((double)arr[i-1][j-1])*box[0][0];
			ret[i-1][j-1]+=((double)arr[i-1][j])*box[0][1];
			ret[i-1][j-1]+=((double)arr[i-1][j+1])*box[0][2];
			ret[i-1][j-1]+=((double)arr[i][j-1])*box[1][0];
			ret[i-1][j-1]+=((double)arr[i][j])*box[1][1];
			ret[i-1][j-1]+=((double)arr[i][j+1])*box[1][2];
			ret[i-1][j-1]+=((double)arr[i+1][j-1])*box[2][0];
			ret[i-1][j-1]+=((double)arr[i+1][j])*box[2][1];
			ret[i-1][j-1]+=((double)arr[i+1][j+1])*box[2][2];
		}
	}
	
	return ret;
}
void outImg(int** field, int x, int y){
	printf("Drawing: \n");
	int i=0;
	for(;i<x;i++){
		int j=0;
		for(;j<y;j++){
			char out=' ';
			/*if(field[i][j]>0 && field[i][j]<3)out = '.';
			if(field[i][j]>=3 && field[i][j]<20 )out = '+';
			if(field[i][j]>=20 && field[i][j]<80 )out = 'X';
			if(field[i][j]>=80)out = '*';
			printf("%c ",out);*/
			printf("%d ", field[i][j]);
		}
		printf("\n");
	}
}
void outImgD(double** field, int x, int y){
	printf("Drawing: \n");
	int i=0;
	for(;i<x;i++){
		int j=0;
		for(;j<y;j++){
			printf("%f ", field[i][j]);
		}
		printf("\n");
	}
}

