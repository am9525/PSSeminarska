/*
 * DisplayImage.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: aljaz
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "opencv2/opencv.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/core/core.hpp>
#include "opencv2/imgcodecs.hpp"
#include <opencv2/highgui/highgui.hpp>


using namespace std;
using namespace cv;

//GLOBAL VARS
int numThetas = 180;
int X = 40;
int Y = 40;
double PI = 3.14159265;
double diag = 0.0;
int transformSizeX = 0;
int transformSizeY = 0;
int minNumOfPixelsOnLine = 4;
double outOfBounds = 1000000;

//PROTOTYPES

//arithmetic-logic functions
int** houghField(int** img, int x, int y);
double** getLines(double** htr, double diag);
double** mergeLines(double** lines);
double** calculateNormalizedDistances(double** vertices);
void mergeTwoLines(double* a, double* b);
double** removeOutOfFieldLines(double** lines);
double calculateAngle(double theta);
double calculateX0Intersection(double theta, double r);
double* calculateLine(double theta, double r);
double* findIntersect(double* line0, double* line1);
void sortLines(double** unsorted);
double** blur(double **arr, int x, int y);
int** fixEdges(int** img, int x, int y);
void rotateLines(double** lines);
double** removeUndesiredLines(double** lines);
double** generateIntersectTable(double** lines);

//I/O functions
int** get2dFrom1d(uchar* a, int x, int y);
void outLines(double** lines);
void outArr(int** field, int x, int y);
void outImgI(int** field, int x, int y);
void outImgD(double** field, int x, int y);
void outIntersectTable(double** table);
void outDistances(double** distances);

//shifters&helpers
double** int2double(int** arr, int x, int y);
double** polar2cartesian(double**lines);
int** generateImg();
int** makeIntArr(int x, int y);
/*
int main() {
	diag = sqrt(X * X + Y * Y);
	transformSizeX = 2 + ((int) (2 * diag));
	transformSizeY = numThetas + 2;
	int** img = generateImg();
	outImgI(img, X,Y);
	int **hTform = houghField(img, X, Y);
	hTform = fixEdges(hTform, transformSizeX, transformSizeY);
	double** hTransform = int2double(hTform, transformSizeX, transformSizeY);
	hTransform = blur(hTransform, transformSizeX, transformSizeY);
	double** linesPotentialPolar = getLines(hTransform, diag);
	double**foundLinesPolar = mergeLines(linesPotentialPolar);
	rotateLines(foundLinesPolar);
	foundLinesPolar = mergeLines(linesPotentialPolar);
	foundLinesPolar = removeUndesiredLines(foundLinesPolar);
	sortLines(foundLinesPolar);
	double** foundLinesCartesian = polar2cartesian(foundLinesPolar);
	outLines(foundLinesCartesian);
	double** intersectTable = generateIntersectTable(foundLinesCartesian);
	outIntersectTable(intersectTable);

	double** distances = calculateNormalizedDistances(intersectTable);
	outDistances(distances);
}*/

int** generateImg() {
	int **img = makeIntArr(X, Y);
	for (int i = 5; i <= 20; i++) {
		img[5][i] = 1;
		img[20][i] = 1;
		img[i][5] = 1;
		img[i][20] = 1;
	}
	//outArr(img, X,Y);
	return img;
}

int** houghField(int** img, int x, int y) {

	int** htr = makeIntArr(transformSizeX, transformSizeY);

	int i = 0;
	for (; i < x; i++) {
		int j = 0;
		for (; j < y; j++) {
			if (img[i][j] == 0)
				continue;
			int theta = 0;
			for (; theta < numThetas; theta++) {
				double angle = PI * (((double) theta) / ((double) numThetas));
				double r = (i + 1) * cos(angle) + (j + 1) * sin(angle);
				r += diag;
				htr[1 + (int) r][1 + theta]++;
			}
		}
	}
	/*outImgI(htr, transformSizeX, transformSizeY);
	 fflush (stdout);*/

	return htr;
}

int** fixEdges(int** img, int x, int y) {
	int i = 0;
	for (; i < x; i++) {
		img[i][0] = img[i][1];
		img[i][y - 1] = img[i][y - 2];
	}
	i = 0;
	for (; i < y; i++) {
		img[0][i] = img[1][i];
		img[x - 1][i] = img[x - 2][i];
	}
	return img;
}

double** calculateNormalizedDistances(double** v) {
	int n = (int) v[0][0];
	int numOfEdges = n * (n - 1) / 2;
	double* distances = (double*) malloc(sizeof(double) * numOfEdges);
	int i = 1, edgeIter = 0;
	for (; i <= n; i++) {
		int j = 1;
		for (; j < i; j++) {
			distances[edgeIter++] = sqrt(
					(v[i][0] - v[j][0]) * (v[i][0] - v[j][0])
							+ (v[i][1] - v[j][1]) * (v[i][1] - v[j][1]));
		}
	}

	double**normalizedDistances = (double**) malloc(
			sizeof(double*) * (numOfEdges + 1));
	normalizedDistances[0] = (double*) malloc(sizeof(double));
	normalizedDistances[0][0] = numOfEdges + 1;
	i = 1;
	for (; i <= numOfEdges; i++)
		normalizedDistances[i] = (double*) malloc(sizeof(double) * numOfEdges);
	i = 1;
	for (; i <= numOfEdges; i++) {
		double maxVal = distances[i - 1];
		int j = 0;
		for (; j < numOfEdges; j++) {
			normalizedDistances[i][j] = distances[j] / maxVal;
		}
		/*j=1;
		 for(;j<=numOfEdges;j++){
		 int k=1;
		 for(;k<j;k++){
		 if(normalizedDistances[i][k]>normalizedDistances[i][j]){
		 double tmp=normalizedDistances[i][j];
		 normalizedDistances[i][j]=normalizedDistances[i][k];
		 normalizedDistances[i][k]=tmp;
		 }
		 }
		 }*/
	}
	return normalizedDistances;
}

double** generateIntersectTable(double** lines) {
	double** iTable = (double**) malloc(sizeof(double*));
	iTable[0] = (double*) malloc(sizeof(double));
	iTable[0][0] = 0;
	int i = 1;
	for (; i < lines[0][0]; i++) {
		int j = 1;
		for (; j < i; j++) {
			int iTabIter = iTable[0][0] + 1;
			iTable = (double**) realloc(iTable,
					sizeof(double*) * (iTabIter + 1));
			iTable[iTabIter] = (double*) malloc(sizeof(double) * 2);
			double* intersectionPoint = findIntersect(lines[i], lines[j]);
			iTable[iTabIter][0] = intersectionPoint[0];
			iTable[iTabIter][1] = intersectionPoint[1];
			if (iTable[iTabIter][0] > X + 10 || iTable[iTabIter][0] < -10
					|| iTable[iTabIter][1] > Y + 10
					|| iTable[iTabIter][1] < -10) {
				iTable[0][0]--;
			}
			iTable[0][0]++;
		}
	}
	return iTable;
}

double* findIntersect(double* line0, double* line1) {
	double* vals = (double*) malloc(sizeof(double) * 2);

	if (line0[0] != outOfBounds && line1[0] != outOfBounds) {
		double xintersect = (line1[1] - line0[1]) / (line0[0] - line1[0]);
		double yintersect = line1[1] + xintersect * line1[0];
		vals[0] = xintersect;
		vals[1] = yintersect;
	} else if (line1[0] != outOfBounds) {
		vals[0] = line0[1];
		vals[1] = line1[1] + line0[1] * line1[0];
	} else if (line0[0] != outOfBounds) {
		vals[0] = line1[1];
		vals[1] = line0[1] + line1[1] * line0[0];
	} else {
		vals[0] = outOfBounds;
		vals[1] = outOfBounds;
	}
	return vals;
}

double** getLines(double** htr, double diag) {
	//lines are visible in the houghTransform field because they are the maximal value there
	int k = 0;
	double** lines = (double**) malloc(sizeof(double*));
	lines[0] = (double*) malloc(sizeof(double));
	lines[0][0] = 0;
	//outImgD(blurred, transformSizeX,transformSizeY);
	int i = 1;
	for (; i < (1 + (int) (2 * diag)); i++) {
		int j = 1;
		for (; j < 1 + numThetas; j++) {
			if (htr[i][j] < (double) minNumOfPixelsOnLine)
				continue;
			if (htr[i][j] > htr[i - 1][j - 1] && htr[i][j] > htr[i - 1][j]
					&& htr[i][j] > htr[i - 1][j + 1]
					&& htr[i][j] > htr[i][j - 1] && htr[i][j] > htr[i][j + 1]
					&& htr[i][j] > htr[i + 1][j - 1]
					&& htr[i][j] > htr[i + 1][j]
					&& htr[i][j] > htr[i + 1][j + 1]) {

				lines = (double**) realloc(lines,
						(int) (lines[0][0] + 2) * sizeof(double*));
				int theta = j - 1;
				int r = i - diag - 1;
				lines[(int) lines[0][0] + 1] = (double*) malloc(
						sizeof(double) * 3);
				lines[(int) lines[0][0] + 1][0] = (double) theta
						/ (double) numThetas;
				lines[(int) lines[0][0] + 1][1] = (double) r;
				lines[(int) lines[0][0] + 1][2] = htr[i][j];
				lines[0][0]++;
			}
		}
	}
	lines[0][0]++;
	return lines;
}

double* calculateLine(double theta, double r) { //y= ax+b => ret[0]=a;ret[1]=b
	double* vals = (double*) malloc(sizeof(double) * 2);
	vals[0] = calculateAngle(theta);
	vals[1] = calculateX0Intersection(theta, r);
	return vals;
}

double calculateAngle(double theta) {
	double angle = tan(PI * theta);
	return angle;
}

double calculateX0Intersection(double theta, double r) {
	double cosAng = cos(PI * theta);
	double intersect = r / cosAng;
	return intersect;
}

double** int2double(int** arr, int x, int y) {
	double ** ret = (double**) malloc(sizeof(double*) * (x));
	int i = 1;
	ret[0] = (double*) malloc(sizeof(double) * (y));
	for (; i < x - 1; i++) {
		ret[i] = (double*) malloc(sizeof(double) * (y));
		int j = 1;
		for (; j < y - 1; j++) {
			ret[i][j] = (double) arr[i][j];
		}
	}
	ret[x - 1] = (double*) malloc(sizeof(double) * (y));

	return ret;
}

double** polar2cartesian(double**lines) {
	int i = 1;
	for (; i < lines[0][0]; i++) {
		if (abs(lines[i][0] - 0.5) < 0.1 || abs(lines[i][0] + 0.5) < 0.1) {
			lines[i][0] = outOfBounds; //vertical line!!!! this is only a placeholder for some other
		} else {
			double* cartesianLine = calculateLine(lines[i][0], lines[i][1]);
			lines[i][0] = cartesianLine[0];
			lines[i][1] = cartesianLine[1];
		}

	}
	return lines;
}

double** blur(double **arr, int x, int y) {
	double box[3][3] = { { 1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0 }, { 1.0 / 8.0, 1.0
			/ 4.0, 1.0 / 8.0 }, { 1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0 } }; //SOME KIND OF GAUSSIAN?
	double ** ret = (double**) malloc(sizeof(double*) * (x));
	int i = 1;
	ret[0] = (double*) malloc(sizeof(double) * (y));
	for (; i < x - 1; i++) {
		ret[i] = (double*) malloc(sizeof(double) * (y));
		int j = 1;
		for (; j < y - 1; j++) {
			ret[i][j] = arr[i - 1][j - 1] * box[0][0];
			ret[i][j] += arr[i - 1][j] * box[0][1];
			ret[i][j] += arr[i - 1][j + 1] * box[0][2];
			ret[i][j] += arr[i][j - 1] * box[1][0];
			ret[i][j] += arr[i][j] * box[1][1];
			ret[i][j] += arr[i][j + 1] * box[1][2];
			ret[i][j] += arr[i + 1][j - 1] * box[2][0];
			ret[i][j] += arr[i + 1][j] * box[2][1];
			ret[i][j] += arr[i + 1][j + 1] * box[2][2];
		}
	}
	ret[x - 1] = (double*) malloc(sizeof(double) * (y));
	return ret;
}

void sortLines(double** unsorted) {
	double numOfLines = unsorted[0][0];
	int i = 1;
	for (; i < numOfLines; i++) {
		int j = i + 1;
		for (; j < numOfLines; j++) {
			if (unsorted[i][2] < unsorted[j][2]) {
				double *tmp = unsorted[i];
				unsorted[i] = unsorted[j];
				unsorted[j] = tmp;
			}
		}
	}
}

void rotateLines(double** lines) {
	double numOfLines = lines[0][0];
	int i = 1;
	for (; i < numOfLines; i++) {
		if (lines[i][0] > 0.5 && lines[i][1] < 0) {
			lines[i][0] -= 1.0;
			lines[i][1] *= -1.0;
		}
	}
}

double** mergeLines(double** lines) {
	double numOfLines = lines[0][0];
	int i = 1;
	for (; i < numOfLines; i++) {
		if (lines[i][0] == 0 && lines[i][1] == -1) {
			continue;
		}
		int j = i + 1;
		for (; j < numOfLines; j++) {
			if (lines[j][0] == 0 && lines[j][1] == -1)
				continue;
			double absA = abs(lines[i][0] - lines[j][0]);
			double absB = abs(lines[i][1] - lines[j][1]);

			if (absA < 0.08 && absB <= 2) {
				//are simmilar!!! Merge them in to one line!
				mergeTwoLines(lines[i], lines[j]);
				i--;
				break;
			}
		}
	}
	double** cleanedLines = removeOutOfFieldLines(lines);
	//outLines(cleanedLines);
	return cleanedLines;
}

void mergeTwoLines(double* a, double* b) {
	double sumIntensity = a[2] + b[2];
	double aFactor = a[2] / sumIntensity;
	double bFactor = b[2] / sumIntensity;
	a[0] = a[0] * aFactor + b[0] * bFactor;
	a[1] = a[1] * aFactor + b[1] * bFactor;
	a[2] = sqrt(a[2] * a[2] + b[2] * b[2]);
	b[0] = 0;
	b[1] = -1;
	b[2] = 0;
}

double** removeOutOfFieldLines(double** lines) {
	double numOfLines = lines[0][0];
	int i = 1, newNumOfLines = 0;
	for (; i < numOfLines; i++) {
		if (lines[i][0] != 0 && lines[i][0] != -1 && lines[i][2] != 0)
			newNumOfLines++;
	}
	newNumOfLines++;
	double** linesNew = (double**) malloc(newNumOfLines * sizeof(double*));
	linesNew[0] = (double*) malloc(sizeof(double));
	linesNew[0][0] = newNumOfLines;
	i = 1;
	int iNew = 1;
	for (; i < numOfLines; i++) {
		if (lines[i][0] != 0 && lines[i][0] != -1 && lines[i][2] != 0)
			linesNew[iNew++] = lines[i];
	}
	return linesNew;
}

double** removeUndesiredLines(double** lines) {
	double numOfLines = lines[0][0], sum = 0.0;
	int i = 1;
	for (; i < numOfLines; i++) {
		sum += lines[i][2];
	}
	double avgIntensity = sum / (numOfLines - 1.0);
	i = 1;
	for (; i < numOfLines; i++) {
		if (lines[i][2] < avgIntensity) {
			lines[i][0] = 0;
			lines[i][1] = -1;
			lines[i][2] = 0;
		}
	}
	return removeOutOfFieldLines(lines);
}

int** makeIntArr(int x, int y) {
	int **ret = (int**) malloc(x * sizeof(int*));
	int i = 0;
	for (; i < x; i++) {
		ret[i] = (int*) malloc(y * sizeof(int));
		int j = 0;
		for (; j < y; j++) {
			ret[i][j] = 0;
		}
	}
	return ret;
}

void outLines(double** lines) {
	//printf("FORMAT OF OUTPUT:\n [i]: A:[a], B:[b], V:[v]; \n	where i is the number of the line coordinates,\n	Form: y=ax+b \n	and v denotes the number of pixels on that line\n\n");
	int i = 1;
	printf("Number of lines:%d\n", (int) lines[0][0] - 1);
	for (; i < lines[0][0]; i++) {
		printf("%d: A:%f, B:%f, V:%f\n", i, lines[i][0], lines[i][1],
				lines[i][2]);
	}
	printf("\n");
}

void outArr(int** field, int x, int y) {
	int i = 0;
	for (; i < x; i++) {
		int j = 0;
		for (; j < y; j++) {
			printf("%d ", field[i][j]);
		}
		printf("\n");
	}
}

void outImgI(int** field, int x, int y) {
	double max = 0;
	int i = 0;
	for (; i < x; i++) {
		int j = 0;
		for (; j < y; j++) {
			if (field[i][j] > max)
				max = field[i][j];
		}
	}

	printf("int[][] img={");
	i = 0;
	for (; i < x; i++) {
		int j = 0;
		printf("{");
		for (; j < y; j++) {
			printf("%d, ", (int) (255 * ((double) field[i][j] / max)));
		}
		printf("},\n");
	}
	printf("};\n");
}

void outImgD(double** field, int x, int y) {
	double max = 0;
	int i = 0;
	for (; i < x; i++) {
		int j = 0;
		for (; j < y; j++) {
			if (field[i][j] > max)
				max = field[i][j];
		}
	}
	i = 0;
	printf("int[][] img={");
	for (; i < x; i++) {
		int j = 0;
		printf("{");
		for (; j < y; j++) {
			printf("%d, ", (int) (255.0 * field[i][j] / max));
		}
		printf("},\n");
	}
	printf("};\n");
}

void outIntersectTable(double** table) {
	printf("Calculated vertices of figure:\n");
	double numIntersects = table[0][0] + 1;
	int i = 1;
	for (; i < numIntersects; i++) {
		printf("[%d, %d], \n", (int) table[i][0], (int) table[i][1]);
	}
	printf("\n");

}

int** get2dFrom1d(uchar* a, int rows, int cols){
	int** table=(int**)malloc(sizeof(int*)*rows);
	int i=0;
	for(;i<rows;i++){
		table[i]=(int*)malloc(sizeof(int)*cols);
		int j=0;
		for(;j<cols;j++){
			table[i][j]=(int)a[(j*3)+(i*cols*3)];
		}
	}
	return table;
}

void outDistances(double** distances) {
	int numOfDists = distances[0][0];
	int i = 1;
	printf("Calculated distances of figure:\n");
	for (; i < numOfDists; i++) {
		int j = 0;
		for (; j < numOfDists - 1; j++) {
			printf("%f,	", distances[i][j]);
		}
		printf("\n");
	}
}




void RgbToHsv(Mat * OrigImage, Mat * HSVImage)
{
	*HSVImage = OrigImage->clone();
    unsigned char s,v;
    unsigned char rgbMin, rgbMax;
    for(int i = 0; i < OrigImage->rows * OrigImage->cols * 3; i+=3){

    	//finds the minimum value between r,g,b
    	// i = b, i+1 = g, i+2 = r
        rgbMin = OrigImage->data[i+2] < OrigImage->data[i+1]  ? (OrigImage->data[i+2]< OrigImage->data[i] ? OrigImage->data[i+2] : OrigImage->data[i]) : (OrigImage->data[i+1] < OrigImage->data[i] ? OrigImage->data[i+1] : OrigImage->data[i]);
        rgbMax = OrigImage->data[i+2] > OrigImage->data[i+1] ? (OrigImage->data[i+2] > OrigImage->data[i] ? OrigImage->data[i+2] : OrigImage->data[i]) : (OrigImage->data[i+1] > OrigImage->data[i] ? OrigImage->data[i+1] : OrigImage->data[i]);

        v = rgbMax;
        //if color is black
        if (v == 0){
            HSVImage->data[i] = 0;//h
            HSVImage->data[i+1] = 0;//s
            HSVImage->data[i+2] = v;//v
            continue;
        }
        //saturation
        //if color is white
        s = 255* long(rgbMax - rgbMin) / v;
        if (s == 0){
            HSVImage->data[i] = 0;//h
            HSVImage->data[i+1] = s;//s
            HSVImage->data[i+2] = v;//v
            continue;
        }

        if (rgbMax == OrigImage->data[i+2])
        	HSVImage->data[i]  = 0 + 43 * (OrigImage->data[i+1] - OrigImage->data[i]) / (rgbMax - rgbMin);
		else if (rgbMax == OrigImage->data[i+1])
			HSVImage->data[i] = (uchar)(85 + 43 * (OrigImage->data[i] - OrigImage->data[i+2]) / (rgbMax - rgbMin));
		else
			HSVImage->data[i] = (uchar)(171 + 43 * (OrigImage->data[i+2] - OrigImage->data[i+1]) / (rgbMax - rgbMin));
        HSVImage->data[i+1] = s;
        HSVImage->data[i+2] = v;
    }

}
//converts char array to int array
void imgToIntArray(Mat * image, int * manArray){
	for(int i = 0; i < image->rows*image->cols*3; i++){
		//cout << "i: "<< i << " data: " << (int)image->data[i] << endl;
		manArray[i] = (int)image->data[i];
	}
}

int * diliate(Mat * image, int k){
	int inRow = image->cols;
	int inCol = image->rows;
	int manhattanArray[image->rows*image->cols*3];
	imgToIntArray(image, manhattanArray);

	for(int i = 0; i < inRow*inCol*3; i+=3){
		//ce smo najdli bel pixel
		if(image->data[i] == 255 && image->data[i+1] == 255 && image->data[i+2] == 255){
			//first pass and pixel was on, it gets a zero
			manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] =  0;
		}
		else{
            // pixel was off
            // It is at most the sum of the lengths of the array
            // away from a pixel that is on
			manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] = inRow+inCol;

			//ce i ni enak prvemu stolpcu,
			//pogledamo ali je manj oddaljen od najvecje mozne oddaljenosti
			//pogledamo gor
			if(i >= (inRow*3)){
				//cout << "hello" << endl;
				manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] = min((int)manhattanArray[i], (int)manhattanArray[i-inRow*3]+1);
			}
			// nismo v prvi vrstici
			if(i%(inRow*3)){
				manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2]  = min((int)manhattanArray[i], (int)manhattanArray[i-3]+1);
			}

			if(manhattanArray[i] <= k)
				image->data[i] = image->data[i+1] = image->data[i+2] = 255;

			else
				image->data[i] = image->data[i+1] = image->data[i+2] = 0;
		}
	}
	for(int i = inRow*inCol*3-3; i >= 0; i-=3){
		//cout << "i: "<< i << endl;
		if(i+inRow*3 < inRow*inCol*3){
			//cout << "hello"<< endl;
			manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] = min((int)manhattanArray[i], (int)manhattanArray[i+inRow*3]+1);
		}
		if(i%(inRow*3) != (inRow*3-3)){
			//cout << (int)image->data[i] << " compare "<< (int)image->data[i+3]+1<< endl;
			manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] = min((int)manhattanArray[i], (int)manhattanArray[i+3]+1);
		}
		if(manhattanArray[i] <= k)
			image->data[i] = image->data[i+1] = image->data[i+2] = 255;

		else
			image->data[i] = image->data[i+1] = image->data[i+2] = 0;
	}

    return manhattanArray;
}
int * erode(Mat * image, int k){
	int inRow = image->cols;
	int inCol = image->rows;
	int manhattanArray[image->rows*image->cols*3];
	imgToIntArray(image, manhattanArray);

	for(int i = 0; i < inRow*inCol*3; i+=3){
		//ce smo najdli bel pixel
		//cout << "i: "<< (int)i << "data: "<<(int)image->data[i];
		if(image->data[i] == 0 && image->data[i+1] == 0 && image->data[i+2] == 0){
			//first pass and pixel was on, it gets a zero
			manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] =  0;
		}
		else{
            // pixel was off
            // It is at most the sum of the lengths of the array
            // away from a pixel that is on
			manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] = inRow+inCol;

			//ce i ni enak prvemu stolpcu,
			//pogledamo ali je manj oddaljen od najvecje mozne oddaljenosti
			//pogledamo gor
			if(i >= (inRow*3)){
				//cout << "hello" << endl;
				manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] = min((int)manhattanArray[i], (int)manhattanArray[i-inRow*3]+1);
			}
			// nismo v prvi vrstici
			if(i%(inRow*3)){
				manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2]  = min((int)manhattanArray[i], (int)manhattanArray[i-3]+1);
			}
			if(manhattanArray[i] <= k)
				image->data[i] = image->data[i+1] = image->data[i+2] = 0;

			else
				image->data[i] = image->data[i+1] = image->data[i+2] = 255;
		}
	}
	for(int i = inRow*inCol*3-3; i >= 0; i-=3){
		//cout << "i: "<< i << endl;
		if(i+inRow*3 < inRow*inCol*3){
			//cout << "hello"<< endl;
			manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] = min((int)manhattanArray[i], (int)manhattanArray[i+inRow*3]+1);
		}
		if(i%(inRow*3) != (inRow*3-3)){
			//cout << (int)image->data[i] << " compare "<< (int)image->data[i+3]+1<< endl;
			manhattanArray[i] = manhattanArray[i+1] = manhattanArray[i+2] = min((int)manhattanArray[i], (int)manhattanArray[i+3]+1);
		}
		if(manhattanArray[i] <= k)
			image->data[i] = image->data[i+1] = image->data[i+2] = 0;
		else
			image->data[i] = image->data[i+1] = image->data[i+2] = 255;
	}
    return manhattanArray;
}

void imgInRange(Mat * image, int LowH, int LowS, int LowV, int HighH, int HighS, int HighV){
    //cout << "helo"<< endl;
	for(int i = 0; i < image->rows*image->cols*3; i+=3){
		//cout << "i: " << i << endl;
		cout << (int)image->data[i] << " " <<  (int)image->data[i+1] << " "<< (int)image->data[i+2] << endl;
		if(image->data[i] >= LowH && image->data[i] <= HighH && image->data[i+1] >= LowS && image->data[i+1] <= HighS && image->data[i+2] >= LowV && image->data[i+2] <= HighV){
			image->data[i] = 255;
			image->data[i+1] = 255;
			image->data[i+2] = 255;
		}else{
			image->data[i] = 0;
			image->data[i+1] = 0;
			image->data[i+2] = 0;
		}
	}
}

void calMoments(Mat * image, Mat * origImage){
	int zerothMoment = 0;
	for(int i = 0; i < image->rows*image->cols*3; i+=3)
		if(image->data[i] == 255)
			zerothMoment += 1;
	//calulate 10 and 01 moments
	int moment10 = 0;
	int moment01 = 0;
	for(int i = 0; i < image->rows*image->cols*3; i+=3){
		//cout << i << endl;
		if(image->data[i] == 255){
			//cout << "x :" << (i/3)%image->cols << " y: " << (i/3)/ image->cols << endl;
			moment01+= (i/3)%image->cols;
			moment10 += (i/3)/image->cols;
		}
	}
	if(zerothMoment > 1000){
		int xCor = moment01 / zerothMoment;
		int yCor= moment10 / zerothMoment;
		circle(*origImage, Point(xCor, yCor), 32.0, Scalar( 0, 255, 0 ), 1, 8 );
	}



}
int main( int argc, char** argv ){

	Mat image, HSVimage, imgThresholded;
	image = imread( "pixel21.png", 1 );

	int iLowH = 0;
	int iHighH = 0;
	int iLowS = 255;
	int iHighS = 255;

	int iLowV = 255;
	int iHighV = 255;


	if((!image.data)){
	  printf( "No image data \n" );
	  return -1;
	}
	HSVimage = image.clone();

	//cvtColor(imgOriginal, imgHSV, COLOR_BGR2HSV); //Convert the captured frame from BGR to HSV
	RgbToHsv(&image, &HSVimage);
	imgThresholded = HSVimage.clone();

	//inRange(imgHSV, Scalar(iLowH, iLowS, iLowV), Scalar(iHighH, iHighS, iHighV), imgThresholded); //Threshold the image
	imgInRange(&imgThresholded, iLowH, iLowS, iLowV, iHighH, iHighS, iHighV);
	//morphological opening (remove small objects from the foreground)
	//dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
	//erode(&imgThresholded, 10);
	//diliate(&imgThresholded, 10);

	//morphological closing (fill small holes in the foreground)
	//diliate(&imgThresholded, 2);
	//erode(&imgThresholded, 2);

	calMoments(&imgThresholded, &image);
	////////////////////////////////
	X = imgThresholded.rows;
	Y = imgThresholded.cols;
	diag = sqrt(X * X + Y * Y);
	transformSizeX = 2 + ((int) (2 * diag));
	transformSizeY = numThetas + 2;
	int** img = get2dFrom1d(image.data, X, Y);
	outImgI(img, X,Y);
	int **hTform = houghField(img, X, Y);
	hTform = fixEdges(hTform, transformSizeX, transformSizeY);
	double** hTransform = int2double(hTform, transformSizeX, transformSizeY);
	hTransform = blur(hTransform, transformSizeX, transformSizeY);
	double** linesPotentialPolar = getLines(hTransform, diag);
	double**foundLinesPolar = mergeLines(linesPotentialPolar);
	rotateLines(foundLinesPolar);
	foundLinesPolar = mergeLines(linesPotentialPolar);
	foundLinesPolar = removeUndesiredLines(foundLinesPolar);
	sortLines(foundLinesPolar);
	double** foundLinesCartesian = polar2cartesian(foundLinesPolar);
	outLines(foundLinesCartesian);
	double** intersectTable = generateIntersectTable(foundLinesCartesian);
	outIntersectTable(intersectTable);

	double** distances = calculateNormalizedDistances(intersectTable);
	outDistances(distances);
	//////////////////////////////
	//imshow( "Display Image5", imgThresholded);
	imshow("Thresholded Image", imgThresholded); //show the thresholded image
	imshow("Original", image); //show the original image


/*
	//declare Mat image files
	Mat imgOriginal, imgHSV, imgThresholded ;
	VideoCapture cap(0); //capture the video from web cam

	if ( !cap.isOpened() )  // if not success, exit program
	{
		 cout << "Cannot open the web cam" << endl;
		 return -1;
	}

	namedWindow("Control", CV_WINDOW_AUTOSIZE); //create a window called "Control"

	int iLowH = 141;
	int iHighH = 179;

	int iLowS = 14;
	int iHighS = 255;

	int iLowV = 5;
	int iHighV = 255;

	//Create trackbars in "Control" window
	cvCreateTrackbar("LowH", "Control", &iLowH, 255); //Hue (0 - 255)
	cvCreateTrackbar("HighH", "Control", &iHighH, 255);

	cvCreateTrackbar("LowS", "Control", &iLowS, 255); //Saturation (0 - 255)
	cvCreateTrackbar("HighS", "Control", &iHighS, 255);

	cvCreateTrackbar("LowV", "Control", &iLowV, 255); //Value (0 - 255)
	cvCreateTrackbar("HighV", "Control", &iHighV, 255);

	while (true)
	{

		bool bSuccess = cap.read(imgOriginal); // read a new frame from video

		// If you do not care about backward compatibility
		// You can use the following instead for OpenCV 3
		// double fps = video.get(CAP_PROP_FPS);

		if (!bSuccess) //if not success, break loop
		{
		   cout << "Cannot read a frame from video stream" << endl;
		   break;
		}



		//cvtColor(imgOriginal, imgHSV, COLOR_BGR2HSV); //Convert the captured frame from BGR to HSV
		RgbToHsv(&imgOriginal, &imgHSV);

		imgThresholded = imgHSV.clone();

		//inRange(imgHSV, Scalar(iLowH, iLowS, iLowV), Scalar(iHighH, iHighS, iHighV), imgThresholded); //Threshold the image
		imgInRange(&imgThresholded, iLowH, iLowS, iLowV, iHighH, iHighS, iHighV);
		//morphological opening (remove small objects from the foreground)
		erode(&imgThresholded, 3);
		diliate(&imgThresholded, 3);

		//morphological closing (fill small holes in the foreground)
		diliate(&imgThresholded, 3);
		erode(&imgThresholded, 3);


		calMoments(&imgThresholded, &imgOriginal);

		imshow("Thresholded Image", imgThresholded); //show the thresholded image
		imshow("Original", imgOriginal); //show the original image

		if (waitKey(30) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
		{
		  cout << "esc key is pressed by user" << endl;
		  break;
		}
	}*/
	waitKey(0);
	return 0;
}

