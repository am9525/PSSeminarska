/*
 * DisplayImage.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: aljaz
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <omp.h>
#include <math.h>
#include <assert.h>
#include "opencv2/opencv.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/core/core.hpp>
#include "opencv2/imgcodecs.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <time.h>
#include <fstream>
#include "Hough.h"

#define NTHREADS 4
using namespace std;
using namespace cv;
//global variables for calMoments variables
int zerothMoment = 0;
int moment10 = 0;
int moment01 = 0;
//manhattan distance variables
int infinity;
int m,n;
int * g;
int *dt;
omp_lock_t lock;

double ** int2double(int** a, double** b, int rows, int cols);
FuncArgs createGenericArgs();
//GLOBAL VARS
int numThetas = 240;
int orgRows = 70;
int orgCols = 70;
double diag = 0.0;
int minNumOfPixelsOnLine = 4;
const int numThreads = 4;
int tfRows;
int tfCols;

void RgbToHsv(Mat * OrigImage, Mat * HSVImage, int threadID)
{
	int start = (int)( ( (OrigImage->rows * OrigImage->cols * 3) / (double)NTHREADS)*threadID);
	int stop = (int)( ( (OrigImage->rows * OrigImage->cols * 3) / (double)NTHREADS)*(threadID+1));

    unsigned char s,v;
    unsigned char rgbMin, rgbMax;
	for(int i = start; i < stop; i+=3){
		//cout << omp_get_thread_num() << endl;
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
//function used for manhattan distance
int f(int x, int i, int g_i){
	return abs(x-i)+g_i;
}
//function used for manhattan distance
int sep(int i, int u, int g_i, int g_u){
	if(g_u >= (g_i+u-i))
		return infinity;
	if(g_i > (g_u+u-i))
		return -infinity;
	return floor((g_u-g_i+u+i)/2);
}
void printArr(int * arr, int m, int n){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++)
			cout << arr[j+i*n] << " ";
		cout << endl;
	}
	cout << endl;
}
void erode(Mat * image, int k, int threadID){
	//one of the threads is responsible for allocation
	if(threadID == 0){
		//alloacate arrays for maeijster distanc
		g = new int[m*n];
		dt = new int[m*n];
	}

	#pragma omp barrier

	//first phase
	for(int x = threadID; x < m; x+=NTHREADS){
		if(image->data[(x+0*m)*3] == 0)
			g[x+0*m] = 0;
		else
			g[x+0*m] = infinity;
		//1 scan
		for(int y = 1; y < n; y++){
			if(image->data[(x+y*m)*3] == 0)
				g[x+y*m] = 0;
			else
				g[x+y*m] = 1 + g[x+(y-1)*m];
		}

		//Scan 2
		for(int y = n-2; y >= 0; y--){
			if(g[x+(y+1)*m] < g[x+y*m])
				g[x+y*m] = 1 + g[x+(y+1)*m];
		}

	}

	#pragma omp barrier
	//second phase
	int * s = new int[m];
	int * t = new int[m];

	for(int y = threadID; y < n; y+=NTHREADS){
		int q = 0;
		int w  = 0;
		s[0] = 0;
		t[0] = 0;
		//Scan 3
		for(int u = 1; u < m; u++){
			while(q >= 0 && f(t[q], s[q], g[s[q]+y*m]) > f(t[q], u, g[u+y*m]))
				q--;
			if(q < 0){
				q = 0;
				s[0] = u;
			}
			else{
				w = 1 + sep(s[q], u, g[s[q]+y*m], g[u+y*m]);
				if(w < m){
					q++;
					s[q] = u;
					t[q] = w;
				}
			}
		}

		//scan 4
		for(int u = m-1; u >= 0; u--){
			int d = f(u, s[q], g[s[q]+y*m]);
			dt[u+y*m] = d;
			if(d <= k)
				image->data[(u+y*m)*3] = image->data[(u+y*m)*3+1] = image->data[(u+y*m)*3+2] = 0;
			else
				image->data[(u+y*m)*3] = image->data[(u+y*m)*3+1] = image->data[(u+y*m)*3+2] = 255;
			if(u == t[q])
				q--;
		}
	}
	//cout << "hello" << endl;
	#pragma omp barrier
	if(threadID == 0){
		free(dt);
		free(g);
	}
	free(s);
	free(t);
}
void diliate(Mat * image, int k, int threadID){
	//one of the threads is responsible for allocation
	if(threadID == 0){
		//alloacate arrays for maeijster distanc
		g = new int[m*n];
		dt = new int[m*n];
	}
	#pragma omp barrier

	//first phase
	for(int x = threadID; x < m; x+=NTHREADS){
		if(image->data[(x+0*m)*3] == 255)
			g[x+0*m] = 0;
		else
			g[x+0*m] = infinity;
		//1 scan
		for(int y = 1; y < n; y++){
			if(image->data[(x+y*m)*3] == 255)
				g[x+y*m] = 0;
			else
				g[x+y*m] = 1 + g[x+(y-1)*m];
		}
		//Scan 2
		for(int y = n-2; y >= 0; y--){
			if(g[x+(y+1)*m] < g[x+y*m])
				g[x+y*m] = 1 + g[x+(y+1)*m];
		}
	}

	//second phase
	#pragma omp barrier
	int * s = new int[m];
	int * t = new int[m];
	for(int y = threadID; y < n; y+=NTHREADS){
		int q = 0;
		int w  = 0;
		s[0] = 0;
		t[0] = 0;
		//Scan 3
		for(int u = 1; u < m; u++){
			while(q >= 0 && f(t[q], s[q], g[s[q]+y*m]) > f(t[q], u, g[u+y*m]))
				q--;
			if(q < 0){
				q = 0;
				s[0] = u;
			}
			else{
				w = 1 + sep(s[q], u, g[s[q]+y*m], g[u+y*m]);
				if(w < m){
					q++;
					s[q] = u;
					t[q] = w;
				}
			}
		}
		//scan 4
		for(int u = m-1; u >= 0; u--){
			int d = f(u, s[q], g[s[q]+y*m]);
			dt[u+y*m] = d;
			if(d <= k)
				image->data[(u+y*m)*3] = image->data[(u+y*m)*3+1] = image->data[(u+y*m)*3+2] = 255;
			else
				image->data[(u+y*m)*3] = image->data[(u+y*m)*3+1] = image->data[(u+y*m)*3+2] = 0;
			if(u == t[q])
				q--;
		}
	}
	#pragma omp barrier
	if(threadID == 0){
		free(dt);
		free(g);
	}
	free(s);
	free(t);
}
void getContour(Mat * image, int k, int threadID){
	//one of the threads is responsible for allocation
	if(threadID == 0){
		//alloacate arrays for maeijster distanc
		g = new int[m*n];
		dt = new int[m*n];
	}
	#pragma omp barrier

	//first phase
	for(int x = threadID; x < m; x+=NTHREADS){
		if(image->data[(x+0*m)*3] == 0)
			g[x+0*m] = 0;
		else
			g[x+0*m] = infinity;
		//1 scan
		for(int y = 1; y < n; y++){
			if(image->data[(x+y*m)*3] == 0)
				g[x+y*m] = 0;
			else
				g[x+y*m] = 1 + g[x+(y-1)*m];
		}

		//Scan 2
		for(int y = n-2; y >= 0; y--){
			if(g[x+(y+1)*m] < g[x+y*m])
				g[x+y*m] = 1 + g[x+(y+1)*m];
		}

	}

	#pragma omp barrier
	//second phase
	int * s = new int[m];
	int * t = new int[m];

	for(int y = threadID; y < n; y+=NTHREADS){
		int q = 0;
		int w  = 0;
		s[0] = 0;
		t[0] = 0;
		//Scan 3
		for(int u = 1; u < m; u++){
			while(q >= 0 && f(t[q], s[q], g[s[q]+y*m]) > f(t[q], u, g[u+y*m]))
				q--;
			if(q < 0){
				q = 0;
				s[0] = u;
			}
			else{
				w = 1 + sep(s[q], u, g[s[q]+y*m], g[u+y*m]);
				if(w < m){
					q++;
					s[q] = u;
					t[q] = w;
				}
			}
		}

		//scan 4
		for(int u = m-1; u >= 0; u--){
			int d = f(u, s[q], g[s[q]+y*m]);
			dt[u+y*m] = d;
			if(d != k)
				image->data[(u+y*m)*3] = image->data[(u+y*m)*3+1] = image->data[(u+y*m)*3+2] = 0;
			else
				image->data[(u+y*m)*3] = image->data[(u+y*m)*3+1] = image->data[(u+y*m)*3+2] = 255;
			if(u == t[q])
				q--;
		}
	}
	//cout << "hello" << endl;
	#pragma omp barrier
	if(threadID == 0){
		free(dt);
		free(g);
	}
	free(s);
	free(t);
}
void imgInRange(Mat * image, int LowH, int LowS, int LowV, int HighH, int HighS, int HighV, int threadID){
	int start = (int)( ( (image->rows * image->cols * 3) / (double)NTHREADS)*threadID);
	int stop = (int)( ( (image->rows * image->cols * 3) / (double)NTHREADS)*(threadID+1));

	for(int i = start; i < stop; i+=3){
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

void calMoments(Mat * image, int threadID){

	int start = (int)( ( (image->rows * image->cols * 3) / (double)NTHREADS)*threadID);
	int stop = (int)( ( (image->rows * image->cols * 3) / (double)NTHREADS)*(threadID+1));

	int t_zerothMoment = 0;
	for(int i = start; i < stop; i+=3)
		if(image->data[i] == 255)
			t_zerothMoment += 1;
	//sinchronize threads and write to the global value
	omp_set_lock(&lock);
	zerothMoment += t_zerothMoment;
	omp_unset_lock(&lock);

	int t_moment10 = 0;
	int t_moment01 = 0;
	for(int i = start; i < stop; i+=3){
		if(image->data[i] == 255){
			t_moment01 += (i/3)%image->cols;
			t_moment10 += (i/3)/image->cols;
		}
	}
	omp_set_lock(&lock);
	moment10 += t_moment10;
	moment01 += t_moment01;
	omp_unset_lock(&lock);
}

double ** int2double(int** a, double** ret, int rows, int cols) {
	int i = 0;
	for (; i < rows; i++) {
		int j = 0;
		for (; j < cols; j++) {
			ret[i][j] = (double)(a[i][j]);
		}
	}

	return ret;
}

FuncArgs createGenericArgs(){
	FuncArgs fa;
	fa.tfCols = tfCols;
	fa.tfRows = tfRows;
	fa.diag = diag;
	fa.numThetas = numThetas;
	fa.minNumPixelsOnLine = minNumOfPixelsOnLine;
	fa.orgRows = orgRows;
	fa.orgCols = orgCols;

	return fa;
}

timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}
int main( int argc, char** argv ){



	VideoCapture cap(0);
	cap.set(CV_CAP_PROP_FRAME_WIDTH,180);
	cap.set(CV_CAP_PROP_FRAME_HEIGHT,120);
	//declare Mat image files
	Mat imgOriginal, imgThresholded;
	cap.read(imgOriginal);
	m = imgOriginal.cols;
	n = imgOriginal.rows;
	infinity = m+n;
	omp_set_num_threads(NTHREADS);
	if ( !cap.isOpened() )  // if not success, exit program
	{
		 cout << "Cannot open the web cam" << endl;
		 return -1;
	}

	namedWindow("Control", CV_WINDOW_AUTOSIZE); //create a window called "Control"

	int iLowH = 0;
	int iHighH = 255;

	int iLowS = 105;
	int iHighS = 193;

	int iLowV = 189;
	int iHighV = 240;

	//Create trackbars in "Control" window
	cvCreateTrackbar("LowH", "Control", &iLowH, 255); //Hue (0 - 255)
	cvCreateTrackbar("HighH", "Control", &iHighH, 255);

	cvCreateTrackbar("LowS", "Control", &iLowS, 255); //Saturation (0 - 255)
	cvCreateTrackbar("HighS", "Control", &iHighS, 255);

	cvCreateTrackbar("LowV", "Control", &iLowV, 255); //Value (0 - 255)
	cvCreateTrackbar("HighV", "Control", &iHighV, 255);
	omp_init_lock(&lock);
	//Hough transform paramters
	int** img;
	double ** fixedEdges;
	double** blurHtr;
	double** orgImg;
	double** htr;
	double step;
	orgRows = imgOriginal.rows;
	orgCols = imgOriginal.cols;
	FuncArgs threadArgs[numThreads];
	FuncArgs linArgs = createGenericArgs();
	diag = sqrt((double)(orgRows * orgRows + orgCols * orgCols));
	tfRows = 2 + ((int)(2 * diag));
	tfCols = 2 + numThetas;

	blurHtr = makeDoublearray(tfRows, tfCols);
	htr = makeDoublearray(tfRows, tfCols);
	img=makeIntarray(orgRows, orgCols);
	orgImg = makeDoublearray(orgRows, orgCols);


	#pragma omp parallel
	{
		//get id of thread
		int threadID = omp_get_thread_num();

		while(true)
		{
			//make one thread capture images

			#pragma omp single
			{
				zerothMoment = 0;
				moment01 = 0;
				moment10 = 0;
				bool bSuccess = cap.read(imgOriginal); // read a new frame from video
				if (!bSuccess) //if not success, break loop
				{
				   cout << "Cannot read a frame from video stream" << endl;
				 // break;
				}
				imgThresholded = imgOriginal.clone();
			}
			#pragma omp barrier
			RgbToHsv(&imgOriginal, &imgThresholded, threadID);
			#pragma omp barrier
			//get thresholded image
			imgInRange(&imgThresholded, iLowH, iLowS, iLowV, iHighH, iHighS, iHighV, threadID);
			//morphological opening (remove small objects from the foreground)
			erode(&imgThresholded, 3, threadID);
			diliate(&imgThresholded, 3, threadID);
			//morphological closing (fill small holes in the foreground)
			diliate(&imgThresholded, 3, threadID);
			getContour(&imgThresholded, 3, threadID);


			calMoments(&imgThresholded, threadID);
			#pragma omp barrier
			//let thread 1 draw
			#pragma omp single
			{
				if(zerothMoment > 100){

					int xCor = moment01 / zerothMoment;
					int yCor= moment10 / zerothMoment;

					circle(imgOriginal, Point(xCor, yCor), 3.0, Scalar( 0, 0, 255 ), 6, 8 );
				}
			}

			#pragma omp barrier

#pragma omp single
			{

				get2dFrom1d(img, imgThresholded.data, orgRows, orgCols);
				orgImg = int2double(img, orgImg, orgRows, orgCols);
				//zeroth thread does all linear work
				zeroArray(htr, tfRows, tfCols);
				step = 1e-10 + (double)orgRows / (double)numThreads;

			}

#pragma omp barrier
			threadArgs[threadID] = createGenericArgs();
			threadArgs[threadID].arrIn = orgImg;
			threadArgs[threadID].arrOut = htr;
			threadArgs[threadID].threadStartRow = (int)((double)threadID*step);
			threadArgs[threadID].threadEndRow = (int)((double)(threadID + 1)*step);

			houghField(&threadArgs[threadID]);

#pragma omp barrier


#pragma omp single
			{
				//changing step size. we have been bound to the
				step = 1e-10 + (double)tfRows / (double)numThreads;

				linArgs.arrIn = htr;
				linArgs.threadStartRow = 0;
				linArgs.threadEndRow = linArgs.tfRows;
				//fixing edges of created transformation
				fixEdges(&linArgs);
				fixedEdges = linArgs.arrOut;
				zeroArray(blurHtr, tfRows, tfCols);
			}
#pragma omp barrier


			threadArgs[threadID] = createGenericArgs();
			threadArgs[threadID].arrIn = fixedEdges;
			threadArgs[threadID].arrOut = blurHtr;

			threadArgs[threadID].threadStartRow = (int)((double)threadID*step);
			threadArgs[threadID].threadEndRow = (int)((double)(threadID + 1)*step);

			//blurring the transformation to dampen the smaller peaks that may have occurred
			//blur(&hough[i]);
			blur(&threadArgs[threadID]);

#pragma omp barrier


#pragma omp single
			{
				linArgs = createGenericArgs();
				linArgs.threadStartRow = 0;
				linArgs.threadEndRow = tfRows;
				linArgs.arrIn = blurHtr;
				linArgs.arrOut = htr;
				getDesiredLines(&linArgs);

				mergeLines(&linArgs);
				rotateLines(&linArgs);
				mergeLines(&linArgs);
				linArgs.lines = removeUndesiredLines(&linArgs);
				//outputting gotten lines
				sortLines(&linArgs);

				polar2cartesian(&linArgs);
				outLines(&linArgs);

				for(int i=0; i<linArgs.numLines;i++)
					free(linArgs.lines[i]);
				free(linArgs.lines);

			}
#pragma omp barrier
			//show the image
#pragma omp single
			{
				imshow("Original",imgOriginal);
				imshow("Processed", imgThresholded);
				waitKey(10);

			}
		}
	}

	return 0;
}

