/*
 * Tracking.cpp
 *
 *  Created on: Nov 13, 2016
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
#include <pthread.h>
#include <queue>

using namespace std;
using namespace cv;

#define NTHREADS 4
//thread variables
pthread_t threads[NTHREADS];
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t barr;
int inLine = 0;
//global capture object, opens camera 1
VideoCapture cap(0);
Mat imgOriginal;
Mat imgHSV;
Mat imgThresholded;

//global control pannel variables
int iLowH = 0;
int iHighH = 0;
int iLowS = 255;
int iHighS = 255;
int iLowV = 255;
int iHighV = 255;

//global calMoments variables
int zerothMoment = 0;
int moment10 = 0;
int moment01 = 0;
//manhattan distance variables
int m;
int n;
int * g ;
int infinity = 0;
//koncen array razdalj
int * dt;

void RgbToHsv(Mat * OrigImage, Mat * HSVImage, int rank){

	int start = (int)( ( (OrigImage->rows * OrigImage->cols * 3) / (double)NTHREADS)*rank);
	int stop = (int)( ( (OrigImage->rows * OrigImage->cols * 3) / (double)NTHREADS)*(rank+1));

	unsigned char s,v;
	unsigned char rgbMin, rgbMax;
	for(int i = start; i < stop; i+=3){

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

void imgInRange(Mat * image, int LowH, int LowS, int LowV, int HighH, int HighS, int HighV, int rank){

	int start = (int)( ( (image->rows * image->cols * 3) / (double)NTHREADS)*rank);
	int stop = (int)( ( (image->rows * image->cols * 3) / (double)NTHREADS)*(rank+1));


	for(int i= start; i < stop; i+=3){

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
void erode(Mat * image, int k, int rank){
	//if(rank == 0){
	if(rank == 0){
		m = image->cols;
		n = image->rows;
		infinity = image->cols+image->rows;
		dt = new int[m*n];
		g = new int[m*n];
	}

	pthread_barrier_wait(&barr);
	for(int x = rank; x < m; x+=NTHREADS){
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
	pthread_barrier_wait(&barr);
	//second phase
	int * s = new int[m];
	int * t = new int[m];
	for(int y = rank; y < n; y+=NTHREADS){

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
}
void diliate(Mat * image, int k, int rank){
	//if(rank == 0){
	if(rank == 0){
		m = image->cols;
		n = image->rows;
		infinity = image->cols+image->rows;
		dt = new int[m*n];
		g = new int[m*n];
	}

	pthread_barrier_wait(&barr);
	for(int x = rank; x < m; x+=NTHREADS){
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
	pthread_barrier_wait(&barr);
	//second phase
	int * s = new int[m];
	int * t = new int[m];
	for(int y = rank; y < n; y+=NTHREADS){
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
}
void calMoments(Mat * image, int rank){

	int start = (int)( ( (image->rows * image->cols * 3) / (double)NTHREADS)*rank);
	int stop = (int)( ( (image->rows * image->cols * 3) / (double)NTHREADS)*(rank+1));

	int t_zerothMoment = 0;
	for(int i = start; i < stop; i+=3)
		if(image->data[i] == 255)
			t_zerothMoment += 1;
	//sinchronize threads and write to the global value
	pthread_mutex_lock(&lock);
	zerothMoment += t_zerothMoment;
	pthread_mutex_unlock(&lock);

	int t_moment10 = 0;
	int t_moment01 = 0;
	for(int i = start; i < stop; i+=3){
		if(image->data[i] == 255){
			t_moment01 += (i/3)%image->cols;
			t_moment10 += (i/3)/image->cols;
		}
	}
	pthread_mutex_lock(&lock);
	moment10 += t_moment10;
	moment01 += t_moment01;
	pthread_mutex_unlock(&lock);
}
void getContours(Mat * image, int rank) {
	if(rank == 0){
	 //go from top to bottom
		bool first = true;
			for (int j= 0; j < image->cols*3; j+=3) {
			first = true;
			for (int i = 0; i < image->rows; i+=3) {
				if (first == true && (image->data[i*image->cols+j] == 255 || image->data[i*image->cols+j] == 2)) {
					image->data[i*image->cols+j] = 2;
					first = false;
				}
			}
		}
	}
 //go from left to right
	if(rank == 1){
		for (int i = 0; i < image->rows*3; i+=3) {
		bool first = true;
			for (int j = 0; j < image->cols*3; j+=3) {
				if (first == true && (image->data[i*image->cols+j] == 255 || image->data[i*image->cols+j] == 2)) {
					image->data[i*image->cols+j] = 2;
					first = false;
				}
			}
		}
	}
 //go from bottom to top
	if(rank == 2){
		for (int j = image->cols*3-3; j >= 0; j-=3) {
			bool first = true;
			for (int i = image->rows*3-3; i >= 0; i-=3) {
				//cout << "checking i: " << i << " and j: " << j <<  endl;
				if (first == true && (image->data[i*image->cols+j] == 255 || image->data[i*image->cols+j] == 2)) {
					image->data[i*image->cols+j] = 2;
					first = false;
				}
			}
		}
	}
 //go right to left
	if(rank == 3){
		for (int i = image->rows*3-3; i >= 0; i-=3) {
			bool first = true;
			for (int j = image->cols*3-3; j >= 0; j-=3) {
				if (first == true && (image->data[i*image->cols+j] == 255 || image->data[i*image->cols+j] == 2)) {
					image->data[i*image->cols+j] = 2;
					first = false;
				}
			}
		}
	}
	if(rank == 0){
 //create contours
		for (int i = 0; i < image->rows*3; i+=3) {
			for (int j = 0; j < image->cols*3; j+=3) {
				if ( image->data[i*image->cols+j] != 2 ) {
					image->data[i*image->cols+j] = 0;
					image->data[i*image->cols+j+1] = 0;
					image->data[i*image->cols+j+2] = 0;
				}
				else{
					image->data[i*image->cols+j] = 255;
					image->data[i*image->cols+j+1] = 255;
					image->data[i*image->cols+j+2] = 255;
				}
			}
		}
	}

}
void * processImage(void * arg){

	int myRank = (int)arg;
	//Image Objects

	while (true){
		//pthread_mutex_lock(&lock);
		if(myRank == 0){
			zerothMoment = 0;
			moment01 = 0;
			moment10 = 0;
			bool bSuccess = cap.read(imgOriginal); // read a new frame from video
			//imgOriginal = imread( "triangle.png", 1 );
			if (!bSuccess) //if not success, break loop
			{
			   cout << "Cannot read a frame from video stream" << endl;
			   break;
			}
		}

		if(myRank == 0)
			imgHSV = imgOriginal.clone();
		//convert rgb to HSV
		pthread_barrier_wait(&barr);
		RgbToHsv(&imgOriginal, &imgHSV, myRank);
		//convert ot threshold
		pthread_barrier_wait(&barr);
		if(myRank == 0)
			imgThresholded = imgHSV.clone();
		pthread_barrier_wait(&barr);

		imgInRange(&imgThresholded, iLowH, iLowS, iLowV, iHighH, iHighS, iHighV, myRank);
		//deklaracija spremenljivk za dilate

		erode(&imgThresholded, 3, myRank);
		diliate(&imgThresholded, 3, myRank);

		pthread_barrier_wait(&barr);
		diliate(&imgThresholded, 3, myRank);
		erode(&imgThresholded, 3, myRank);

		pthread_barrier_wait(&barr);
		calMoments(&imgThresholded, myRank);
		pthread_barrier_wait(&barr);
		getContours(&imgThresholded, myRank);
		pthread_barrier_wait(&barr);
		if(myRank == 0 && zerothMoment > 1000){

			int xCor = moment01 / zerothMoment;
			int yCor= moment10 / zerothMoment;

			circle(imgOriginal, Point(xCor, yCor), 3.0, Scalar( 0, 0, 255 ), 6, 8 );
		}
		//pthread_mutex_lock(&lock);
		if(myRank == 0){


			imshow("Original", imgOriginal); //show the original image
			imshow("Thresholded", imgThresholded);
		}
		//imshow("Thresholded", *threshQueue.front());
		//if(myRank == 0)
			//HSVQueue.pop();
		//threshQueue.pop();
		//imshow("Thresholded", imgThresholded);
		if(myRank == 0){
			if (waitKey(30) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
			{
			  cout << "esc key is pressed by user" << endl;
			  break;
			}
		}


	}

	return NULL;
}
int main(){

    pthread_barrier_init(&barr, NULL, NTHREADS);
/*
	//imgOriginal = imread( "triangle.png", 1 );
	pthread_t t[NTHREADS];

	if((!imgOriginal .data)){
	  printf( "No image data \n" );
	  return -1;
	}

	imgHSV = imgOriginal.clone();
	//convert rgb to HSV
	for(int i = 0; i < NTHREADS; i++)
		pthread_create(&t[i], NULL, RgbToHsv, (void *)i);

	for (int i = 0; i<NTHREADS; i++)
		pthread_join(t[i], NULL);
	//finsehd converting rgb to HSV
	imgThresholded = imgHSV.clone();
	//make thresholded image
	for(int i = 0; i < NTHREADS; i++)
		pthread_create(&t[i], NULL, imgInRange, (void *)i);

	for (int i = 0; i<NTHREADS; i++)
		pthread_join(t[i], NULL);
	//finsehd making thresol image

	//convert char array to int array
	manhattanArray = new int[imgThresholded.rows*imgThresholded.cols*3];
	for(int i = 0; i < NTHREADS; i++)
		pthread_create(&t[i], NULL, imgToIntArray, (void *)i);

	for (int i = 0; i<NTHREADS; i++)
		pthread_join(t[i], NULL);

	//finish coverting char array to int array
	//erode(&imgThresholded, 3);
	//diliate(&imgThresholded, 3);

	//calMoments
	for(int i = 0; i < NTHREADS; i++)
		pthread_create(&t[i], NULL, calMoments, (void *)i);

	for (int i = 0; i<NTHREADS; i++)
		pthread_join(t[i], NULL);

	if(zerothMoment > 1000){
		int xCor = moment01 / zerothMoment;
		int yCor= moment10 / zerothMoment;
		circle(imgOriginal, Point(xCor, yCor), 32.0, Scalar( 0, 255, 0 ), 1, 8 );
	}
	//calMoments(&imgThresholded, &imgOriginal);
	//printArray();
	imshow("Original", imgOriginal); //show the original image
	imshow("Thresholded", imgThresholded);
*/


	namedWindow("Control", CV_WINDOW_AUTOSIZE); //create a window called "Control"


	//Create trackbars in "Control" window
	cvCreateTrackbar("LowH", "Control", &iLowH, 255); //Hue (0 - 255)
	cvCreateTrackbar("HighH", "Control", &iHighH, 255);

	cvCreateTrackbar("LowS", "Control", &iLowS, 255); //Saturation (0 - 255)
	cvCreateTrackbar("HighS", "Control", &iHighS, 255);

	cvCreateTrackbar("LowV", "Control", &iLowV, 255); //Value (0 - 255)
	cvCreateTrackbar("HighV", "Control", &iHighV, 255);
	//capture the video from web cam
	if ( !cap.isOpened() )  // if not success, exit program
	{
		 cout << "Cannot open the web cam" << endl;
		 return -1;
	}

	for(int i = 0; i < NTHREADS; i++){
		pthread_create(&threads[i], NULL, processImage, (void *)i);
	}

	for (int i = 0; i<NTHREADS; i++)
		pthread_join(threads[i], NULL);

	return 0;
}
