/*
 * DisplayImage.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: aljaz
 */

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

int * manhattan3(int inRow, int inCol, Mat * image, int * manhattanArray){
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
	}

    return manhattanArray;
}
int * manhattan2(int inRow, int inCol, Mat * image, int * manhattanArray){
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
	}
    return manhattanArray;
}
void diliate(Mat * image, int k){

	int manhattanArray[image->rows*image->cols*3];

	manhattan3(image->cols, image->rows, image, manhattanArray);
	for(int i = 0; i < image->rows*image->cols*3; i++){
		if(manhattanArray[i] <= k){
			image->data[i] = 255;
		}
		else{
			image->data[i] = 0;

		}
	}
}
void erode(Mat * image, int k){

	int manhattanArray[image->rows*image->cols*3];
	manhattan2(image->cols, image->rows, image, manhattanArray);
	for(int i = 0; i < image->rows*image->cols*3; i++){
		if(manhattanArray[i] <= k)
			image->data[i] = 0;
		else
			image->data[i] = 255;
	}

}
void imgInRange(Mat * image, int LowH, int LowS, int LowV, int HighH, int HighS, int HighV){
    //cout << "helo"<< endl;
	for(int i = 0; i < image->rows*image->cols*3; i+=3){
		//cout << "i: " << i << endl;
		//cout << (int)image->data[i] << " " <<  (int)image->data[i+1] << " "<< (int)image->data[i+2] << endl;
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

	int xCor = moment01 / zerothMoment;
	int yCor= moment10 / zerothMoment;


	circle(*origImage, Point(xCor, yCor), 32.0, Scalar( 0, 255, 0 ), 1, 8 );
}
int main( int argc, char** argv ){
/*
  Mat image, HSVimage, imgThresholded;
  image = imread( "/home/aljaz/Pictures/triangle.png", 1 );

 int iLowH = 0;
 int iHighH = 85;
 int iLowS = 255;
 int iHighS = 255;

 int iLowV = 255;
 int iHighV = 255;


  if((!image.data))
    {
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
 erode(&imgThresholded, 2);
 diliate(&imgThresholded, 2);

 //morphological closing (fill small holes in the foreground)
 diliate(&imgThresholded, 2);
 erode(&imgThresholded, 2);

 calMoments(&imgThresholded, &image);
  //imshow( "Display Image5", imgThresholded);
  imshow("Thresholded Image", imgThresholded); //show the thresholded image
  imshow("Original", image); //show the original image

	*/
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
	}
	waitKey(0);
	return 0;
}

