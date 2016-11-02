/*
 * DisplayImage.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: aljaz
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/core/core.hpp>
#include "opencv2/imgcodecs.hpp"
#include <opencv2/highgui/highgui.hpp>


using namespace std;
using namespace cv;


void RgbToHsv(Mat * OrigImage, Mat * HSVImage)
{

    unsigned char s,v;
    unsigned char rgbMin, rgbMax;
    *HSVImage = OrigImage->clone();
    for(int i = 0; i < OrigImage->rows * OrigImage->cols * 3; i+=3){

    	//finds the minimum value between r,g,b
    	// i = b, i+1 = g, i+2 = r
        rgbMin = OrigImage->data[i+2] < OrigImage->data[i+1]  ? (OrigImage->data[i+2]< OrigImage->data[i] ? OrigImage->data[i+2] : OrigImage->data[i]) : (OrigImage->data[i+1] < OrigImage->data[i] ? OrigImage->data[i+1] : OrigImage->data[i]);
        rgbMax = OrigImage->data[i+2] > OrigImage->data[i+1] ? (OrigImage->data[i+2] > OrigImage->data[i] ? OrigImage->data[i+2] : OrigImage->data[i]) : (OrigImage->data[i+1] > OrigImage->data[i] ? OrigImage->data[i+1] : OrigImage->data[i]);

        v = rgbMax;
        //if color is black
        if (v == 0)
        {
            HSVImage->data[i] = 0;//h
            HSVImage->data[i+1] = 0;//s
            HSVImage->data[i+2] = v;//v
            continue;
        }
        //saturation
        //if color is white
        s = 255* long(rgbMax - rgbMin) / v;

        if (s == 0)
        {
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

        HSVImage->data[i+1] = s;//s
        HSVImage->data[i+2] = v;//v
    }

}
int mod(int stevilo, int modus){
	if(stevilo == 0)
		return 0;
	if(abs(stevilo) < modus)
		return modus-abs(stevilo);
	else
		return abs(stevilo%modus);
}

void imgToIntArray(Mat * image, int * manArray){

	for(int i = 0; i < image->rows*image->cols*3; i++){
		//cout << "i: "<< i << " data: " << (int)image->data[i] << endl;
		manArray[i] = (int)image->data[i];
	}
}
int * manhattan3(int inRow, int inCol, Mat * image, int * manhattanArray){
	imgToIntArray(image, manhattanArray);

/*
	for(int k = 0; k < inRow*inCol*3; k+=3){
		//cout << "k: "<< (int)k << "| " ;
		if(k%(inRow*3) == 0)
			cout << endl ;
		cout << setfill(' ') << setw(3) << manhattanArray[k] << " ";
		//cout << (int)image->data[k+1] << " ";
		//cout << (int)image->data[k+2] << " ";
    }

    cout << endl<<endl;*/

	for(int i = 0; i < inRow*inCol*3; i+=3){
		//ce smo najdli bel pixel
		//cout << "i: "<< (int)i << "data: "<<(int)image->data[i];
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

/*
	for(int k = 0; k < inRow*inCol*3; k+=3){
		//cout << "k: "<< (int)k << "| " ;
		if(k%(inRow*3) == 0)
			cout << endl ;
		cout << setfill(' ') << setw(3) << manhattanArray[k] << " ";
		//cout << (int)image->data[k+1] << " ";
		//cout << (int)image->data[k+2] << " ";
    }

    cout << endl<<endl;*/

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
	}/*
	for(int k = 0; k < inRow*inCol*3; k+=3){
		//cout << "k: "<< (int)k << "| " ;
		if(k%(inRow*3) == 0)
			cout << endl ;
		cout << setfill(' ') << setw(3) << manhattanArray[k] << " ";
    }

    cout << endl<<endl;*/
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

int main( int argc, char** argv ){
/*
  Mat image, HSVimage, imgThresholded;
  image = imread( "/home/aljaz/Pictures/circle3.png", 1 );

 int iLowH = 0;
 int iHighH = 0;
 int iLowS = 0;
 int iHighS = 0;

 int iLowV = 255;
 int iHighV = 255;

	//cout << image;
  HSVimage = image.clone();

  //cvtColor(image, HSVimage , COLOR_BGR2HSV); //Convert the captured frame from BGR to HSV
  //cout << HSVimage << endl;
  RgbToHsv(&image, &HSVimage);
  imgThresholded = image.clone();
 // inRange(HSVimage, Scalar(iLowH, iLowS, iLowV), Scalar(iHighH, iHighS, iHighV), imgThresholded);

  imgInRange(&imgThresholded, iLowH, iLowS, iLowV, iHighH, iHighS, iHighV);

 // diliate(&imgThresholded, 2);
  //dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_RECT, Size(5, 5)) );
  erode(&imgThresholded, 9);
  if((!image.data))
    {
      printf( "No image data \n" );
      return -1;
    }

  namedWindow( "Display Image", CV_WINDOW_AUTOSIZE );
  namedWindow( "Display Image3", CV_WINDOW_AUTOSIZE );
  //imshow( "Display Image5", imgThresholded);
  imshow( "Display Image", HSVimage);
  imshow( "Display Image3", imgThresholded);*/

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
 cvCreateTrackbar("LowH", "Control", &iLowH, 179); //Hue (0 - 179)
 cvCreateTrackbar("HighH", "Control", &iHighH, 179);

  cvCreateTrackbar("LowS", "Control", &iLowS, 255); //Saturation (0 - 255)
 cvCreateTrackbar("HighS", "Control", &iHighS, 255);

  cvCreateTrackbar("LowV", "Control", &iLowV, 255); //Value (0 - 255)
 cvCreateTrackbar("HighV", "Control", &iHighV, 255);
  while (true)
  {
      Mat imgOriginal;

      bool bSuccess = cap.read(imgOriginal); // read a new frame from video

       if (!bSuccess) //if not success, break loop
      {
           cout << "Cannot read a frame from video stream" << endl;
           break;
      }

  Mat imgHSV = imgOriginal.clone();

 //vtColor(imgOriginal, imgHSV, COLOR_BGR2HSV); //Convert the captured frame from BGR to HSV
 RgbToHsv(&imgOriginal, &imgHSV);
 Mat imgThresholded = imgHSV.clone();

 //inRange(imgHSV, Scalar(iLowH, iLowS, iLowV), Scalar(iHighH, iHighS, iHighV), imgThresholded); //Threshold the image
 imgInRange(&imgThresholded, iLowH, iLowS, iLowV, iHighH, iHighS, iHighV);
//morphological opening (remove small objects from the foreground)
 erode(imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
 //dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
 erode(&imgThresholded, 2);
 diliate(&imgThresholded, 2);

 //morphological closing (fill small holes in the foreground)
//dilate( imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );
//erode(imgThresholded, imgThresholded, getStructuringElement(MORPH_ELLIPSE, Size(5, 5)) );

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

