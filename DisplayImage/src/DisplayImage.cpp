/*
 * DisplayImage.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: aljaz
 */

#include <iostream>
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

int main( int argc, char** argv ){
  Mat image, HSVimage, HSVimage2;
  image = imread( "/home/aljaz/image2.jpg", 1 );

  HSVimage = image.clone();
  HSVimage2 = image.clone();
  cvtColor(image, HSVimage2 , COLOR_BGR2HSV); //Convert the captured frame from BGR to HSV
  //cout << HSVimage << endl;
  RgbToHsv(&image, &HSVimage);
  if( (argc != 2) | (!image.data))
    {
      printf( "No image data \n" );
      return -1;
    }
  namedWindow( "Display Image", CV_WINDOW_AUTOSIZE );
  namedWindow( "Display Image2", CV_WINDOW_AUTOSIZE );
  namedWindow( "Display Image3", CV_WINDOW_AUTOSIZE );
  imshow( "Display Image3", HSVimage2);
  imshow( "Display Image2", image);
  imshow( "Display Image", HSVimage);

  waitKey(0);

  return 0;
}

