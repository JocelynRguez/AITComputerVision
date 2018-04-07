// Source file for image class



// Include files

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <algorithm>
#include <vector>
#include <iostream>
using namespace std;


////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////

vector<Feature> harrisFeatures;

R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width),
    height(image.height)

{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	// R2Point p1(1.2,3.5);
	// R2Point p2(2.1,2.2);
	// R2Point p3(0.2,1.6);
	// R2Point p4(0.0,0.25);
	// R2Point p5(-0.2,4.2);

  vector<R2Point> points;
  vector<R2Point> vectorA;
  vector<R2Point> vectorB;

  R2Point p1(0.0,0.0);
  R2Point p2(1.0,0.0);
  R2Point p3(1.0,1.0);
  R2Point p4(0.0,1.0);

  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);
  points.push_back(p4);

  R2Point a1(1.0, 0.0);
  R2Point a2(2.0, 0.0);
  R2Point a3(2.0, 1.0);
  R2Point a4(1.0, 1.0);


  vectorA.push_back(a1);
  vectorA.push_back(a2);
  vectorA.push_back(a3);
  vectorA.push_back(a4);

  R2Point b1(1.0, 2.0);
  R2Point b2(1.0, 1.0);
  R2Point b3(3.0, 1.0);
  R2Point b4(3.0, 2.0);


  vectorB.push_back(b1);
  vectorB.push_back(b2);
  vectorB.push_back(b3);
  vectorB.push_back(b4);



	// build the 5x6 matrix of equations
	//double** linEquations = dmatrix(1,5,1,6);

  //8X9 matrix
  double** linEquationsA = dmatrix(1,8,1,9);
  double** linEquationsB = dmatrix(1,8,1,9);


  for(int i = 0; i < 4; i++){
      //first row
      linEquationsA[2*i +1][1] = 0.0;
    	linEquationsA[2*i +1][2] = 0.0;
    	linEquationsA[2*i +1][3] = 0.0;
    	linEquationsA[2*i +1][4] = (-1.0)*points[i][0]; //-w*x_a
    	linEquationsA[2*i +1][5] = (-1.0)*points[i][1]; //-w*y_a
    	linEquationsA[2*i +1][6] = -1.0; //-w*w
      linEquationsA[2*i +1][7] = vectorA[i][1]*points[i][0]; //y_b*x_a
      linEquationsA[2*i +1][8] = vectorA[i][1]*points[i][1]; //y_b*y_a
      linEquationsA[2*i +1][9] = vectorA[i][1]; //y_b*w

      //second row
      linEquationsA[2*i +2][1] = points[i][0]; //w*x_a
      linEquationsA[2*i +2][2] = points[i][1]; //w*y_a
      linEquationsA[2*i +2][3] = 1.0; //w*w
      linEquationsA[2*i +2][4] = 0.0;
      linEquationsA[2*i +2][5] = 0.0;
      linEquationsA[2*i +2][6] = 0.0;
      linEquationsA[2*i +2][7] = -(vectorA[i][0]*points[i][0]); //-x_b*x_a
      linEquationsA[2*i +2][8] = -(vectorA[i][0]*points[i][1]); //-x_b*y_a
      linEquationsA[2*i +2][9] = (-1.0)*vectorA[i][0]; //-x_b*w

  }

  for(int i = 0; i < 4; i++){
      //first row
      linEquationsB[2*i +1][1] = 0.0;
    	linEquationsB[2*i +1][2] = 0.0;
    	linEquationsB[2*i +1][3] = 0.0;
    	linEquationsB[2*i +1][4] = -points[i][0]; //-w*x_a
    	linEquationsB[2*i +1][5] = -points[i][1]; //-w*y_a
    	linEquationsB[2*i +1][6] = -1.0; //-w*w
      linEquationsB[2*i +1][7] = vectorB[i][1]*points[i][0]; //y_b*x_a
      linEquationsB[2*i +1][8] = vectorB[i][1]*points[i][1]; //y_b*y_a
      linEquationsB[2*i +1][9] = vectorB[i][1]; //y_b*w

      //second row
      linEquationsB[2*i +2][1] = points[i][0]; //w*x_a
      linEquationsB[2*i +2][2] = points[i][1]; //w*y_a
      linEquationsB[2*i +2][3] = 1.0; //w*w
      linEquationsB[2*i +2][4] = 0.0;
      linEquationsB[2*i +2][5] = 0.0;
      linEquationsB[2*i +2][6] = 0.0;
      linEquationsB[2*i +2][7] = -vectorB[i][0]*points[i][0]; //-x_b*x_a
      linEquationsB[2*i +2][8] = -vectorB[i][0]*points[i][1]; //-x_b*y_a
      linEquationsB[2*i +2][9] = -vectorB[i][0]; //-x_b*w
  }

  // vector<int> A = dmatrix(1,8,1,9);
  // double** B = dmatrix(1,8,1,9);


  fprintf(stderr, "Matrix for A points: \n");
  for(int k = 1; k < 9; k++){
    fprintf(stderr, "%f %f %f %f %f %f %f %f %f\n", linEquationsA[k][1], linEquationsA[k][2], linEquationsA[k][3], linEquationsA[k][4], linEquationsA[k][5], linEquationsA[k][6], linEquationsA[k][7], linEquationsA[k][8], linEquationsA[k][9]);

  }

  fprintf(stderr, "Matrix for B points: \n");
  for(int k = 1; k < 9; k++){
    fprintf(stderr, "%f %f %f %f %f %f %f %f %f\n", linEquationsB[k][1], linEquationsB[k][2], linEquationsA[k][3], linEquationsB[k][4], linEquationsB[k][5], linEquationsB[k][6], linEquationsA[k][7], linEquationsB[k][8], linEquationsB[k][9]);

  }
	// linEquations[1][1] = p1[0]*p1[0];
	// linEquations[1][2] = p1[0]*p1[1];
	// linEquations[1][3] = p1[1]*p1[1];
	// linEquations[1][4] = p1[0];
	// linEquations[1][5] = p1[1];
	// linEquations[1][6] = 1.0;
  //
	// linEquations[2][1] = p2[0]*p2[0];
	// linEquations[2][2] = p2[0]*p2[1];
	// linEquations[2][3] = p2[1]*p2[1];
	// linEquations[2][4] = p2[0];
	// linEquations[2][5] = p2[1];
	// linEquations[2][6] = 1.0;
  //
	// linEquations[3][1] = p3[0]*p3[0];
	// linEquations[3][2] = p3[0]*p3[1];
	// linEquations[3][3] = p3[1]*p3[1];
	// linEquations[3][4] = p3[0];
	// linEquations[3][5] = p3[1];
	// linEquations[3][6] = 1.0;
  //
	// linEquations[4][1] = p4[0]*p4[0];
	// linEquations[4][2] = p4[0]*p4[1];
	// linEquations[4][3] = p4[1]*p4[1];
	// linEquations[4][4] = p4[0];
	// linEquations[4][5] = p4[1];
	// linEquations[4][6] = 1.0;
  //
	// linEquations[5][1] = p5[0]*p5[0];
	// linEquations[5][2] = p5[0]*p5[1];
	// linEquations[5][3] = p5[1]*p5[1];
	// linEquations[5][4] = p5[0];
	// linEquations[5][5] = p5[1];
	// linEquations[5][6] = 1.0;

	printf("\n Four points: \n");
	printf("Point #1: %f,%f\n",p1[0],p1[1]);
	printf("Point #2: %f,%f\n",p2[0],p2[1]);
	printf("Point #3: %f,%f\n",p3[0],p3[1]);
	printf("Point #4: %f,%f\n",p4[0],p4[1]);

	// compute the SVD
	//double singularValues[7]; // 1..6
  double singularValuesA[10]; // 1..9
  double singularValuesB[10]; // 1..9

	//double** nullspaceMatrix = dmatrix(1,6,1,6);
  double** nullspaceMatrixA = dmatrix(1,9,1,9);
  double** nullspaceMatrixB = dmatrix(1,9,1,9);
	svdcmp(linEquationsA, 8, 9, singularValuesA, nullspaceMatrixA);
  svdcmp(linEquationsB, 8, 9, singularValuesB, nullspaceMatrixB);

	// get the result
	fprintf(stderr, "\nSingular values A: %f, %f, %f, %f, %f, %f, %f, %f, %f\n", singularValuesA[1],singularValuesA[2],singularValuesA[3],singularValuesA[4],singularValuesA[5],singularValuesA[6], singularValuesA[7], singularValuesA[8], singularValuesA[9]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for(int i=2;i<10;i++) if(singularValuesA[i]<singularValuesA[smallestIndex]) smallestIndex=i;

  fprintf(stderr, "smallestIndexA %d\n", smallestIndex);

  // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  fprintf(stderr, "H Values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n\n",nullspaceMatrixA[1][smallestIndex],nullspaceMatrixA[2][smallestIndex],nullspaceMatrixA[3][smallestIndex],nullspaceMatrixA[4][smallestIndex],nullspaceMatrixA[5][smallestIndex],nullspaceMatrixA[6][smallestIndex],nullspaceMatrixA[7][smallestIndex],nullspaceMatrixA[8][smallestIndex],nullspaceMatrixA[9][smallestIndex]);


  fprintf(stderr, "H Matrix for A points: \n");
  for(int j = 1; j < 9; j = j+3){
      fprintf(stderr, "%f, %f, %f \n", nullspaceMatrixA[j][smallestIndex],nullspaceMatrixA[j
      +1][smallestIndex], nullspaceMatrixA[j+2][smallestIndex]);
  }


  // get the result
	fprintf(stderr, "\nSingular values B: %f, %f, %f, %f, %f, %f, %f, %f, %f\n", singularValuesB[1],singularValuesB[2],singularValuesB[3],singularValuesB[4],singularValuesB[5],singularValuesB[6], singularValuesB[7], singularValuesB[8], singularValuesB[9]);

	// find the smallest singular value:
	int smallestIndexB = 1;
	for(int i=2;i<10;i++) if(singularValuesB[i]<singularValuesB[smallestIndexB]) smallestIndexB=i;

  fprintf(stderr, "smallestIndexB: %d\n", smallestIndexB);

  fprintf(stderr, "H Values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n\n",nullspaceMatrixB[1][smallestIndexB],nullspaceMatrixB[2][smallestIndexB],nullspaceMatrixB[3][smallestIndexB],nullspaceMatrixB[4][smallestIndexB],nullspaceMatrixB[5][smallestIndexB],nullspaceMatrixB[6][smallestIndexB],nullspaceMatrixB[7][smallestIndexB],nullspaceMatrixB[8][smallestIndexB],nullspaceMatrixB[9][smallestIndexB]);

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  fprintf(stderr, "H Matrix for B points: \n");
  for(int j = 1; j < 10; j=j+3){
    	fprintf(stderr, "%f, %f, %f \n", nullspaceMatrixB[j][smallestIndexB],nullspaceMatrixB[j+1][smallestIndexB], nullspaceMatrixB[j+2][smallestIndexB]);
  }


//	make sure the solution is correct:
    // for(int n = 1; n < 9; n ++){
    //   printf("Equation #%d result: %f\n", n,	linEquationsA[n][1]*nullspaceMatrixA[1][smallestIndex] +
    // 	linEquationsA[n][2]*nullspaceMatrixA[2][smallestIndex] +
    // 	linEquationsA[n][3]*nullspaceMatrixA[3][smallestIndex] +
    // 	linEquationsA[n][4]*nullspaceMatrixA[4][smallestIndex] +
    // 	linEquationsA[n][5]*nullspaceMatrixA[5][smallestIndex] +
    // 	linEquationsA[n][6]*nullspaceMatrixA[6][smallestIndex] +
    //   linEquationsA[n][7]*nullspaceMatrixA[7][smallestIndex] +
    //   linEquationsA[n][8]*nullspaceMatrixA[8][smallestIndex] +
    //   linEquationsA[n][9]*nullspaceMatrixA[9][smallestIndex]);
    //
    // }

	// printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrixA[1][smallestIndex] +
	// 									p1[0]*p1[1]*nullspaceMatrixA[2][smallestIndex] +
	// 									p1[1]*p1[1]*nullspaceMatrixA[3][smallestIndex] +
	// 									p1[0]*nullspaceMatrixA[4][smallestIndex] +
	// 									p1[1]*nullspaceMatrixA[5][smallestIndex] +
	// 									nullspaceMatrixA[6][smallestIndex]);
  //
	// printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrixA[1][smallestIndex] +
	// 									p2[0]*p2[1]*nullspaceMatrixA[2][smallestIndex] +
	// 									p2[1]*p2[1]*nullspaceMatrixA[3][smallestIndex] +
	// 									p2[0]*nullspaceMatrixA[4][smallestIndex] +
	// 									p2[1]*nullspaceMatrixA[5][smallestIndex] +
	// 									nullspaceMatrixA[6][smallestIndex]);
  //
	// printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrixA[1][smallestIndex] +
	// 									p3[0]*p3[1]*nullspaceMatrixA[2][smallestIndex] +
	// 									p3[1]*p3[1]*nullspaceMatrixA[3][smallestIndex] +
	// 									p3[0]*nullspaceMatrixA[4][smallestIndex] +
	// 									p3[1]*nullspaceMatrixA[5][smallestIndex] +
	// 									nullspaceMatrixA[6][smallestIndex]);
  //
	// printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrixA[1][smallestIndex] +
	// 									p4[0]*p4[1]*nullspaceMatrixA[2][smallestIndex] +
	// 									p4[1]*p4[1]*nullspaceMatrixA[3][smallestIndex] +
	// 									p4[0]*nullspaceMatrixA[4][smallestIndex] +
	// 									p4[1]*nullspaceMatrixA[5][smallestIndex] +
	// 									nullspaceMatrixA[6][smallestIndex]);

	// printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrixA[1][smallestIndex] +
	// 									p5[0]*p5[1]*nullspaceMatrixA[2][smallestIndex] +
	// 									p5[1]*p5[1]*nullspaceMatrixA[3][smallestIndex] +
	// 									p5[0]*nullspaceMatrixA[4][smallestIndex] +
	// 									p5[1]*nullspaceMatrixA[5][smallestIndex] +
	// 									nullspaceMatrixA[6][smallestIndex]);

	// R2Point test_point(0.34,-2.8);
  //
	// printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] +
	// 										test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] +
	// 										test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] +
	// 										test_point[0]*nullspaceMatrix[4][smallestIndex] +
	// 										test_point[1]*nullspaceMatrix[5][smallestIndex] +
	// 										nullspaceMatrix[6][smallestIndex]);

	return;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();

    }
  }
}

void R2Image::
SobelX(void)
{
  (*this).ChangeSaturation(0.0);

  float kernelx[3][3] =  {{1, 0, -1},
                          {2, 0, -2},
                          {1, 0, -1}};

  R2Image tempImg(*this);
  R2Pixel* tempPixel = new R2Pixel();

	// Apply the Sobel oprator to the image in X direction
  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {

        *tempPixel = (kernelx[0][0]*Pixel(i-1, j-1)) +
                    (kernelx[0][1]*Pixel(i, j-1)) +
                    (kernelx[0][2]*Pixel(i+1, j-1)) +
                    (kernelx[1][0]*Pixel(i-1, j)) +
                    (kernelx[1][1]*Pixel(i, j)) +
                    (kernelx[1][2]*Pixel(i+1, j)) +
                    (kernelx[2][0]*Pixel(i-1, j+1)) +
                    (kernelx[2][1]*Pixel(i, j+1)) +
                    (kernelx[2][2]*Pixel(i+1, j+1));

      tempImg.Pixel(i,j) = *tempPixel;
     //Pixel(i,j).Clamp();

    }
  }

  for(int x = 0; x < width; x++){
    for(int y = 0; y < height; y++){
      Pixel(x,y) = tempImg.Pixel(x,y);
      //Pixel(x,y).Clamp();
    }

  }



}

void R2Image::
SobelY(void)
{

  (*this).ChangeSaturation(0.0);
	// Apply the Sobel oprator to the image in Y direction
  float kernely[3][3] = {{-1, -2, -1},
                          {0, 0, 0},
                          {1, 2, 1}};

  R2Image tempImg(*this);

  R2Pixel* tempPixel = new R2Pixel();

  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {

      //change tempImg
      *tempPixel = (kernely[0][0]*Pixel(i-1, j-1)) +
                  (kernely[0][1]*Pixel(i, j-1)) +
                  (kernely[0][2]*Pixel(i+1, j-1)) +
                  (kernely[1][0]*Pixel(i-1, j)) +
                  (kernely[1][1]*Pixel(i, j)) +
                  (kernely[1][2]*Pixel(i+1, j)) +
                  (kernely[2][0]*Pixel(i-1, j+1)) +
                  (kernely[2][1]*Pixel(i, j+1)) +
                  (kernely[2][2]*Pixel(i+1, j+1));
      //temp image
      tempImg.Pixel(i,j) = *tempPixel;
      //Pixel(i,j).Clamp();
    }
  }

  for(int x = 0; x < width; x++){
    for(int y = 0; y < height; y++){
      Pixel(x,y) = tempImg.Pixel(x,y);
      //Pixel(x,y).Clamp();
    }

  }


}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  double for_red = .3;
  double for_green = .6;
  double for_blue = .1;


  double sq;
  double red_val;
  double blue_val;
  double green_val;

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      red_val = Pixel(i,j).Red();
      blue_val = Pixel(i,j).Blue();
      green_val = Pixel(i,j).Green();

      sq = sqrt(red_val*red_val*for_red +
        blue_val*blue_val*for_blue +
        green_val*green_val*for_green);

      Pixel(i,j).SetRed(sq+(red_val-sq)*factor);
      Pixel(i,j).SetGreen(sq+(green_val-sq)*factor);
      Pixel(i,j).SetBlue(sq+(blue_val-sq)*factor);
    }
  }


}
// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  // Gaussian blur of the image. Separable solution is preferred
  R2Image tempImg(*this);

  //kernel width
  int kwidth = sigma*6 + 1;
  int ksize = sigma*3;
  //fprintf(stderr, "kwidth: (%d), ksize: (%d)\n", kwidth, ksize);

  //needs to ultimitalely add up to 1
  double weight = 0.0;
  double sum = 0.0;

  //need to create kernel
  double kernel[kwidth];
  double pi = 3.14159265358979323846;
  double square = pow(sigma, 2);

  //compute the weights
  for(int w = 0; w <= kwidth/2; w++){
    weight = (1/(sigma*sqrt(2*pi)))*(exp(-(pow((-kwidth/2+w),2)/(2*square))));
    kernel[w] = weight;
    kernel[kwidth-w-1] = weight;
    //kernel[kwidth-1-w] = weight;
  }

  //get sum of kernel weights
  for(int s = 0; s < kwidth; s++){
    sum += kernel[s];
  }

  //normalize kernel
  double endsum = 0.0;
  for(int l = 0; l < kwidth; l++){
    kernel[l] = kernel[l]/sum;
    endsum += kernel[l];
    //fprintf(stderr, "kernel(%d): (%f)\n", l, kernel[l]);
  }
  //fprintf(stderr, "end(%f)\n", endsum);


  //y direction
  for(int i = ksize; i < width-ksize; i++){
    for(int j = ksize; j < height-ksize; j++){
      //another iteration, over kernel??
      R2Pixel* tempPixel = new R2Pixel();
      for(int ly = -ksize; ly <= ksize; ly++){
        *tempPixel += Pixel(i,j+ly)*kernel[ly+ksize];
      }
      tempImg.Pixel(i,j) = *tempPixel;
      //tempImg.Pixel(i,j).Clamp();
    }
  }

  //x direction
  for(int i = ksize; i < width-ksize; i++){
    for(int j = ksize; j < height-ksize; j++){

      R2Pixel* tempPixel = new R2Pixel();
      for(int lx = -ksize; lx <= ksize; lx++){
        *tempPixel += tempImg.Pixel(i+lx,j)*kernel[lx+ksize];
      }
      Pixel(i,j) = *tempPixel;
      //Pixel(i,j).Clamp();
    }
  }




  //Look at Separable filter psuedo code

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
}


void R2Image::
Median(int median)
{
    // The main idea of the median filter is to run through the signal entry by entry, replacing each entry with the median of neighboring entries.
    //Blur(2.0);

    int size = median/2;
    R2Pixel neighbors[median];
    //double middle = 0.0;
    //double sorted[median];

    double red = 0.0;
    double green = 0.0;
    double blue = 0.0;
    double sum = 0.0;

    R2Image tempImg(*this);
    //y direction
    for(int x = 0; x < width; x++){
      for(int y = size; y< height-size; y++){
        //R2Pixel* tempPixel = new R2Pixel();
        for(int i = -size; i <= size; i++){
          red = Pixel(x, y+i).Red();
          green = Pixel(x, y+i).Green();
          blue = Pixel(x, y+i).Blue();
          sum = red+green+blue;

          neighbors[i+size] = Pixel(x, y+i);
          //fprintf(stderr, "Val: (%g)\n", neighbors[i]);
        }

        //sort through array

        //std::sort(neighbors, neighbors + median);
        tempImg.Pixel(x,y) = neighbors[size-1];
        //fprintf(stderr, "Val: (%g)\n", neighbors[median/2]);
      }
    }

    //x direction
    for(int x = median; x < width-median; x++){
      for(int y = 0; y< height; y++){


        for(int i = -size; i <= size; i++){
          red = tempImg.Pixel(x+i, y).Red();
          green = tempImg.Pixel(x+i, y).Green();
          blue = tempImg.Pixel(x+i, y).Blue();
          neighbors[i+size] = tempImg.Pixel(x+i, y);
          //fprintf(stderr, "Val: (%g)\n", neighbors[i]);
        }

        //std::sort(neighbors, neighbors + median);
        Pixel(x,y) = neighbors[size];
        //fprintf(stderr, "Val: (%g)\n", *neighbors);

      }
    }



  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "Median not implemented\n");
}


void R2Image::
Bilateral()
{
    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "Bilateral not implemented\n");
}


vector<Feature> R2Image::
filterHarris(vector<Feature> values){
    //int prev = bpixels[0];
    //double kernel[][] = {}

    int currX;
    int currY;
    int offset = 5;

    R2Image tempImg(*this);

  for(int index = 0; index < values.size(); index++){
    currX = values[index].centerX;
    currY = values[index].centerY;

      for(int s = -offset; s <= offset; s++){

        Pixel(currX+s,currY+offset).SetRed(0);
        Pixel(currX+s,currY+offset).SetGreen(1);
        Pixel(currX+s,currY+offset).SetBlue(0);

        Pixel(currX+offset,currY+s).SetRed(0);
        Pixel(currX+offset,currY+s).SetGreen(1);
        Pixel(currX+offset,currY+s).SetBlue(0);
        //
        Pixel(currX-s,currY-offset).SetRed(0);
        Pixel(currX-s,currY-offset).SetGreen(1);
        Pixel(currX-s,currY-offset).SetBlue(0);
        //
        Pixel(currX-offset,currY-s).SetRed(0);
        Pixel(currX-offset,currY-s).SetGreen(1);
        Pixel(currX-offset,currY-s).SetBlue(0);
      }
  }

        return values;

      }

vector<Feature> featuresVec;
vector<Feature> matchingFeatures;

void R2Image::
Harris(double sigma)
{
    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges
  R2Image Ixx(*this);
  Ixx.SobelX();
  R2Image Iyy(*this);
  Iyy.SobelY();
  R2Image Ixy(*this);

  //R2Pixel corners[width*height];

  //Ix * Iy
  for(int i = 0; i < width; i++){
    for(int j = 0; j < height; j++){
        Ixy.Pixel(i,j) = Ixx.Pixel(i,j)*Iyy.Pixel(i,j);
    }
  }

  //square to get Ixx
  for(int i = 0; i < width; i++){
    for(int j = 0; j < height; j++){
        Ixx.Pixel(i,j) *= Ixx.Pixel(i,j);
    }
  }

  //square to get Iyy
  for(int i = 0; i < width; i++){
    for(int j = 0; j < height; j++){
        Iyy.Pixel(i,j) *= Iyy.Pixel(i,j);
    }
  }

  Ixx.Blur(sigma);
  Iyy.Blur(sigma);
  Ixy.Blur(sigma);

  float red;
  float green;
  float blue;
  vector<Feature> values;

  R2Image tempImg(*this);
  //R2Pixel* tempPixel = new R2Pixel();

  for(int x = 0; x < width; x++){
    for(int y = 0; y < height; y++){
      tempImg.Pixel(x,y) = (Ixx.Pixel(x,y)*Iyy.Pixel(x,y))
      -(Ixy.Pixel(x,y)*Ixy.Pixel(x,y))
      -0.04*((Ixx.Pixel(x,y)+Iyy.Pixel(x,y))
      *(Ixx.Pixel(x,y)+Iyy.Pixel(x,y)));

      red = tempImg.Pixel(x,y).Red();
      green = tempImg.Pixel(x,y).Green();
      blue = tempImg.Pixel(x,y).Blue();


      tempImg.Pixel(x,y).SetRed(red + 0.5);
      tempImg.Pixel(x,y).SetBlue(blue + 0.5);
      tempImg.Pixel(x,y).SetGreen(green + 0.5);

      tempImg.Pixel(x,y).Clamp();
       //fprintf(stderr, "R: (%f) G:(%f) B:(%f) sum:(%f)\n", Pixel(x,y).Red(), Pixel(x,y).Blue(), Pixel(x,y).Green(), rgbsum);

          //fprintf(stderr, "Sum (%f), count(%d)\n", rgbsum, count);
          //fprintf(stderr, "R: (%f) G:(%f) B:(%f) Luminance (%f) count(%d) sum(%f)\n", Pixel(x,y).Red(), Pixel(x,y).Blue(), Pixel(x,y).Green(), Pixel(x,y).Luminance(), count, rgbsum);

          Feature temp(x, y, tempImg.Pixel(x,y));
          values.push_back(temp);

    }
  }

    sort(values.rbegin(), values.rend());

    int count = 0;
    int index = 0;
    int window = 10;
    bool empty;
    int currX;
    int currY;
  //  R2Pixel *redPixel(1, 0, 0, 0);

    while(count < 150){
      empty = true;
      currX = values[index].centerX;
      currY = values[index].centerY;

      for(int i = -window; i <= window; i++){
        for(int j = -window; j <= window; j++){
          //not going outside pic
            if(currX + i > 0 && currX + i < width &&
              currY + j > 0 && currY + j < height ){
               if(tempImg.Pixel(currX + i, currY + j).Red() == 1.0){
                 empty = false;
               }
              }
        }
      }

      if(empty){
        featuresVec.push_back(values[index]);
        tempImg.Pixel(currX, currY).SetRed(1.0);
        count++;
      }

      index++;
    }


    filterHarris(featuresVec);
    fprintf(stderr, "height (%d) width (%d) pixels (%d)\n", height, width, height*width);
    //(*this).features = filterHarris(brightness);

    // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
    fprintf(stderr, "Harris(%g) implemented\n", sigma);
}




void R2Image::
Sharpen()
{
  // Sharpen an image using a linear filter. Use a kernel of your choosing.
  float kernel[5][5] = {{-0.03,  -0.03, -0.03, -0.03, -0.03},
                        {-0.03, -0.005, -0.005, -0.005, -0.03},
                        {-0.03, -0.005, 1.6   , -0.005, -0.03},
                        {-0.03, -0.005, -0.005, -0.005, -0.03},
                        {-0.03, -0.03, -0.03, -0.03, -0.03}
                      };

  R2Image tempImg(*this);

  R2Pixel* tempPixel = new R2Pixel();

  for (int i = 0; i < width-2; i++) {
    for (int j = 0;  j < height-2; j++) {

      //change tempImg
      *tempPixel = kernel[0][0]*tempImg.Pixel(i-2, j-2) +
                    kernel[0][1]*tempImg.Pixel(i-1, j-2) +
                    kernel[0][2]*tempImg.Pixel(i, j-2) +
                    kernel[0][3]*tempImg.Pixel(i+1, j-2) +
                    kernel[0][4]*tempImg.Pixel(i+2, j-2) + ///
                    kernel[1][0]*tempImg.Pixel(i-2, j-1) +
                    kernel[1][1]*tempImg.Pixel(i-1, j-1) +
                    kernel[1][2]*tempImg.Pixel(i, j-1) +
                    kernel[1][3]*tempImg.Pixel(i+1, j-1) +
                    kernel[1][4]*tempImg.Pixel(i+2, j-1) + ///
                    kernel[2][0]*tempImg.Pixel(i-2, j) +
                    kernel[2][1]*tempImg.Pixel(i-1, j) +
                    kernel[2][2]*tempImg.Pixel(i, j) +
                    kernel[2][3]*tempImg.Pixel(i+1, j) +
                    kernel[2][4]*tempImg.Pixel(i+2, j) + ///
                    kernel[3][0]*tempImg.Pixel(i-2, j+1) +
                    kernel[3][1]*tempImg.Pixel(i-1, j+1) +
                    kernel[3][2]*tempImg.Pixel(i, j+1) +
                    kernel[3][3]*tempImg.Pixel(i+1, j+1) +
                    kernel[3][4]*tempImg.Pixel(i+2, j+1) + ///
                    kernel[4][0]*tempImg.Pixel(i-2, j+2) +
                    kernel[4][1]*tempImg.Pixel(i-1, j+2) +
                    kernel[4][2]*tempImg.Pixel(i, j+2) +
                    kernel[4][3]*tempImg.Pixel(i+1, j+2) +
                    kernel[4][4]*tempImg.Pixel(i+2, j+2);

      //temp image
      Pixel(i,j) = *tempPixel;
    }
  }

}

void R2Image::
line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
	if(x0>x1)
	{
		int x=y1;
		y1=y0;
		y0=x;

		x=x1;
		x1=x0;
		x0=x;
	}
     int deltax = x1 - x0;
     int deltay = y1 - y0;
     float error = 0;
     float deltaerr = 0.0;
	 if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
           // note that this division needs to be done in a way that preserves the fractional part
     int y = y0;
     for(int x=x0;x<=x1;x++)
	 {
		 Pixel(x,y).Reset(r,g,b,1.0);
         error = error + deltaerr;
         if(error>=0.5)
		 {
			 if(deltay>0) y = y + 1;
			 else y = y - 1;

             error = error - 1.0;
		 }
	 }
	 if(x0>3 && x0<width-3 && y0>3 && y0<height-3)
	 {
		 for(int x=x0-3;x<=x0+3;x++)
		 {
			 for(int y=y0-3;y<=y0+3;y++)
			 {
				 Pixel(x,y).Reset(r,g,b,1.0);
			 }
		 }
	 }
}


void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage"

  //search through values in a 6x6 window
  int window = 10;
  // double searchWidth = .2*width;
  // double searchHeight = .2*height;

  R2Image firstImage(*this);
  firstImage.Harris(2.0);

  R2Image secondImage(*otherImage);
  R2Image thirdImage(*this);
  // thirdImage.ChangeSaturation(0.0);
  // secondImage.ChangeSaturation(0.0);
  //Harris(2.0);
  int currX;
  int currY;
  // R2Pixel greenPix(0, 1, 0, 1);
  // R2Pixel redPix(1, 0, 0, 1);
  //
  // vector<double> ssdResults;

  double diff;
  double ssd;

  double finalSSD;
  int finalX;
  int finalY;


  double sum1;
  double sum2;



  //go though 150 features & compute avg of window
   for(int index = 0; index < featuresVec.size(); index++ ){

    currX = featuresVec[index].centerX;
    currY = featuresVec[index].centerY;

    // fprintf(stderr, "X: %d, Y: %d \n", currX, currY);

    finalX = currX;
    finalY = currY;
    finalSSD = 5000000000;

    //go through %20 of image size (search area)
    for(int x = (int) currX-(width*0.2); x <= (int) currX+(width*0.2); x++) {
        for(int y = (int) currY-(height*0.2); y <= (int) currY+(height*0.2); y++){

          //fprintf(stderr, "X: %d, Y: %d \n", dX, dY);

          ssd = 0.0;
         if(x >= 0 && x < width && y >=0 && y < height){
           //fprintf(stderr, "%s\n", inBounds ? "true" : "false" );

           //go through window to compare ssd from orig and other image
           for(int i = -window; i <= window; i++){
             for(int j = -window; j <= window ; j++){

                sum1 = thirdImage.Pixel(currX+i, currY+j).Red()+thirdImage.Pixel(currX+i, currY+j).Green() + thirdImage.Pixel(currX+i, currY+j).Blue();

                 sum2 = secondImage.Pixel(x+i, y+j).Red() +secondImage.Pixel(x+i, y+j).Green() + secondImage.Pixel(x+i, y+j).Blue();

                 //compute ssd
                 diff = sum1 - sum2;
                 ssd += (diff*diff);
                 //fprintf(stderr, "SSD: %f\n", ssd);

             }
           }

           //fprintf(stderr, "CURR SSD: %f , FINAL SSD: %f \n", ssd, finalSSD);
           //update coordinates if distance is smaller
           if(ssd <= finalSSD){
             finalSSD = ssd;
             finalX = x;
             finalY = y;

           }
         }

        }
    }



    //fprintf(stderr, "prevX: %d prevY: %d finalX: %d finalY: %d \n", currX, currY, finalX, finalY);
    Feature tempFeature(finalX, finalY, secondImage.Pixel(finalX,finalY));
    matchingFeatures.push_back(tempFeature);
    //create line
    line(currX, finalX, currY, finalY, 0, 1, 0);
  }

/*
 Reject outliers
*/

  int count = 0;
  //number of iterations needed
  int iterations = 150;
  int maxInliers = 0;
  int currentInliers = 0;

  double threshhold = 3;

  //keep track of vectors to reject
  vector<int> badFeatures;
  vector<int> finalBadFeatures;

  Feature temp1;
  Feature temp2;


  int index1;


  while(count < iterations){

    index1 = rand()%150;
    temp1 = featuresVec[index1];
    temp2 = matchingFeatures[index1];

    int xOrig = temp1.centerX;
    int yOrig = temp1.centerY;

    int xMatch = temp2.centerX;
    int yMatch = temp2.centerY;

    //vector
    int xVec = xOrig - xMatch;
    int yVec = yOrig - yMatch;

    //compute vectors for other pairs
    for(int i = 0; i < 150; i++){

      int x1 = featuresVec[i].centerX;
      int y1 = featuresVec[i].centerY;

      int x2 = matchingFeatures[i].centerX;
      int y2 = matchingFeatures[i].centerY;

      //vector
      int xVec2 = x1 - x2;
      int yVec2 = y1 - y2;

      //difference vector coordinates
      int xDiff = xVec - xVec2;
      int yDiff = yVec - yVec2;

      //magnitude of difference vector
      double difference = sqrt(pow(xDiff,2) + pow(yDiff,2));
      fprintf(stderr, "Difference: %f \n", difference);

      if( difference < threshhold){
        currentInliers++;
      }else{
        badFeatures.push_back(i);
      }

    }

    //update maxInliers
    if(maxInliers < currentInliers){
      maxInliers = currentInliers;
      finalBadFeatures = badFeatures;
      fprintf(stderr, "Max Inliers: %d \n", maxInliers);

    }else{
      badFeatures.clear();
    }

    currentInliers = 0;
    count++;
  }

  //color rejecting lines
  for(int k = 0; k < finalBadFeatures.size(); k++){
    int startX = featuresVec[finalBadFeatures[k]].centerX;
    int startY = featuresVec[finalBadFeatures[k]].centerY;
    int stopX = matchingFeatures[finalBadFeatures[k]].centerX;
    int stopY = matchingFeatures[finalBadFeatures[k]].centerY;

    line(startX, stopX, startY, stopY, 1, 0, 0);
  }

	return;
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.

  //hw8 4 points in A image



	fprintf(stderr, "fit other image using a homography and blend imageB over imageA\n");
	return;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp);
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);

  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);

  // Check info header
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }

  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }

  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }

    // Close file
    fclose(fp);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" {
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 100, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}
