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


vector<double> R2Image::
svdTest(vector<R2Point> startPoints, vector<R2Point> endPoints)
//change to take in 4 points....or better, 2 vectors of points for features from image A and image B
{


  int N = startPoints.size();
  double** linEquations = dmatrix(1,2*N,1,9);

  for(int i = 0; i < N; i++){
      //first row
      linEquations[2*i +1][1] = 0.0;
    	linEquations[2*i +1][2] = 0.0;
    	linEquations[2*i +1][3] = 0.0;
    	linEquations[2*i +1][4] = (-1.0)*startPoints[i][0]; //-w*x_a
    	linEquations[2*i +1][5] = (-1.0)*startPoints[i][1]; //-w*y_a
    	linEquations[2*i +1][6] = -1.0; //-w*w
      linEquations[2*i +1][7] = endPoints[i][1]*startPoints[i][0]; //y_b*x_a
      linEquations[2*i +1][8] = endPoints[i][1]*startPoints[i][1]; //y_b*y_a
      linEquations[2*i +1][9] = endPoints[i][1]; //y_b*w

      //second row
      linEquations[2*i +2][1] = startPoints[i][0]; //w*x_a
      linEquations[2*i +2][2] = startPoints[i][1]; //w*y_a
      linEquations[2*i +2][3] = 1.0; //w*w
      linEquations[2*i +2][4] = 0.0;
      linEquations[2*i +2][5] = 0.0;
      linEquations[2*i +2][6] = 0.0;
      linEquations[2*i +2][7] = -(endPoints[i][0]*startPoints[i][0]); //-x_b*x_a
      linEquations[2*i +2][8] = -(endPoints[i][0]*startPoints[i][1]); //-x_b*y_a
      linEquations[2*i +2][9] = (-1.0)*endPoints[i][0]; //-x_b*w

  }

  // fprintf(stderr, "Matrix for A points: \n");
  // for(int k = 1; k < 9; k++){
  //   fprintf(stderr, "%f %f %f %f %f %f %f %f %f\n", linEquations[k][1], linEquations[k][2], linEquations[k][3], linEquations[k][4], linEquations[k][5], linEquations[k][6], linEquations[k][7], linEquations[k][8], linEquations[k][9]);
  //
  // }


	// printf("\n Four points: \n");
	// printf("Point #1: %f,%f\n",startPoints[0][0],startPoints[0][1]);
	// printf("Point #2: %f,%f\n",startPoints[1][0],startPoints[1][1]);
	// printf("Point #3: %f,%f\n",startPoints[2][0],startPoints[2][1]);
	// printf("Point #4: %f,%f\n",startPoints[3][0],startPoints[3][1]);
  //
  // printf("\n Other points: \n");
	// printf("Point #1: %f,%f\n",endPoints[0][0],endPoints[0][1]);
	// printf("Point #2: %f,%f\n",endPoints[1][0],endPoints[1][1]);
	// printf("Point #3: %f,%f\n",endPoints[2][0],endPoints[2][1]);
	// printf("Point #4: %f,%f\n",endPoints[3][0],endPoints[3][1]);

	// compute the SVD
  double singularValues[10]; // 1..9

  double** nullspaceMatrix = dmatrix(1,9,1,9);
	svdcmp(linEquations, 2*N, 9, singularValues, nullspaceMatrix);

	// get the result
	// fprintf(stderr, "\nSingular values A: %f, %f, %f, %f, %f, %f, %f, %f, %f\n", singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6], singularValues[7], singularValues[8], singularValues[9]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

  //fprintf(stderr, "smallestIndexA %d\n", smallestIndex);

  // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  // fprintf(stderr, "H Values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex],nullspaceMatrix[7][smallestIndex],nullspaceMatrix[8][smallestIndex],nullspaceMatrix[9][smallestIndex]);


  // fprintf(stderr, "H Matrix for A points: \n");
  // for(int j = 1; j < 9; j = j+3){
  //     fprintf(stderr, "%f, %f, %f \n", nullspaceMatrix[j][smallestIndex],nullspaceMatrix[j
  //     +1][smallestIndex], nullspaceMatrix[j+2][smallestIndex]);
  // }


	// make sure the solution is correct:
  //   for(int n = 1; n < 9; n ++){
  //     printf("Equation #%d result: %f\n", n,	linEquationsA[n][1]*nullspaceMatrixA[1][smallestIndex] +
  //   	linEquationsA[n][2]*nullspaceMatrixA[2][smallestIndex] +
  //   	linEquationsA[n][3]*nullspaceMatrixA[3][smallestIndex] +
  //   	linEquationsA[n][4]*nullspaceMatrixA[4][smallestIndex] +
  //   	linEquationsA[n][5]*nullspaceMatrixA[5][smallestIndex] +
  //   	linEquationsA[n][6]*nullspaceMatrixA[6][smallestIndex] +
  //     linEquationsA[n][7]*nullspaceMatrixA[7][smallestIndex] +
  //     linEquationsA[n][8]*nullspaceMatrixA[8][smallestIndex] +
  //     linEquationsA[n][9]*nullspaceMatrixA[9][smallestIndex]);
  //
  //   }
  vector<double> hValues;
  for(int i = 1; i < 10; i++){
       hValues.push_back(nullspaceMatrix[i][smallestIndex]);
  }

	return hValues;
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

      // for(int s = -offset; s <= offset; s++){
      //
      //   Pixel(currX+s,currY+offset).SetRed(0);
      //   Pixel(currX+s,currY+offset).SetGreen(1);
      //   Pixel(currX+s,currY+offset).SetBlue(0);
      //
      //   Pixel(currX+offset,currY+s).SetRed(0);
      //   Pixel(currX+offset,currY+s).SetGreen(1);
      //   Pixel(currX+offset,currY+s).SetBlue(0);
      //   //
      //   Pixel(currX-s,currY-offset).SetRed(0);
      //   Pixel(currX-s,currY-offset).SetGreen(1);
      //   Pixel(currX-s,currY-offset).SetBlue(0);
      //   //
      //   Pixel(currX-offset,currY-s).SetRed(0);
      //   Pixel(currX-offset,currY-s).SetGreen(1);
      //   Pixel(currX-offset,currY-s).SetBlue(0);
      // }
  }

        return values;

      }

vector<Feature> featuresVec;

void R2Image::
Harris(double sigma)
{
    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges
  R2Image *Ixx = new R2Image(*this);
  Ixx->SobelX();
  R2Image *Iyy = new R2Image(*this);
  Iyy->SobelY();
  R2Image *Ixy = new R2Image(*this);

  //R2Pixel corners[width*height];

  //Ix * Iy
  for(int i = 0; i < width; i++){
    for(int j = 0; j < height; j++){
        Ixy->Pixel(i,j) = Ixx->Pixel(i,j)*Iyy->Pixel(i,j);
    }
  }

  //square to get Ixx
  for(int i = 0; i < width; i++){
    for(int j = 0; j < height; j++){
        Ixx->Pixel(i,j) *= Ixx->Pixel(i,j);
    }
  }

  //square to get Iyy
  for(int i = 0; i < width; i++){
    for(int j = 0; j < height; j++){
        Iyy->Pixel(i,j) *= Iyy->Pixel(i,j);
    }
  }

  Ixx->Blur(sigma);
  Iyy->Blur(sigma);
  Ixy->Blur(sigma);

  float red;
  float green;
  float blue;
  vector<Feature> values;

  R2Image *tempImg = new R2Image(*this);
  //R2Pixel* tempPixel = new R2Pixel();

  for(int x = 0; x < width; x++){
    for(int y = 0; y < height; y++){
      tempImg->Pixel(x,y) = (Ixx->Pixel(x,y)*Iyy->Pixel(x,y))
      -(Ixy->Pixel(x,y)*Ixy->Pixel(x,y))
      -0.04*((Ixx->Pixel(x,y)+Iyy->Pixel(x,y))
      *(Ixx->Pixel(x,y)+Iyy->Pixel(x,y)));

      red = tempImg->Pixel(x,y).Red();
      green = tempImg->Pixel(x,y).Green();
      blue = tempImg->Pixel(x,y).Blue();

      tempImg->Pixel(x,y).SetRed(red + 0.5);
      tempImg->Pixel(x,y).SetBlue(blue + 0.5);
      tempImg->Pixel(x,y).SetGreen(green + 0.5);

      tempImg->Pixel(x,y).Clamp();

      Feature temp(x, y, tempImg->Pixel(x,y));
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

    while(count < 250){
      empty = true;
      currX = values[index].centerX;
      currY = values[index].centerY;

      for(int i = -window; i <= window; i++){
        for(int j = -window; j <= window; j++){
          //not going outside pic
            if(currX + i > 0 && currX + i < width &&
              currY + j > 0 && currY + j < height ){
               if(tempImg->Pixel(currX + i, currY + j).Red() == 1.0){
                 empty = false;
               }
              }
        }
      }

      if(empty){
        featuresVec.push_back(values[index]);
        tempImg->Pixel(currX, currY).SetRed(1.0);
        count++;
      }

      index++;
    }


    filterHarris(featuresVec);
    fprintf(stderr, "height (%d) width (%d) pixels (%d)\n", height, width, height*width);
    //(*this).features = filterHarris(brightness);

    // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
    fprintf(stderr, "Harris(%g) implemented\n", sigma);

    delete Ixx;
    delete Iyy;
    delete Ixy;
    delete tempImg;

    return;
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

//draw a line
void R2Image::
line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
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

// 	if(x0>x1)
// 	{
// 		int x=y1;
// 		y1=y0;
// 		y0=x;
//
// 		x=x1;
// 		x1=x0;
// 		x0=x;
// 	}
//      int deltax = x1 - x0;
//      int deltay = y1 - y0;
//      float error = 0;
//      float deltaerr = 0.0;
// 	 if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
//            // note that this division needs to be done in a way that preserves the fractional part
//      int y = y0;
//
//    for(int x=x0;x<=x1;x++)
// {
//   Pixel(x,y).Reset(r,g,b,1.0);
//        error = error + deltaerr;
//        if(error>=0.5)
//   {
//     if(deltay>0) y = y + 1;
//     else y = y - 1;
//
//            error = error - 1.0;
//   }
// }

// delete Ixx;
// delete Iyy;
// delete Ixy;

}


void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage"

  //search through values in a 6x6 window
  int window = 6;
  vector<Feature> matchingFeatures;


  R2Image firstImage(*this);
  firstImage.Harris(2.0);

  R2Image secondImage(*otherImage);
  R2Image thirdImage(*this);

  int currX;
  int currY;

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

    finalX = currX;
    finalY = currY;
    finalSSD = 5000000000;

    //go through %20 of image size (search area)
    for(int x = (int) currX-(width*0.05); x <= (int) currX+(width*0.05); x++) {
        for(int y = (int) currY-(height*0.05); y <= (int) currY+(height*0.05); y++){

          //fprintf(stderr, "X: %d, Y: %d \n", dX, dY);

          ssd = 0.0;
         if(x >= 0 && x < width && y >=0 && y < height){
           //fprintf(stderr, "%s\n", inBounds ? "true" : "false" );

           //go through window to compare ssd from orig and other image
           for(int i = -window; i <= window; i++){
             for(int j = -window; j <= window ; j++){

                // tempPixel = thirdImage.Pixel(currX+i, currY+j) - secondImage.Pixel(x+i, y+j);
                // colorsum = tempPixel.Red() + tempPixel.Green() + tempPixel.Blue();

                // sum1 = thirdImage.Pixel(currX+i, currY+j).Red()+thirdImage.Pixel(currX+i, currY+j).Green() + thirdImage.Pixel(currX+i, currY+j).Blue();
                sum1 = thirdImage.Pixel(currX+i, currY+j).Luminance();
                //
                // sum2 = secondImage.Pixel(x+i, y+j).Red() +secondImage.Pixel(x+i, y+j).Green() + secondImage.Pixel(x+i, y+j).Blue();
                sum2 = secondImage.Pixel(x+i, y+j).Luminance();
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
    otherImage->line(currX, finalX, currY, finalY, 0, 1, 0);
  }

/*
 Reject outliers
*/

//   int count = 0;
//   //number of iterations needed
//   int iterations = 150;
//   int maxInliers = 0;
//   int currentInliers = 0;
//
//   double threshhold = 3;
//
//   //keep track of vectors to reject
//   vector<int> badFeatures;
//   vector<int> finalBadFeatures;
//   vector<int> goodFeatures;
//   vector<int> finalGoodFeatures;
//
//   Feature temp1;
//   Feature temp2;
//
//
//   int index1;
//
//
//   while(count < iterations){
//
//     index1 = rand()%150;
//     temp1 = featuresVec[index1];
//     temp2 = matchingFeatures[index1];
//
//     int xOrig = temp1.centerX;
//     int yOrig = temp1.centerY;
//
//     int xMatch = temp2.centerX;
//     int yMatch = temp2.centerY;
//
//     //vector
//     int xVec = xOrig - xMatch;
//     int yVec = yOrig - yMatch;
//
//     //compute vectors for other pairs
//     for(int i = 0; i < 150; i++){
//
//       int x1 = featuresVec[i].centerX;
//       int y1 = featuresVec[i].centerY;
//
//       int x2 = matchingFeatures[i].centerX;
//       int y2 = matchingFeatures[i].centerY;
//
//       //vector
//       int xVec2 = x1 - x2;
//       int yVec2 = y1 - y2;
//
//       //difference vector coordinates
//       int xDiff = xVec - xVec2;
//       int yDiff = yVec - yVec2;
//
//       //magnitude of difference vector
//       double difference = sqrt(pow(xDiff,2) + pow(yDiff,2));
//       //fprintf(stderr, "Difference: %f \n", difference);
//
//       if( difference < threshhold){
//         currentInliers++;
//         goodFeatures.push_back(i);
//       }else{
//         badFeatures.push_back(i);
//       }
//
//     }
//
// //     H Matrix for A points:
// // -0.004756, 0.001282, -0.808207
// // -0.001169, -0.004307, 0.588845
// // -0.000001, 0.000000, -0.004356
//
//     //update maxInliers
//     if(maxInliers < currentInliers){
//       maxInliers = currentInliers;
//       finalGoodFeatures = goodFeatures;
//       finalBadFeatures = badFeatures;
//       fprintf(stderr, "Max Inliers: %d \n", maxInliers);
//
//     }else{
//       badFeatures.clear();
//     }
//
//     currentInliers = 0;
//     count++;
//   }
//
//   //color rejecting lines
//   for(int k = 0; k < finalBadFeatures.size(); k++){
//     int startX = featuresVec[finalBadFeatures[k]].centerX;
//     int startY = featuresVec[finalBadFeatures[k]].centerY;
//     int stopX = matchingFeatures[finalBadFeatures[k]].centerX;
//     int stopY = matchingFeatures[finalBadFeatures[k]].centerY;
//
//     otherImage->line(startX, stopX, startY, stopY, 1, 0, 0);
//   }


	return;
}

vector<Feature> matchingFeatures;
vector<Feature> prevImgFeatures;

vector<Feature> R2Image::
findMatchingFeatures(R2Image* otherImage){

  matchingFeatures.clear();

  int window = 5;

  R2Image secondImage(*otherImage);
  R2Image thirdImage(*this);

  int currX;
  int currY;

  double diff;
  double ssd;

  double finalSSD;
  int finalX;
  int finalY;


  double sum1;
  double sum2;



  vector<Feature> featuresToTrack;
  if(prevImgFeatures.size() != 0){
    featuresToTrack = prevImgFeatures;
    fprintf(stderr, "FEATURES TO TRACK: %lu  PREV FEATURES: %lu\n", featuresToTrack.size(), prevImgFeatures.size());
  }else{
    featuresToTrack = featuresVec;
  }



  //go though 150 features & compute avg of window
   for(int index = 0; index < featuresToTrack.size(); index++ ){

    currX = featuresToTrack[index].centerX;
    currY = featuresToTrack[index].centerY;

    // fprintf(stderr, "X: %d, Y: %d \n", currX, currY);

    finalX = currX;
    finalY = currY;
    finalSSD = 5000000000;

    int swidth = 0.08*width;
    int sheight = 0.08*height;


    //go through %20 of image size (search area)
    for(int x = (int) currX-(swidth); x <= (int) currX+(swidth); x++) {
        for(int y = (int) currY-(sheight); y <= (int) currY+(sheight); y++){
          //fprintf(stderr, "WINDOW (%d, %d)\n", x, y );
         ssd = 0.0;
         //otherSum = 0.0
         //R2Pixel *otherPixel = new R2Pixel();
         if((x - window) >= 0 && (x+ window) < width && (y - window) >=0 && (y+window) < height){
           //R2Pixel *otherPixel = new R2Pixel();
           //go through window to compare ssd from orig and other image
           for(int i = -window; i <= window; i++){
             for(int j = -window; j <= window ; j++){

               sum1 = thirdImage.Pixel(currX+i, currY+j).Red()+thirdImage.Pixel(currX+i, currY+j).Green() + thirdImage.Pixel(currX+i, currY+j).Blue();
               //sum1 = thirdImage.Pixel(currX+i, currY+j).Luminance();
                //sum2 = secondImage.Pixel(x+i, y+j).Luminance();
               sum2 = secondImage.Pixel(x+i, y+j).Red() +secondImage.Pixel(x+i, y+j).Green() + secondImage.Pixel(x+i, y+j).Blue();
                //  *otherPixel += (secondImage.Pixel(x+i, y+j) - thirdImage.Pixel(currX+1, currY+j)) * (secondImage.Pixel(x+i, y+j) - thirdImage.Pixel(currX+1, currY+j));
                //  //compute ssd
                // otherSum = otherPixel->Red() + otherPixel->Blue() + otherPixel->Green();

                 //
                 diff = sum2 - sum1;
                 ssd += (diff*diff);


             }
           }

           //otherSum = otherPixel->Red() + otherPixel->Blue() + otherPixel->Green();

           if(ssd < finalSSD){
             finalSSD = ssd;
             finalX = x;
             finalY = y;
             //fprintf(stderr, "%f\n", finalSSD );
           }
         }

        }
    }

    Feature tempFeature(finalX, finalY, secondImage.Pixel(finalX,finalY));
    matchingFeatures.push_back(tempFeature);
    //create line
    otherImage->line(currX, finalX, currY, finalY, 0, 1, 0);
    //fprintf(stderr, "Feature #%d\n", index );
  }

    return matchingFeatures;
}


vector<Feature> R2Image::
blendOtherImageHomography(R2Image *otherImage)
{

  int count = 0;
  double threshhold = 7.0;

  vector<int> badFeatures;
  vector<int> finalBadFeatures;
  vector<R2Point> goodFeaturesStart;
  vector<R2Point> goodFeaturesEnd;
  vector<R2Point> finalGoodFeaturesStart;
  vector<R2Point> finalGoodFeaturesEnd;
  vector<int> goodFeatures;
  vector<int> finalGoodFeatures;
  vector<Feature> matchingFeatures;


  //R2Image firstImage(*this);
  // firstImage.Harris(2.0);
  matchingFeatures = findMatchingFeatures(otherImage);
  fprintf(stderr, "Matching Features Size: %lu\n",  matchingFeatures.size());

  R2Image secondImage(*otherImage);
  R2Image thirdImage(*this);

  int randIndex1;
  int randIndex2;
  int randIndex3;
  int randIndex4;

  Feature r1;
  Feature r2;
  Feature r3;
  Feature r4;

  Feature m1;
  Feature m2;
  Feature m3;
  Feature m4;

  vector<R2Point> startPoints;
  vector<R2Point> endPoints;
  vector<double> bestHMatrix;

  double computation = 0.0;
  //double computation2 = 0.0;
  int inliers = 0;
  int maxInliers = 0;

  int size = featuresVec.size();
  fprintf(stderr, "FeaturesVec size: %lu\n", featuresVec.size());

  while(count < 1000){

    //fprintf(stderr, "Loop %d\n", count);
    //use to get 4 random points
    randIndex1 = rand()%size;
    randIndex2 = rand()%size;
    randIndex3 = rand()%size;
    randIndex4 = rand()%size;

    r1 = featuresVec[randIndex1];
    r2 = featuresVec[randIndex2];
    r3 = featuresVec[randIndex3];
    r4 = featuresVec[randIndex4];

    m1 = matchingFeatures[randIndex1];
    m2 = matchingFeatures[randIndex2];
    m3 = matchingFeatures[randIndex3];
    m4 = matchingFeatures[randIndex4];

    R2Point p1(r1.centerX,r1.centerY);
    R2Point p2(r2.centerX,r2.centerY);
    R2Point p3(r3.centerX,r3.centerY);
    R2Point p4(r4.centerX,r4.centerY);

    R2Point p5(m1.centerX,m1.centerY);
    R2Point p6(m2.centerX,m2.centerY);
    R2Point p7(m3.centerX,m3.centerY);
    R2Point p8(m4.centerX,m4.centerY);

    startPoints.push_back(p1);
    startPoints.push_back(p2);
    startPoints.push_back(p3);
    startPoints.push_back(p4);

    endPoints.push_back(p5);
    endPoints.push_back(p6);
    endPoints.push_back(p7);
    endPoints.push_back(p8);

    vector<double> hValues = svdTest(startPoints, endPoints);

    for(int i = 0; i < matchingFeatures.size(); i++){
      vector<double> pointA;
      pointA.push_back(featuresVec[i].centerX); //x_a
      pointA.push_back(featuresVec[i].centerY); //y_a
      pointA.push_back(1); //z_a

      vector<double> pointB;
      pointB.push_back(matchingFeatures[i].centerX); //x_b
      pointB.push_back(matchingFeatures[i].centerY); //y_b

      double xVal = 0.0;
      double yVal = 0.0;
      double zVal = 0.0;

      xVal = hValues[0]*pointA[0] + hValues[1]*pointA[1] + hValues[2]*pointA[2];
      yVal = hValues[3]*pointA[0] + hValues[4]*pointA[1] + hValues[5]*pointA[2];
      zVal = hValues[6]*pointA[0] + hValues[7]*pointA[1] + hValues[8]*pointA[2];

      // fprintf(stderr, "A Point: <%f, %f>\n", pointA[0], pointA[1] );
      // fprintf(stderr, "Vector:  <%f, %f, %f>\n", xVal, yVal, zVal);

      double xResult = xVal/zVal;
      double yResult = yVal/zVal;

      // fprintf(stderr, "Result Point: <%f, %f>\n", xResult, yResult);
      // fprintf(stderr, "B Point: <%f, %f>\n", pointB[0], pointB[1]);

      computation = sqrt(pow((xResult - pointB[0]),2) + pow(yResult - pointB[1],2));

      // fprintf(stderr, "\nComputation of Distance: %f\n", computation );
      if(computation < threshhold){
        //fprintf(stderr, "%s\n", );
        inliers++;
        goodFeatures.push_back(i);
        goodFeaturesStart.push_back(R2Point(pointA[0], pointA[1]));
        goodFeaturesEnd.push_back(R2Point(pointB[0], pointB[1]));
      } else {
        badFeatures.push_back(i);
      }

      computation = 0.0;
      pointA.clear();
      pointB.clear();
    }

    if(inliers > maxInliers){
      maxInliers = inliers;
      finalGoodFeatures = goodFeatures;
      finalGoodFeaturesStart = goodFeaturesStart;
      finalGoodFeaturesEnd = goodFeaturesEnd;
      finalBadFeatures = badFeatures;
      bestHMatrix = hValues;
      fprintf(stderr, "\nMaxInliers: %d\n", maxInliers );

      fprintf(stderr, "\n\nBest H Matrix:\n");
      fprintf(stderr, "%f %f %f\n", bestHMatrix[0], bestHMatrix[1], bestHMatrix[2]);
      fprintf(stderr, "%f %f %f\n", bestHMatrix[3], bestHMatrix[4], bestHMatrix[5]);
      fprintf(stderr, "%f %f %f\n\n", bestHMatrix[6], bestHMatrix[7], bestHMatrix[8]);


      fprintf(stderr, "H Matrix Normalized\n");
      for(int k = 0; k < 8; k++){
        fprintf(stderr, "%f ", bestHMatrix[k]/bestHMatrix[8]);
      }
      fprintf(stderr, "%f ", bestHMatrix[8]);

    }




    inliers = 0;
    count++;
    startPoints.clear();
    endPoints.clear();
    goodFeatures.clear();
    goodFeaturesStart.clear();
    goodFeaturesEnd.clear();
    badFeatures.clear();
  }

  //fprintf(stderr, "final maxInliers: %d Vector size: %lu %lu\n", maxInliers, finalGoodFeaturesStart.size(), finalGoodFeaturesEnd.size());

  vector<double> best = svdTest(finalGoodFeaturesStart, finalGoodFeaturesEnd);

  fprintf(stderr, "\n\nInverse H Matrix:\n");
  fprintf(stderr, "%f %f %f\n", best[0], best[1], best[2]);
  fprintf(stderr, "%f %f %f\n", best[3], best[4], best[5]);
  fprintf(stderr, "%f %f %f\n", best[6], best[7], best[8]);

  fprintf(stderr, "\n\nFinal Best H Matrix:\n");
  fprintf(stderr, "%f %f %f\n", bestHMatrix[0], bestHMatrix[1], bestHMatrix[2]);
  fprintf(stderr, "%f %f %f\n", bestHMatrix[3], bestHMatrix[4], bestHMatrix[5]);
  fprintf(stderr, "%f %f %f\n", bestHMatrix[6], bestHMatrix[7], bestHMatrix[8]);


  fprintf(stderr, "finalGoodFeatures size: %lu\n", finalGoodFeatures.size());
  fprintf(stderr, "finalBadFeatures siez: %lu\n", finalBadFeatures.size());

  // vector<Feature> finalFeatures;

  for(int k = 0; k < finalGoodFeatures.size(); k++){
    fprintf(stderr, "Features %d: (%d, %d)\n", k,  featuresVec[finalGoodFeatures[k]].centerX,
    featuresVec[finalGoodFeatures[k]].centerY);
    int startX = featuresVec[finalGoodFeatures[k]].centerX;
    int startY = featuresVec[finalGoodFeatures[k]].centerY;
    int stopX = matchingFeatures[finalGoodFeatures[k]].centerX;
    int stopY = matchingFeatures[finalGoodFeatures[k]].centerY;

    //finalFeatures.push_back(matchingFeatures[finalGoodFeatures[k]]);

    otherImage->line(stopX, startX, stopY, startY, 0, 1, 0);
 }

/****
BEGIN
  // for(int j = 0; j < finalBadFeatures.size(); j++){
  //   int startX = featuresVec[finalBadFeatures[j]].centerX;
  //   int startY = featuresVec[finalBadFeatures[j]].centerY;
  //   int stopX = matchingFeatures[finalBadFeatures[j]].centerX;
  //   int stopY = matchingFeatures[finalBadFeatures[j]].centerY;
  //
  //   otherImage->line(stopX, startX, stopY, startY, 1, 0, 0);
  // }

  // Best H Matrix:
  // -0.005324 0.001369 -0.845843
  // -0.001127 -0.004889 0.533353
  // -0.000000 0.000000 -0.005295

  // 1 0 x
  // 0 1 y
  // 0 0 1




    R2Image *warpFrom = new R2Image(*otherImage);

   // double** bestMatrix = dmatrix(1,3,1,3);
   //  bestMatrix[1][1] = bestHMatrix[0];
   //  bestMatrix[1][2] = bestHMatrix[1];
   //  bestMatrix[1][3] = bestHMatrix[2];
   //
   //  bestMatrix[2][1] = bestHMatrix[3];
   //  bestMatrix[2][2] = bestHMatrix[4];
   //  bestMatrix[2][3] = bestHMatrix[5];
   //
   //  bestMatrix[3][1] = bestHMatrix[6];
   //  bestMatrix[3][2] = bestHMatrix[7];
   //  bestMatrix[3][3] = bestHMatrix[8];


    // for(int i = 1; i < 4; i++){
    //   bestMatrix[1][i] = bestHMatrix[i];
    //   bestMatrix[2][i] = bestHMatrix[i+3];
    //   bestMatrix[3][i] = bestHMatrix[i+6];
    // }


    // fprintf(stderr, "\nH Matrix:\n");
    // fprintf(stderr, "%f %f %f\n", bestMatrix[1][1], bestMatrix[1][2], bestMatrix[1][3]);
    // fprintf(stderr, "%f %f %f\n", bestMatrix[2][1], bestMatrix[2][2], bestMatrix[2][3]);
    // fprintf(stderr, "%f %f %f\n", bestMatrix[3][1], bestMatrix[3][2], bestMatrix[3][3]);
    //

    //
    // float det =
    // bestMatrix[1][1] * (bestMatrix[2][2]*bestMatrix[3][3] - bestMatrix[3][2] * bestMatrix[2][3])
    // - bestMatrix[1][2] * (bestMatrix[2][1]*bestMatrix[3][3] - bestMatrix[2][3]*bestMatrix[3][1]
    // + bestMatrix[1][3]*(bestMatrix[2][1]*bestMatrix[3][2] - bestMatrix[2][2]*bestMatrix[3][1]));


    double** invMatrix = dmatrix(1,3,1,3);
     invMatrix[1][1] = best[0];
     invMatrix[1][2] = best[1];
     invMatrix[1][3] = best[2];

     invMatrix[2][1] = best[3];
     invMatrix[2][2] = best[4];
     invMatrix[2][3] = best[5];

     invMatrix[3][1] = best[6];
     invMatrix[3][2] = best[7];
     invMatrix[3][3] = best[8];


    // invMatrix = bestMatrix;
    // invMatrix[1][1] = (bestMatrix[2][2]*bestMatrix[3][3] - bestMatrix[3][2]*bestMatrix[2][3])/det;
    //
    // invMatrix[1][2] = (bestMatrix[1][3]*bestMatrix[3][2] - bestMatrix[1][2]*bestMatrix[3][3])/det;
    //
    // invMatrix[1][3] = (bestMatrix[1][2]*bestMatrix[2][3] - bestMatrix[1][3]*bestMatrix[2][2])/det;
    //
    // invMatrix[2][1] = (bestMatrix[2][3]*bestMatrix[3][1] - bestMatrix[2][1]*bestMatrix[3][3])/det;
    //
    // invMatrix[2][2] = (bestMatrix[1][1]*bestMatrix[3][3] - bestMatrix[1][3]*bestMatrix[3][1])/det;
    //
    // invMatrix[2][3] = (bestMatrix[2][1]*bestMatrix[1][3] - bestMatrix[1][1]*bestMatrix[2][3])/det;
    //
    // invMatrix[3][1] = (bestMatrix[2][1]*bestMatrix[3][2] - bestMatrix[3][1]*bestMatrix[2][2])/det;
    //
    // invMatrix[3][2] = (bestMatrix[3][1]*bestMatrix[1][2] - bestMatrix[1][1]*bestMatrix[3][2])/det;
    //
    // invMatrix[3][3] = (bestMatrix[1][1]*bestMatrix[2][2] - bestMatrix[2][1]*bestMatrix[1][2])/det;

    // fprintf(stderr, "\nH Inverses:\n");
    // fprintf(stderr, "%f %f %f\n", invMatrix[1][1], invMatrix[1][2], invMatrix[1][3]);
    // fprintf(stderr, "%f %f %f\n", invMatrix[2][1], invMatrix[2][2], invMatrix[2][3]);
    // fprintf(stderr, "%f %f %f\n", invMatrix[3][1], invMatrix[3][2], invMatrix[3][3]);


    float xVal2;
    float yVal2;
    float zVal2;

    float xResult2;
    float yResult2;

    for(int i = 0; i < width; i++){
      for(int j = 0; j < height; j++){
        //fprintf(stderr, "Error at: (%d, %d)\n", i, j );
        xVal2 = 0.0;
        yVal2 = 0.0;
        zVal2 = 0.0;

        xVal2 = invMatrix[1][1]*i + invMatrix[1][2]*j + invMatrix[1][3];
        yVal2 = invMatrix[2][1]*i + invMatrix[2][2]*j + invMatrix[2][3];
        zVal2 = invMatrix[3][1]*i + invMatrix[3][2]*j + invMatrix[3][3];

        xResult2 = xVal2/zVal2;
        yResult2 = yVal2/zVal2;

        //fprintf(stderr, "Result Point <%f, %f>\n", xResult2, yResult2);
        int xround = (int) xResult2;
        int yround = (int) yResult2;


        if((xround < 0 || xround >= width) || (yround < 0 || yround >= height) ){
          //fprintf(stderr, "Out of Bounds: <%d, %d>\n", xround, yround);
          otherImage->Pixel(i,j).Reset(1, 1, 1, 1);
        } else {
          //fprintf(stderr, "IN BOUNDS: <%d, %d>\n", xround, yround);


          double redAvg = (Pixel(i, j).Red() + warpFrom->Pixel(xResult2, yResult2).Red())/2;
          double blueAvg = (Pixel(i, j).Blue() + warpFrom->Pixel(xResult2, yResult2).Blue())/2;
          double greenAvg = (Pixel(i, j).Green() + warpFrom->Pixel(xResult2, yResult2).Green())/2;

          otherImage->Pixel(i,j).Reset(redAvg, greenAvg, blueAvg, 1);
        }
      }

    }


  delete warpFrom;
END OF WARP
 ***/
	fprintf(stderr, "fit other image using a homography and blend imageB over imageA\n");
	return matchingFeatures;
}

void R2Image::
FirstFrameProcessing(){
  local_firstImage = new R2Image(*this);
  local_firstImage->Harris(2.0);
  local_latestImage = new R2Image(*this);
  prevImgFeatures = findMatchingFeatures(this);

}

void R2Image::
FrameProcessing(R2Image *otherImage){
  prevImgFeatures = local_latestImage->findMatchingFeatures(otherImage);
  local_latestImage = new R2Image(*otherImage);
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
