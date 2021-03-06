// Include file for image class
#ifndef R2_IMAGE_INCLUDED
#define R2_IMAGE_INCLUDED
#include <vector>
using namespace std;


// this structure is defined at the same level as function definitions

struct svdResult
{
  double** Hmatrix;
  double* hValues;
};


struct Feature
{
    int centerX;
    int centerY;
    R2Pixel HarrisValue;

    Feature( )
    {
        centerX = -1;
        centerY = -1;
    }


    Feature(int x, int y, R2Pixel val)
    {
        centerX = x;
        centerY = y;
        HarrisValue = val;
    }

    bool operator<(const Feature& feature) const {
        double valueIntensity = HarrisValue[0] + HarrisValue[1] + HarrisValue[2];
        double featureIntensity = feature.HarrisValue[0] + feature.HarrisValue[1] + feature.HarrisValue[2];

        return valueIntensity < featureIntensity;
    }

};



typedef enum {
  R2_IMAGE_RED_CHANNEL,
  R2_IMAGE_GREEN_CHANNEL,
  R2_IMAGE_BLUE_CHANNEL,
  R2_IMAGE_ALPHA_CHANNEL,
  R2_IMAGE_NUM_CHANNELS
} R2ImageChannel;

typedef enum {
  R2_IMAGE_POINT_SAMPLING,
  R2_IMAGE_BILINEAR_SAMPLING,
  R2_IMAGE_GAUSSIAN_SAMPLING,
  R2_IMAGE_NUM_SAMPLING_METHODS
} R2ImageSamplingMethod;

typedef enum {
  R2_IMAGE_OVER_COMPOSITION,
  R2_IMAGE_IN_COMPOSITION,
  R2_IMAGE_OUT_COMPOSITION,
  R2_IMAGE_ATOP_COMPOSITION,
  R2_IMAGE_XOR_COMPOSITION,
} R2ImageCompositeOperation;


// Class definition

class R2Image {
 public:


  //std::vector<Feature> features;
  // Constructors/destructor
  R2Image(void);
  R2Image(const char *filename);
  R2Image(int width, int height);
  R2Image(int width, int height, const R2Pixel *pixels);
  R2Image(const R2Image& image);
  ~R2Image(void);

  // Image properties
  int NPixels(void) const;
  int Width(void) const;
  int Height(void) const;


  // Pixel access/update
  R2Pixel& Pixel(int x, int y);
  R2Pixel *Pixels(void);
  R2Pixel *Pixels(int row);
  R2Pixel *operator[](int row);
  const R2Pixel *operator[](int row) const;
  void SetPixel(int x, int y,  const R2Pixel& pixel);

  // Image processing
  R2Image& operator=(const R2Image& image);

  // Per-pixel operations
  void Brighten(double factor);
  void ChangeSaturation(double factor);

  // show how SVD works
  vector<double> svdTest(vector<R2Point> startPoints, vector<R2Point> endPoints);

  //helper functions
  vector<Feature> filterHarris(vector<Feature> values);
  void line(int x0, int x1, int y0, int y1, float r, float g, float b);
  vector<Feature> findMatchingFeatures(R2Image* prevImage, R2Image *currImage);
  void warp(vector<int> goodFeatures, vector<double> bestHMatrix);

  // Linear filtering operations
  void SobelX();
  void SobelY();
  void LoG();
  void Blur(double sigma);
  void Harris(double sigma);
  void Sharpen(void);

  //Extra Credit HW
  void Median(int median);
  void Bilateral();

  // further operations
  void blendOtherImageTranslated(R2Image * otherImage);
  vector<Feature> blendOtherImageHomography(R2Image *prevImage, R2Image *currentImage, R2Image *warpImage);

  //image processing
  void imageProcessing(R2Image *firstImage);

  //video processing
  void FirstFrameProcessing(R2Image *warpImage);
  void FrameProcessing(R2Image *prevImage, R2Image *currentImage, R2Image *warpImage, int i);

  void makeMask(R2Image *currentImage);

  // File reading/writing
  int Read(const char *filename);
  int ReadBMP(const char *filename);
  int ReadPPM(const char *filename);
  int ReadJPEG(const char *filename);
  int Write(const char *filename) const;
  int WriteBMP(const char *filename) const;
  int WritePPM(const char *filename, int ascii = 0) const;
  int WriteJPEG(const char *filename) const;

 private:
  // Utility functions
  void Resize(int width, int height);
  R2Pixel Sample(double u, double v,  int sampling_method);

 private:
  R2Pixel *pixels;
  int npixels;
  int width;
  int height;
  R2Image *local_firstImage;
  R2Image *local_latestImage;

};


// Inline functions

inline int R2Image::
NPixels(void) const
{
  // Return total number of pixels
  return npixels;
}



inline int R2Image::
Width(void) const
{
  // Return width
  return width;
}



inline int R2Image::
Height(void) const
{
  // Return height
  return height;
}



inline R2Pixel& R2Image::
Pixel(int x, int y)
{
  // Return pixel value at (x,y)
  // (pixels start at lower-left and go in row-major order)
  return pixels[x*height + y];
}



inline R2Pixel *R2Image::
Pixels(void)
{
  // Return pointer to pixels for whole image
  // (pixels start at lower-left and go in row-major order)
  return pixels;
}



inline R2Pixel *R2Image::
Pixels(int x)
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}



inline R2Pixel *R2Image::
operator[](int x)
{
  // Return pixels pointer for row at x
  return Pixels(x);
}



inline const R2Pixel *R2Image::
operator[](int x) const
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}



inline void R2Image::
SetPixel(int x, int y, const R2Pixel& pixel)
{
  // Set pixel
  pixels[x*height + y] = pixel;
}



#endif
