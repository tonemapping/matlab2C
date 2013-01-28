/*
 * File:   main.c
 * Author: Josmil9
 *
 * Created on January 24th 2013, 17:02
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define sigma_s 60
#define sigma_r 0.4


int getImageWidth(char* inputfile) {
    //INPUT FILE
    FILE *fp;
    fp = fopen(inputfile, "rt");

    float pixel;
    float line[80];

    //GET IMAGE DIMENSIONS
    fgets(line, 80, fp);
    sscanf(line, "%f", &pixel);
    int m = pixel;
    fgets(line, 80, fp);
    sscanf(line, "%f", &pixel);
    int n = pixel;

    fclose(fp);
    return m;
}

int getImageHeight(char* inputfile) {
    //INPUT FILE
    FILE *fp;
    fp = fopen(inputfile, "rt");

    float pixel;
    float line[80];

    //GET IMAGE DIMENSIONS
    fgets(line, 80, fp);
    sscanf(line, "%f", &pixel);
    int m = pixel;
    fgets(line, 80, fp);
    sscanf(line, "%f", &pixel);
    int n = pixel;

    fclose(fp);
    return n;
}

float** generateImageArray(char* inputfile) {

    //INITIALIZATION

    //INPUT FILE
    FILE *fp;
    fp = fopen(inputfile, "rt");

    float pixel;
    float line[80];

    //GET IMAGE DIMENSIONS
    fgets(line, 80, fp);
    sscanf(line, "%f", &pixel);
    int rows = pixel;
    fgets(line, 80, fp);
    sscanf(line, "%f", &pixel);
    int cols = pixel;
    int imgSize = rows*cols;

    //HOLD PIXEL VALUES INTO ARRAY
    float* values = calloc(rows*cols, sizeof(float));
    float** outputArray = malloc(cols*sizeof(float*));
    for (int i=0; i<rows; ++i)
    {
        outputArray[i] = values + i*cols;
    }
    if (outputArray == '/0') {
        printf("Not enough memory!\n");
        return 0;
    }

    int row = 0;
    int col = 0;
    int it = 0;
    for (it = 0; it < imgSize; it++) {
        fgets(line, 80, fp);
        sscanf(line, "%f", &pixel);
        outputArray[row][col] = pixel;
        if (col == 399) {
            col = 0;
            row++;
        } else {
            col++;
        }
    }

    fclose(fp);
    return outputArray;
}

void generateImageFile(float** channel, char* name, int h, int w) {

    //INITIALIZATION

    //OUTPUT FILE
    FILE *fpout;
    fpout = fopen(name, "w");

    //PRINT ARRAY TO A FILE
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w-1; j++) {
            fprintf(fpout, "%.9f\t", channel[i][j]);
        }
        fprintf(fpout, "\n");
    }

    fclose(fpout);
}

//ARRAY HELPER FUNCTIONS

float** createArray(int rows, int cols) {
    float* values = calloc(rows*cols, sizeof (float));
    float** array = malloc(rows * sizeof (float*));
    for (int i = 0; i < rows; i++) {
        array[i] = values + i*cols;
    }
    for(int row =0; row<rows; row++){
        for(int col=0; col<cols; col++){
            array[row][col]=0;
        }
    }
    return array;
}

void destroyArray(float** arr) {
    free(*arr);
    free(arr);
}

//DIFF DIMENSION 2 FUNCTION AND ADD ZEROES COLUMN
float** diffH(float** array, int rowN, int colN) {
    //CREATE OUTPUT ARRAY
    float** imgDiffH = createArray(rowN, colN);

    //PERFORM DIFF FUNCTION AND ABS THE RESULT
    for (int row = 0; row < rowN; row++) {
        for (int col = 0; col < colN-1; col++) {
            imgDiffH[row][col+1] = fabs(array[row][col + 1] - array[row][col]);
        }
    }

    return imgDiffH;
}

//DIFF DIMENSION 1 FUNCTION AND ABS THE RESULT
float** diffV(float** array, int rowN, int colN) {
    //CREATE OUTPUT ARRAY
    float** imgDiffV = createArray(rowN, colN);

    //PERFORM DIFF FUNCTION
    for (int row = 0; row < rowN-1; row++) {
        for (int col = 0; col < colN; col++) {
            imgDiffV[row+1][col] = fabs(array[row+1][col] - array[row][col]);
        }
    }

    return imgDiffV;
}

//ADD CHANNELS
float** addChannels(float** red, float** green, float** blue, int rowN, int colN) {
    //CREATE OUTPUT ARRAY
    float** result = createArray(rowN, colN);

    //PERFORM ADDITION
    for (int row = 0; row < rowN; row++) {
        for (int col = 0; col < colN; col++) {
            result[row][col] = red[row][col] + green[row][col] + blue[row][col];
        }
    }

    return result;
}

//TRANSPOSE
float** transpose(float** array, int rowN, int colN) {
    //CREATE OUTPUT ARRAY
    float** result = createArray(colN, rowN);

    //PERFORM ADDITION
    for (int row = 0; row < rowN; row++) {
        for (int col = 0; col < colN; col++) {
            result[col][row] = array[row][col];
        }
    }

    return result;
}

//TransformedDomainRecursiveFilter_Horizontal
void TDRF_H(float** channel, float** derivatives, float sigma, int h, int w){
    float a = expf(-sqrtf(2)/sigma);

    //V
    float** V = createArray(h,w);
    for (int row = 0; row < h; row++) {
        for (int col = 0; col < w; col++) {
                V[row][col] = powf(a,derivatives[row][col]);
        }
    }
    //LEFT->RIGHT FILTER
    for (int row = 0; row < h; row++) {
        for (int col = 1; col < w; col++) {
                channel[row][col] = channel[row][col] + V[row][col]*(channel[row][col - 1] - channel[row][col]);
        }
    }
    //RIGHT->LEFT FILTER
    for (int row = 0; row < h; row++) {
        for (int col = w-2; col >= 0; col--) {
                channel[row][col] = channel[row][col] + V[row][col+1]*(channel[row][col + 1] - channel[row][col]);
        }
    }

    destroyArray(V);
}

void filtering(float** channel, float** dIcdx, float** dHdyT, int h, int w, float num_iterations) {
    float N = num_iterations;

    for (int i = 1; i <= num_iterations; i++) {
        float sigma_H = sigma_s;
        float sigma_H_i = sigma_H * sqrtf(3) * exp2f(N - i) / sqrtf(powf(4, N) - 1);

        TDRF_H(channel, dIcdx, sigma_H_i, h, w);
        channel = transpose(channel, h, w);

        TDRF_H(channel, dHdyT, sigma_H_i, w, h);
        channel = transpose(channel, w, h);
    }
}

int main(int argc, char** argv) {

    float** testRedChannel = generateImageArray("red.txt");
    float** testGreenChannel = generateImageArray("green.txt");
    float** testBlueChannel = generateImageArray("blue.txt");

    int h = getImageWidth("red.txt");
    int w = getImageHeight("red.txt");
    printf("%i rows x %i columns\n", h, w);

    ///////////////////////////////////////////////////////////
    //////Compute the l1-norm distance of neighbor pixels./////
    ///////////////////////////////////////////////////////////

    float** dIcdxRed = diffH(testRedChannel, h,w);
    float** dIcdxGreen = diffH(testGreenChannel, h,w);
    float** dIcdxBlue = diffH(testBlueChannel, h,w);

    float** dIcdx = addChannels(dIcdxRed, dIcdxGreen, dIcdxBlue, h,w);
    destroyArray(dIcdxRed);
    destroyArray(dIcdxGreen);
    destroyArray(dIcdxBlue);

    float** dIcdyRed = diffV(testRedChannel, h,w);
    float** dIcdyGreen = diffV(testGreenChannel, h,w);
    float** dIcdyBlue = diffV(testBlueChannel, h,w);

    float** dIcdy = addChannels(dIcdyRed, dIcdyGreen, dIcdyBlue, h,w);
    destroyArray(dIcdyRed);
    destroyArray(dIcdyGreen);
    destroyArray(dIcdyBlue);

    ///////////////////////////////////////////////////////////
    ///////Compute the derivatives of the horizontal///////////
    ///////and vertical domain transforms./////////////////////
    ///////////////////////////////////////////////////////////

    float factor = sigma_s / sigma_r;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            dIcdx[i][j] = 1 + factor * dIcdx[i][j];
        }
    }

    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            dIcdy[i][j] = 1 + factor * dIcdy[i][j];
        }
    }
    ///The vertical pass is performed using a transposed image.
    float** dHdyT = transpose(dIcdy, h, w);
    destroyArray(dIcdy);

    ///////////////////////////////////////////////////////////
    ///////////////////Perform the filtering.//////////////////
    ///////////////////////////////////////////////////////////

    filtering(testRedChannel, dIcdx, dHdyT, h, w, 10);
    filtering(testGreenChannel, dIcdx, dHdyT, h, w, 10);
    filtering(testBlueChannel, dIcdx, dHdyT, h, w, 10);

    ///////////////////////////////////////////////////////////
    ///////////////////GENERATE OUTPUT IMAGE///////////////////
    ///////////////////////////////////////////////////////////

    generateImageFile(testRedChannel, "outred.txt", h, w);
    generateImageFile(testGreenChannel, "outgreen.txt", h, w);
    generateImageFile(testBlueChannel, "outblue.txt", h, w);

    ///////////////////////////////////////////////////////////
    ///////////////////////FREE MEMORY/////////////////////////
    ///////////////////////////////////////////////////////////

    destroyArray(dIcdx);
    destroyArray(dHdyT);
    destroyArray(testRedChannel);
    destroyArray(testGreenChannel);
    destroyArray(testBlueChannel);
    return 0;
}

