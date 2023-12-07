#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

double **genWAM(double **matrix, int rows, int cols);
double **genDDG(double **wam, int n);
double **genGL(double **wam, double **ddg, int n);
double vectorsWeight(double *v, double *u, int cols);

/* jacobi */
#define ROTATION_LIMIT 100
#define EPSILON 0.00001
typedef struct rotationMatrixData
{
    int i;
    int j;
    double phi;
    double c;
    double s;
    double t;
} rotationMatrixData;
typedef struct eigenHolder
{
    double eigValue;
    double *eigenVector;
} eigenHolder;

void getMaxOffDiagonal(double **matrix, int rows, rotationMatrixData *rotationMatrix);
double calcPhi(double **matrix, rotationMatrixData *rotationMatrix);
double calcT(double phi);
double calcC(double t);
int sign(double);
double **getIdMatrix(int rows);
double **matrixProduct(double **matrixA, double **matrixB, int rows);
double **updateAtag(double **matrix, rotationMatrixData *rotationMatrix, int rows);
double **genRotationMatrix(rotationMatrixData *rotationMatrix, int rows);
double offDiagEuclidDistance(double **matrix, int rows);
eigenHolder **matrixToEig(double **matrix, double **eigenVectors, int rows);
eigenHolder **jacobi(double **matrix, int n);

/* printing */
void printMatrix(double **matrix, int rows, int columns);
void printEigs(eigenHolder **eigs, int rows);

/* memory */
void freeMatrix(double **matrix, int rows);
void freeEigenHolder(eigenHolder **eigenHolders, int n);