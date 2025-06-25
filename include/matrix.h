#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

class Matrix
{
    public:
        Matrix(const int fil, const int col);
        Matrix(const int fil, const int col, double v[], int n);
        Matrix(const int v_size);
        Matrix(const Matrix& m);
        ~Matrix();
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2);
        Matrix  operator-(const Matrix& matrix2);
        Matrix  operator*(const Matrix& matrix2);
        Matrix operator/(const double div);
        Matrix operator*(const double mul);
        double &operator()(const int i, const int j) const;
        double &operator()(const int n);
    
        int n_rows();
        int n_cols();
        Matrix fila(const int num_fila);
        Matrix columna(const int num_col);
        void insertCol(Matrix m, int col);
        
        void print();
        double norm();
        double det();

        private:
        void initMatrix();

        private:
        int fil;
        int col;
        double **matrix;
};
        double dot(Matrix & m1, Matrix & m2);
        Matrix zeros(const int n_row, const int n_col);
        Matrix transpose(Matrix m);
        Matrix inv(Matrix m);
        Matrix cross(Matrix & m1, Matrix & m2);
        Matrix eye(const int size);
        double* vectorToArray(Matrix& v);
        Matrix arrayToVector(double a[], int n);
 
    


#endif