#include "../include/matrix.h"


Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}
 
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();
 
    int k = 0;
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

Matrix::Matrix(int v_size)
{
    if (v_size <= 0)
    {
        cout << "Vector create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
    }
    this->fil = 1;
    this->col = v_size;
    this->matrix = (double **)malloc(v_size * sizeof(double *));

    if (this->matrix == NULL)
    {
        cout << "Vector create: error en la matriz\n";
        exit(EXIT_FAILURE);
    }

    this->matrix[0] = (double *)malloc(col * sizeof(double));

    for (int j = 0; j < v_size; j++)
    {
        this->matrix[0][j] = 0;
    }

}
 
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}
 
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}
 
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}
 
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];
 
    return *this;
}
 
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}

Matrix Matrix::operator/(const double div)
{
    Matrix result(this->fil, this->col);
 
    for (int i = 0; i < this->fil ; i++){
            for (int k = 0; k < this->col; k++){
                matrix[i][k] = matrix[i][k]/div;
            }
    }
 
    return result;
}

Matrix Matrix::operator*(const double mul)
{
    Matrix result(fil, col);
 
    for (int i = 0; i < this->fil ; i++){
            for (int k = 0; k < this->col; k++){
                matrix[i][k] = matrix[i][k]*mul;
            }
    }
 
    return result;
}
 
 
double & Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}

double & Matrix::operator()(const int n)
{
    if (n <= 0 || n > this->col)
    {
        cout << "Vector get: error límite de columna superado";
        exit(EXIT_FAILURE);
    }
    return this->matrix[0][n - 1];
}
 
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int Matrix::n_rows(){
    return fil;
}


int Matrix::n_cols(){
    return col;
}

void Matrix::insertCol(Matrix m, int col){
    if(m.n_cols()!=1){
        cout << "Error. Más de una columna insertada";
        exit(EXIT_FAILURE);
    }
    for(int i=0;i<m.n_rows(); i++){
        matrix[i][col]=m(i+1,1);
    }
    
}

Matrix Matrix::fila(int num_fila){
    Matrix fila(this->n_cols());
    for(int i = 1 ; i <= col ; i++){
        fila(i)=matrix[num_fila][i];
    }
    return fila;
}

Matrix Matrix::columna(int num_columna){
    Matrix columna(this->n_rows());
    for(int i = 1 ; i <= fil ; i++){
        columna(i)=matrix[i][num_columna];
    }
    return columna;
}

double Matrix::norm(){
    if(fil>1){
        printf("ERROR");
    }
    double sum=0.0;
    for(int i=0 ; i<col;i++){
        sum=sum+(matrix[0][i])*(matrix[0][i]);
    }
    return sqrt(sum);
}

double Matrix::det(){
    if (this->n_cols()!=this->n_rows() || (this->n_cols() != 1 ||  this->n_rows() != 3))
    {
        cout << "Matrix det: La matriz no es cuadrada";
        exit(EXIT_FAILURE);
    }
    double det=0;
    if(this->n_cols() == 3)
    {
        det = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) 
                    - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) 
                    + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
        
    }else if(this->n_cols() == 1){
        double det = matrix[0][0];
    }
    return det;
}

double dot(Matrix & m1, Matrix & m2){    
        double sum = 0.0;
    if(m1.n_rows()>1 || m2.n_rows()>1 || m1.n_cols()!=m2.n_cols()){
        EXIT_FAILURE;
    }else{
        for(int i = 1; i<= m1.n_cols();i++){
                sum = sum + m1(i)*m2(i);            
        }
    }
    return sum;
}

Matrix cross(Matrix & m1, Matrix & m2){
    if (m1.n_rows() != 1 || m2.n_rows() != 1 || m1.n_cols() != 3 || m2.n_cols() != 3)
    {
        cout << "Matrix cross: Error. No se puede hacer la operación cross con los vectores entrada";
        exit(EXIT_FAILURE);
    }

    Matrix result(3);

    result(1) = m1(2) * m2(3) - m1(3) * m2(2);
    result(2) = -(m1(1) * m2(3) - m1(3)* m2(1));
    result(3) = m1(1) * m2(2)- m1(2) * m2(1);

    return result;
}

Matrix zeros(const int n_row, const int n_col){
    if (n_row <= 0 || n_col <= 0)
    {
        cout << "Matrix zeros: error en row/cols";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(n_row, n_col);
    return *result;
}

Matrix transpose(Matrix m)
{
    Matrix *aux = new Matrix(m.n_cols(), m.n_rows());

    for (int i = 1; i <= m.n_rows(); i++)
    {
        for (int j = 1; j <= m.n_cols(); j++)
        {
            (*aux)(j,i) = m(i,j);
        }
    }

    return *aux;
}

Matrix inv(Matrix m){
    if(m.n_rows()!=m.n_cols()){
        cout << "Error inv. La matriz no es cuadrada.";
        exit(EXIT_FAILURE);
    }
    if(m.n_rows()!=3 || m.n_rows()!=1){
        cout << "Error inv. La matriz no es 3x3 o 1x1.";
        exit(EXIT_FAILURE);
    }
    double determinante = m.det();
    if(determinante==0.0){
        cout << "Error inv. La matriz no es invertible.";
        exit(EXIT_FAILURE);
    }
    if(m.n_rows()==3){
        Matrix inversa = zeros(3,3);
        inversa(1,1) = (m(2,2) * m(3,3) - m(2,3) * m(3,2)) / determinante;
        inversa(1,2) = (m(1,3) * m(3,2) - m(1,2) * m(3,3)) / determinante;
        inversa(1,3) = (m(1,2) * m(2,3) - m(1,3)* m(2,2)) / determinante;

        inversa(2,1) = (m(2,3) * m(3,1) - m(2,1) * m(3,3)) / determinante;
        inversa(2,2) = (m(1,1) * m(3,3) - m(1,3) * m(3,1)) / determinante;
        inversa(2,3) = (m(1,3) * m(2,1) - m(1,1) * m(2,3)) / determinante;

        inversa(3,1) = (m(2,1) * m(3,2) - m(2,2) * m(3,1)) / determinante;
        inversa(3,2) = (m(1,2) * m(3,1) - m(1,1) * m(3,2)) / determinante;
        inversa(3,3) = (m(1,1) * m(2,2) - m(1,2) * m(2,1)) / determinante;

        return inversa;

    }else{
        Matrix inversa = zeros(1,1);
        inversa(1,1)=1/m(1,1);
        return inversa;
    }

    
}

Matrix eye(const int size)
{
    if (size <= 0)
    {
        cout << "Matrix eye: error de tamaño";
        exit(EXIT_FAILURE);
    }

    Matrix result(size, size);

    for (int i = 1; i <= size; i++)
        result(i,i) = 1;

    return result;
}

double* vectorToArray(Matrix& v) {
    double* a = new double[v.n_cols()];
    for (int i = 0; i < v.n_cols(); i++) {
        a[i] = v(1,i+1);
    }
    return a;
}

Matrix arrayToVector(double a[], int n){
    Matrix m (n);
    for(int i = 0; i<n; i++){
        m(1,i+1) = a[i];
    }
    return m;
}