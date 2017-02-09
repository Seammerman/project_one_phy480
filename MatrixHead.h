/*
header file to hold the matrix class for use in later projects as well as the current project 1
template found at (https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op)
Currently has no way to handle incorrect shapes
Must be used to define square matrices (M x M)
*/
class Matrix
{
public:
	Matrix(unsigned row, unsigned col);
	double& operator() (unsigned row, unsigned col);
	double operator()(unsigned row, unsigned col) const;
	//..
	~Matrix(); // Destructor
	Matrix(const Matrix& m); // copy constructor
	Matrix& operator= (const Matrix& m); // Assignment operator
										 // ...
private:
	unsigned rows_, cols_;
	double* data_;
};

inline
Matrix::Matrix(unsigned rows, unsigned cols) :rows_(rows), cols_(cols)
{
	data_ = new double[rows * cols];
}
inline
Matrix::~Matrix() {
	delete[] data_;
}
inline
double& Matrix:: operator()(unsigned row, unsigned col)
{
	return data_[cols_*row + col];
}
inline
double Matrix::operator()(unsigned row, unsigned col) const
{
	return data_[cols_*row + col];
}
