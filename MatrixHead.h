/*
Currently has no way to handle incorrect shapes
Must be used to define square matrices (M x M)

header file to hold the matrix class for use in later projects as well as the current project 1
template found at (https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op)

*/

class Matrix
{
public:
	Matrix(unsigned row, unsigned col);
	double& operator() (unsigned row, unsigned col);
	double operator() (unsigned row, unsigned col) const;
	//..
	~Matrix(); // Destructor
	Matrix(const Matrix& m); // copy constructor
	Matrix& operator= (const Matrix& m); // Assignment operator
	
									 // ...
private:
	unsigned rows_, cols_;
	double* data_;
};


Matrix::Matrix(unsigned rows, unsigned cols) 
	:rows_(rows)
	,cols_(cols)
{
	data_ = new double[rows * cols];
}

Matrix::~Matrix() {
	delete[] data_;
}

double& Matrix::operator()(unsigned row, unsigned col)
{
	return data_[row*cols_ + col];
}

double Matrix::operator()(unsigned row, unsigned col) const
{
	return data_[row*cols_ + col];
}


/*
class Matrix
{
public:
	Matrix(size_t rows, size_t cols);
	double& operator()(size_t i, size_t j);
	double operator()(size_t i, size_t j) const;
	void print()const;
private:
	size_t mRows;
	size_t mCols;
	std::vector<double> mData;
};

Matrix::Matrix(size_t rows, size_t cols)
	: mRows(rows),
	mCols(cols),
	mData(rows * cols)
{
}

double& Matrix::operator()(size_t i, size_t j)
{
	return mData[i * mCols + j];
}

double Matrix::operator()(size_t i, size_t j) const
{
	return mData[i * mCols + j];
}
void Matrix::print()const {
	for (int j = 0; j < mRows; j++) {
		for (int k = 0; k < mCols; k++) {
			std::cout << j <<" , "<< k << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
*/