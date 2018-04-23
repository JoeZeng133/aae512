#pragma once

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>

#define FIR_1(f, i, j) (f(i + 1, j) - f(i - 1, j))
#define FIR_2(f, i, j) (f(i, j + 1) - f(i, j - 1))
#define SEC_1(f, i, j) (f(i + 1, j) + f(i - 1, j) - 2 * f(i, j))
#define SEC_2(f, i, j) (f(i, j + 1) + f(i, j - 1) - 2 * f(i, j))
#define SEC_MIX(f, i, j) (f(i + 1, j + 1) + f(i - 1, j - 1) - f(i + 1, j- 1) - f(i - 1, j + 1))
#define AVE_1(f, i, j) (f(i + 1, j) + f(i - 1, j))
#define AVE_2(f, i, j) (f(i, j - 1) + f(i, j + 1))
#define SQR(x) ((x) * (x))

class Mat {
private:
	std::vector<double> datareal;
	double* databuffer;
	size_t _size;
	size_t N, M;

public:
	Mat() = default;
	Mat(const size_t x, const size_t y) {
		resize(x, y);
	}
    Mat(const Mat& other) {
        datareal = other.datareal;
        databuffer = datareal.data();
        _size = other._size;
        N = other.N;
        M = other.M;
    }
    
	Mat& operator=(const Mat& other) {
		if (this != &other) {
			datareal = other.datareal;
			databuffer = datareal.data();
			_size = other._size;
			N = other.N;
			M = other.M;
		}
		return *this;
	}
	Mat& operator=(Mat&& other) {
		if (this != &other) {
			datareal = std::move(other.datareal);
			databuffer = other.databuffer;
			_size = other._size;
			N = other.N;
			M = other.M;
		}
		return *this;
	}
    
    ~Mat() = default;

	double* data() {
		return databuffer;
	}

	void resize(const size_t x, const size_t y) {
		M = y; N = x;
		_size = M * N;
		datareal.resize(_size);
		databuffer = datareal.data();
	}

	double& operator()(const int x, const int y) const{
		return databuffer[x * M + y];
	}
    
    void fill(const size_t x, const size_t y, const double val) {
        resize(x, y);
        std::fill(datareal.begin(), datareal.end(), val);
    }
    
    void fill(const double val) {
        std::fill(datareal.begin(), datareal.end(), val);
    }

	void zeros(const size_t x, const size_t y) {
		resize(x, y);
		std::fill(datareal.begin(), datareal.end(), 0);
	}

	void ones(const size_t x, const size_t y) {
		resize(x, y);
		std::fill(datareal.begin(), datareal.end(), 1);
	}

	void rand(const size_t x, const size_t y) {
		resize(x, y);
		for (auto& i : datareal) {
			i = double(std::rand()) / RAND_MAX;
		}
	}

	void getSize(int& n, int& m) const{
		n = N;
		m = M;
	}

	size_t getSize() const{
		return _size;
	}
};

double calcError(const Mat& A, const Mat& B);
void tridiag(double *aa, double *dd, double *cc, double *bb, double *x, int imax);
