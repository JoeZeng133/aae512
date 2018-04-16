#pragma once

#include <cstdlib>
#include <vector>
#include <algorithm>

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

	size_t getSize() {
		return _size;
	}
};