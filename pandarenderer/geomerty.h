#pragma once
#include <cmath>
#include <cassert>

#ifdef DOUBLE
#define Float double
#else
#define Float float
#endif // DOUBLE


template<typename T, int n>
class vec
{
public:
	vec()=default;

	T& operator[](int index)
	{
		assert(index >= 0 && index < n);
		return data_[index];
	}

	T operator[](int index) const
	{
		assert(index >= 0 && index < n);
		return data_[index];
	}

	T norm2() { return *this * *this; }
	double norm() { return std::sqrt(norm2()); }

private:
	T data_[n] = { 0 };
};


template<typename T>
class vec<T,2>
{
public:
	vec() = default;
	vec(T x, T y) :x(x), y(y)
	{

	}

	T& operator[](int i)
	{
		assert(i >= 0 && i < 2);
		return i == 0 ? x : y;
	}

	T operator[](int i) const
	{
		assert(i >= 0 && i < 2);
		return i == 0 ? x : y;
	}

	Float norm2() { return *this * *this; }
	Float norm() { return std::sqrt(norm2()); }

	vec& normalize()
	{
		*this = (*this) / norm();
		return *this;
	}

	T x{};
	T y{};
};

template<typename T>
class vec<T,3> {
public:
	vec() = default;
	vec(T x, T y, T z) :x(x), y(y), z(z)
	{

	}

	T& operator[](int i)
	{
		assert(i >= 0 && i < 3);
		return (i == 0 ? x : (i == 1 ? y : z));
	}

	T operator[](int i) const
	{
		assert(i >= 0 && i < 3);
		return (i == 0 ? x : (i == 1 ? y : z));
	}

	Float norm2() { return *this * *this; }
	Float norm() { return std::sqrt(norm2()); }

	vec& normalize()
	{
		*this = (*this) / norm();
		return *this;
	}

	T x{}, y{}, z{};
};

template<typename T, int n>
T operator*(const vec<T, n>& lhs, const vec<T, n>& rhs)
{
	T sum = 0;
	for (int i = 0;i < n;++i)
	{
		sum += lhs[i] * rhs[i];
	}

	return sum;
}

template<typename T, int n>
vec<T, n> operator+(const vec<T, n>& lhs, const vec<T, n>& rhs)
{
	vec<T, n> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs[i] + rhs[i];
	}

	return tmp;
}

template<typename T, int n>
vec<T, n> operator-(const vec<T, n>& lhs, const vec<T, n>& rhs)
{
	vec<T, n> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs[i] - rhs[i];
	}

	return tmp;
}

template<typename T, int n>
vec<T, n> operator*(const vec<T, n>& lhs, Float rhs)
{
	vec<T, n> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs[i] * rhs;
	}

	return tmp;
}
template<typename T, int n>
vec<T, n> operator*(Float lhs, const vec<T, n>& rhs)
{
	vec<T, n> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs * rhs[i];
	}

	return tmp;
}

template<typename T, int n>
vec<T, n> operator/(const vec<T, n>& lhs, Float rhs)
{
	vec<T, n> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs[i] / rhs;
	}

	return tmp;
}

template<typename T, int n1, int n2>
vec<T, n1> embed(const vec<T, n2>& lhs, T fill = 1.0f)
{
	vec<T, n1> tmp;
	for (int i = 0;i < n1;++i)
	{
		tmp[i] = (i < n2 ? lhs[i] : fill);
	}

	return tmp;
}

template<typename T, int n1, int n2>
vec<T, n1> proj(const vec<T, n2>& v)
{
	vec<T, n1> tmp;
	for (int i = 0;i < n1;++i)
	{
		tmp[i] = v[i];
	}

	return tmp;
}


using vec2 = vec<Float,2>;
using vec3 = vec<Float,3>;
using vec2i = vec<int, 2>;
using vec3i = vec<int, 3>;

vec3 cross(const vec3& v1, const vec3& v2);


template<int n> class dt;

template<int nrows,int ncols> 
class mat
{
public:
	vec<Float, ncols> data[nrows] = { {} };

	vec<Float, ncols>& operator[](int i)
	{
		assert(i >= 0 && i < nrows);
		return data[i];
	}

	const vec<Float, ncols>& operator[](int i) const
	{
		assert(i >= 0 && i < nrows);
		return data[i];
	}

	vec<Float, nrows> col(const int idx) const
	{
		assert(idx >= 0 && idx < ncols);
		vec<Float, nrows> tmp;
		for (int i = 0;i < nrows;++i)
		{
			tmp[i] = data[i][idx];
		}

		return tmp;
	}

	void set_col(int idx, const vec<Float, nrows>& v)
	{
		assert(idx >= 0 && idx < ncols);
		for (int i = 0;i < nrows;++i)
		{
			data[i][idx] = v[i];
		}
	}

	static mat<nrows,ncols> identity()
	{
		mat<nrows, ncols> tmp = { 0 };
		for (int i = 0;i < nrows;++i)
		{
			data[i][i] = 1;
		}

		return tmp;
	}

	Float det() const
	{
		return dt<ncols>::det(*this);
	}

	mat<nrows - 1, ncols - 1> get_minor(const int row, const int col) const
	{
		mat<nrows - 1, ncols - 1> ret;
		for (int i = 0;i < nrows - 1;++i)
		{
			for (int j = 0;j < ncols;++j)
			{
				ret[i][j] = data[i < row ? i : i+1][j < col ? j : j + 1];
			}
		}

		return ret;
	}

	Float cofactor(const int row, const int col) const {
		return get_minor(row, col).det() * ((row + col) % 2 ? -1 : 1);
	}

	mat<nrows, ncols> adjugate() const {
		mat<nrows, ncols> ret;
		for (int i = 0;i < nrows;++i)
		{
			for (int j = 0;j < ncols;++j)
			{
				ret[i][j] = cofactor(i, j);
			}
		}

		return ret;
	}

	mat<nrows, ncols> invert_transpose() const
	{
		mat<nrows, ncols> ret = adjugate();
		return ret / (ret[0] * data[0]);
	}

	mat<nrows, ncols> invert() const
	{
		return invert_transpose().transpose();
	}

	mat<nrows, ncols> transpose() const
	{
		mat<ncols, nrows> ret;
		for (int i = 0;i < ncols;++i)
		{
			ret[i] = this->col(i);
		}

		return ret;
	}

};

template<int nrows, int ncols> 
vec<Float,nrows> operator*(const mat<nrows, ncols>& lhs, const vec<Float, ncols>& rhs)
{
	vec<Float, nrows> ret;
	for (int i = 0;i < nrows;++i)
	{
		ret[i] = lhs[i] * rhs;
	}

	return ret;
}

template<int R1, int C1, int C2>
mat<R1, C2> operator*(const mat<R1, C1>& lhs, const mat<C1, C2>& rhs)
{
	mat<R1, C2> ret;
	for (int i = 0;i < R1;++i)
	{
		for (int j = 0;j < C2;++j)
		{
			ret[i][j] = lhs[i] * rhs.col(j);
		}
	}

	return ret;
}

template<int nrows,int ncols>
mat<nrows, ncols> operator*(const mat<nrows, ncols>& lhs, Float val)
{
	mat<nrows, ncols> ret;
	for (int i = 0;i < nrows;++i)
	{
		ret[i] = lhs[i] * val;
	}

	return ret;
}

template<int nrows,int ncols>
mat<nrows, ncols> operator/(const mat<nrows, ncols>& lhs, Float val)
{
	mat<nrows, ncols> ret;
	for (int i = 0;i < nrows;++i)
	{
		ret[i] = lhs[i] / val;
	}

	return ret;
}

template<int nrows,int ncols>
mat<nrows, ncols> operator+(const mat<nrows, ncols>& lhs, const mat<nrows, ncols>& rhs)
{
	mat<nrows, ncols> ret;
	for (int i = 0;i < nrows;++i)
	{
		for (int j = 0;j < ncols;++j)
		{
			ret[i][j] = lhs[i][j] + rhs[i][j];
		}
	}

	return ret;
}

template<int nrows,int ncols>
mat<nrows, ncols> operator-(const mat<nrows, ncols>& lhs, const mat<nrows, ncols>& rhs)
{
	mat<nrows, ncols> ret;
	for (int i = 0;i < nrows;++i)
	{
		for (int j = 0;j < ncols;++j)
		{
			ret[i][j] = lhs[i][j] - rhs[i][j];
		}
	}

	return ret;
}

template<int n>
class dt
{
public:
	static Float det(const mat<n, n>& src)
	{
		Float ret = 0.0;
		for (int i = 0;i < n;++i)
		{
			ret += src[0][i] * src.cofactor(0, i);
		}

		return ret;
	}
};

template<>
class dt<1>
{
public:
	static Float det(const mat<1,1>& src)
	{
		return src[0][0];
	}
};
