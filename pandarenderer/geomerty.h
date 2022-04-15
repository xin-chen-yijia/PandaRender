#pragma once
#include <cmath>
#include <cassert>

template<int n,typename T>
class vec
{
public:
	vec()=default;

	T& operator[](int index)
	{
		assert(index >= 0 && index < n);
		return data[index];
	}

	T operator[](int index) const
	{
		assert(index >= 0 && index < n);
		return data[index];
	}

	T norm2() { return *this * *this; }
	double norm() { return std::sqrt(norm2()); }

	T data[n] = { 0 };
};


template<typename T>
class vec<2,T>
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

	T norm2() { return *this * *this; }
	T norm() { return std::sqrt(norm2()); }

	vec& normalize()
	{
		*this = (*this) / norm();
		return *this;
	}

	T x{};
	T y{};
};

template<typename T>
class vec<3,T> {
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

	T norm2() { return *this * *this; }
	T norm() { return std::sqrt(norm2()); }

	vec& normalize()
	{
		*this = (*this) / norm();
		return *this;
	}

	T x{}, y{}, z{};
};

template<int n,typename T>
T operator*(const vec<n,T>& lhs, const vec<n,T>& rhs)
{
	T sum = 0;
	for (int i = 0;i < n;++i)
	{
		sum += lhs[i] * rhs[i];
	}

	return sum;
}

template<int n,typename T>
vec<n,T> operator+(const vec<n,T>& lhs, const vec<n,T>& rhs)
{
	vec<n,T> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs[i] + rhs[i];
	}

	return tmp;
}

template<int n,typename T>
vec<n,T> operator-(const vec<n,T>& lhs, const vec<n,T>& rhs)
{
	vec<n,T> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs[i] - rhs[i];
	}

	return tmp;
}

template<int n,typename T, typename U>
vec<n,T> operator*(const vec<n,T>& lhs, const U& rhs)
{
	vec<n,T> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs[i] * rhs;
	}

	return tmp;
}
template<int n, typename T>
vec<n,T> operator*(T lhs, const vec<n,T>& rhs)
{
	vec<n,T> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs * rhs[i];
	}

	return tmp;
}

template<int n, typename T, typename U>
vec<n,T> operator/(const vec<n,T>& lhs, const U& rhs)
{
	vec<n,T> tmp;
	for (int i = 0;i < n;++i)
	{
		tmp[i] = lhs[i] / rhs;
	}

	return tmp;
}

template<int n1, int n2, typename T>
vec<n1,T> embed(const vec<n2,T>& lhs, T fill = 1.0f)
{
	vec<n1,T> tmp;
	for (int i = 0;i < n1;++i)
	{
		tmp[i] = (i < n2 ? lhs[i] : fill);
	}

	return tmp;
}

template<int n1, int n2, typename T>
vec<n1,T> proj(const vec<n2,T>& v)
{
	vec<n1,T> tmp;
	for (int i = 0;i < n1;++i)
	{
		tmp[i] = v[i];
	}

	return tmp;
}


using vec2f = vec<2,float>;
using vec3f = vec<3,float>;
using vec4f = vec<4,float>;
using vec2i = vec<2, int>;
using vec3i = vec<3, int>;

vec3f cross(const vec3f& v1, const vec3f& v2);


template<int n,typename T> class dt;

template<int nrows,int ncols,typename T> 
struct mat
{
public:
	vec<ncols, T> data[nrows] = { {} };

	vec<ncols,T>& operator[](int i)
	{
		assert(i >= 0 && i < nrows);
		return data[i];
	}

	const vec<ncols,T>& operator[](int i) const
	{
		assert(i >= 0 && i < nrows);
		return data[i];
	}

	vec<nrows,T> col(const int idx) const
	{
		assert(idx >= 0 && idx < ncols);
		vec<nrows, T> tmp;
		for (int i = 0;i < nrows;++i)
		{
			tmp[i] = data[i][idx];
		}

		return tmp;
	}

	void set_col(int idx, const vec<nrows,T>& v)
	{
		assert(idx >= 0 && idx < ncols);
		for (int i = 0;i < nrows;++i)
		{
			data[i][idx] = v[i];
		}
	}

	static mat<nrows,ncols,T> identity()
	{
		mat<nrows, ncols,T> tmp = { 0 };
		for (int i = 0;i < nrows;++i)
		{
			tmp.data[i][i] = 1;
		}

		return tmp;
	}

	T det() const
	{
		return dt<ncols,T>::det(*this);
	}

	mat<nrows - 1, ncols - 1, T> get_minor(const int row, const int col) const
	{
		mat<nrows - 1, ncols - 1, T> ret;
		for (int i = 0;i < nrows - 1;++i)
		{
			for (int j = 0;j < ncols;++j)
			{
				ret[i][j] = data[i < row ? i : i+1][j < col ? j : j + 1];
			}
		}

		return ret;
	}

	T cofactor(const int row, const int col) const {
		return get_minor(row, col).det() * ((row + col) % 2 ? -1 : 1);
	}

	mat<nrows, ncols, T> adjugate() const {
		mat<nrows, ncols, T> ret;
		for (int i = 0;i < nrows;++i)
		{
			for (int j = 0;j < ncols;++j)
			{
				ret[i][j] = cofactor(i, j);
			}
		}

		return ret;
	}

	mat<nrows, ncols, T> invert_transpose() const
	{
		mat<nrows, ncols, T> ret = adjugate();
		return ret / (ret[0] * data[0]);
	}

	mat<nrows, ncols, T> invert() const
	{
		return invert_transpose().transpose();
	}

	mat<nrows, ncols, T> transpose() const
	{
		mat<ncols, nrows, T> ret;
		for (int i = 0;i < ncols;++i)
		{
			ret[i] = this->col(i);
		}

		return ret;
	}

};

template<int nrows, int ncols,typename T> 
vec<nrows,T> operator*(const mat<nrows, ncols,T>& lhs, const vec<ncols,T>& rhs)
{
	vec<nrows,T> ret;
	for (int i = 0;i < nrows;++i)
	{
		ret[i] = lhs[i] * rhs;
	}

	return ret;
}

template<int R1, int C1, int C2, typename T>
mat<R1, C2,T> operator*(const mat<R1, C1,T>& lhs, const mat<C1, C2,T>& rhs)
{
	mat<R1, C2,T> ret;
	for (int i = 0;i < R1;++i)
	{
		for (int j = 0;j < C2;++j)
		{
			ret[i][j] = lhs[i] * rhs.col(j);
		}
	}

	return ret;
}

template<int nrows,int ncols,typename T, typename U>
mat<nrows, ncols,T> operator*(const mat<nrows, ncols,T>& lhs, const U& val)
{
	mat<nrows, ncols, T> ret;
	for (int i = 0;i < nrows;++i)
	{
		ret[i] = lhs[i] * val;
	}

	return ret;
}

template<int nrows,int ncols, typename T, typename U>
mat<nrows, ncols,T> operator/(const mat<nrows, ncols,T>& lhs, const U& val)
{
	mat<nrows, ncols,T> ret;
	for (int i = 0;i < nrows;++i)
	{
		ret[i] = lhs[i] / val;
	}

	return ret;
}

template<int nrows,int ncols,typename T>
mat<nrows, ncols, T> operator+(const mat<nrows, ncols, T>& lhs, const mat<nrows, ncols,T>& rhs)
{
	mat<nrows, ncols,T> ret;
	for (int i = 0;i < nrows;++i)
	{
		for (int j = 0;j < ncols;++j)
		{
			ret[i][j] = lhs[i][j] + rhs[i][j];
		}
	}

	return ret;
}

template<int nrows,int ncols,typename T>
mat<nrows, ncols,T> operator-(const mat<nrows, ncols,T>& lhs, const mat<nrows, ncols,T>& rhs)
{
	mat<nrows, ncols,T> ret;
	for (int i = 0;i < nrows;++i)
	{
		for (int j = 0;j < ncols;++j)
		{
			ret[i][j] = lhs[i][j] - rhs[i][j];
		}
	}

	return ret;
}

template<int n,typename T>
class dt
{
public:
	static T det(const mat<n, n, T>& src)
	{
		T ret = 0.0;
		for (int i = 0;i < n;++i)
		{
			ret += src[0][i] * src.cofactor(0, i);
		}

		return ret;
	}
};

template<typename T>
class dt<1,T>
{
public:
	static T det(const mat<1,1,T>& src)
	{
		return src[0][0];
	}
};

using mat4x4 = mat<4, 4, float>;
