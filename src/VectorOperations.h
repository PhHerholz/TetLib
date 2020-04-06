#pragma once

#include <vector>
#include <cmath>
#include <iterator>


///////////////////////// constructors /////////////////////////////////

template<class T, int N>
std::array<T, N>
zeros(const T v = T(0))
{
    std::array<T, N> ret;
    ret.fill(v);
    return ret;
}

////////////////// std::array<T, N> operations /////////////////////////

template<class T, unsigned long N>
inline
std::array<T, N>
operator+(const std::array<T, N>& v0, const std::array<T, N>& v1)
{
    std::array<T, N> res;
    
    for(int i = 0; i < N; ++i)
    {
        res[i] = v0[i] + v1[i];
    }
    
    return res;
}


template<class T, unsigned long N>
inline
std::array<T, N>
operator*(const std::array<T, N>& v0, const std::array<T, N>& v1)
{
    std::array<T, N> res;
    
    for(int i = 0; i < N; ++i)
    {
        res[i] = v0[i] * v1[i];
    }
    
    return res;
}

template<class T, unsigned long N>
inline
std::array<T, N>
operator/(const std::array<T, N>& v0, const std::array<T, N>& v1)
{
    std::array<T, N> res;
    
    for(int i = 0; i < N; ++i)
    {
        res[i] = v0[i] / v1[i];
    }
    
    return res;
}


template<class T, unsigned long N>
inline
bool
operator==(std::array<T, N>& v0, const std::array<T, N>& v1)
{
    
    for(int i = 0; i < N; ++i)
    {
        if(v0[i] != v1[i]) return false;
    }
    
    return true;
}

template<class T, unsigned long N>
inline
std::array<T, N>&
operator+=(std::array<T, N>& v0, const std::array<T, N>& v1)
{
    for(int i = 0; i < N; ++i)
    {
        v0[i] += v1[i];
    }
    
    return v0;
}

template<class T, unsigned long N>
inline
std::array<T, N>&
operator-=(std::array<T, N>& v0, const std::array<T, N>& v1)
{
    for(int i = 0; i < N; ++i)
    {
        v0[i] -= v1[i];
    }
    
    return v0;
}


template<class T, unsigned long N>
inline
std::array<T, N>
operator/=(std::array<T, N>& v0, const double d)
{
    for(int i = 0; i < N; ++i)
    {
        v0[i] /= d;
    }
    
    return v0;
}

template<class T, unsigned long N>
inline
std::array<T, N>
operator*=(std::array<T, N>& v0, const double d)
{
    for(int i = 0; i < N; ++i)
    {
        v0[i] *= d;
    }
    
    return v0;
}

template<class T, unsigned long N>
inline
std::array<T, N>
operator*(const std::array<T, N>& v0, const double d)
{
    std::array<T, N> ret;
    
    for(int i = 0; i < N; ++i)
    {
        ret[i] = v0[i] * d;
    }
    
    return ret;
}

template<class T, unsigned long N>
inline
std::array<T, N>
operator*(const double d, const std::array<T, N>& v0)
{
    return v0 * d;
}

template<class T, unsigned long N>
inline
std::array<T, N>
operator/(const std::array<T, N>& v0, const double d)
{
    std::array<T, N> ret;
    
    for(int i = 0; i < N; ++i)
    {
        ret[i] = v0[i] / d;
    }
    
    return ret;
}

template<class T, unsigned long N>
inline
std::ostream&
operator<<( std::ostream &output, const std::array<T, N>& arr )
{
    for(int i = 0; i < N; ++i) output << arr[i] << " ";

    return output;
}
////////////////// std::vector<T> operations /////////////////////////


template<class T>
std::vector<T>&
operator+=(std::vector<T>& v0, const std::vector<T>& v1)
{
    const int n = v0.size();
    
    for(int i = 0; i < n; ++i)
    {
        v0[i] += v1[i];
    }
    
    return v0;
}

template<class T>
bool
operator<(const std::vector<T>& v0, const std::vector<T>& v1)
{
    using namespace std;
    
    auto func = [](const T t0, const T t1){return abs(t0) < abs(t1);};
    return *max_element(v0.begin(), v0.end(), func) < *max_element(v1.begin(), v1.end(), func);
}

template<class T>
bool
operator<(const std::vector<T>& v0, const T t)
{
    using namespace std;
    
    auto func = [](const T t0, const T t1){return abs(t0) < abs(t1);};
    return *max_element(cbegin(v0), cend(v0), func) < t;
}

template<class T>
bool
operator<(const T t, const std::vector<T>& v0)
{
    return !(v0 < t);
}

template<class T>
std::vector<T>
operator+(const std::vector<T>& v0, const std::vector<T>& v1)
{
    const int n = v0.size();
    std::vector<T> res(n);
    
    for(int i = 0; i < n; ++i)
    {
        res[i] = v0[i] + v1[i];
    }
    
    return res;
}

template<class T1, class T2, class TR = decltype(T1() * T2())>
std::vector<TR>
operator*(const std::vector<T1>& v0, const std::vector<T2>& v1)
{
    const int n = v0.size();
    std::vector<TR> res(n);
    
    for(int i = 0; i < n; ++i)
    {
        res[i] = v0[i] * v1[i];
    }
    
    return res;
}

template<class T>
std::vector<T>&
operator*=(std::vector<T>& v0, const double d)
{
    const int n = v0.size();
    
    for(int i = 0; i < n; ++i)
    {
        v0[i] *= d;
    }
    
    return v0;
}

template<class T>
std::vector<T>&
operator*(const std::vector<T>& v0, const double d)
{
    const int n = v0.size();
    std::vector<T> ret(n);
    
    for(int i = 0; i < n; ++i)
    {
        ret[i] = v0[i] * d;
    }
    
    return ret;
}

template<class T>
std::vector<T>&
operator*(const double d, const std::vector<T>& v0)
{
    return v0 * d;
}

template<class T>
std::vector<T>&
operator/=(std::vector<T>& v0, const double d)
{
    const int n = v0.size();
    
    for(int i = 0; i < n; ++i)
    {
        v0[i] /= d;
    }
    
    return v0;
}

template<class T>
std::vector<T>&&
operator/=(std::vector<T>&& v0, const double d)
{
    return move(v0 /= d);
}

template<class T>
std::vector<T>
operator-(const std::vector<T>& v0, const std::vector<T>& v1)
{
    const int n = v0.size();
    std::vector<T> res(n);
    
    for(int i = 0; i < n; ++i)
    {
        res[i] = v0[i] - v1[i];
    }
    
    return res;
}

template<class T>
std::vector<T>
operator-=(std::vector<T>& v0, const std::vector<T>& v1)
{
    const int n = v0.size();
    
    for(int i = 0; i < n; ++i)
    {
        v0[i] -= v1[i];
    }
    
    return v0;
}

template<class T>
double
total(std::vector<T>& d)
{
    T t(0);
    for(auto x : d) t += x;
    
    return t;
}

template<class T>
std::vector<std::vector<T>>
transpose(std::vector<std::vector<T>>& v)
{
    const int n = v.front().size();
    const int m = v.size();
    
    std::vector<std::vector<T>> ret(n, std::vector<T>(m));
    
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < m; ++j)
            ret[i][j] = v[j][i];
    
    return ret;
}

template<class T>
std::vector<T>
pow(std::vector<T> vec, const int k = 2)
{
    for(auto& d : vec) d = pow(d, k);
    return vec;
}

template<class T>
T norm2(const std::vector<T>& vec)
{
    const int n = vec.size();
    
    double r = .0;
    
    for(int i = 0; i < n; ++i) r += pow(vec[i], 2);
    return r;
}


template<class T>
T dist2(const std::vector<T>& vec, const std::vector<T>& vec2)
{
    const int n = vec.size();
    assert(n == vec2.size());
    
    double r = .0;
    
    for(int i = 0; i < n; ++i) r += std::pow(vec[i] - vec2[i], 2);
    return r;
}

template<class T>
T dist(const std::vector<T>& vec, const std::vector<T>& vec2)
{
    return sqrt(dist2(vec, vec2));
}


