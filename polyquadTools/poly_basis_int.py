
from ctypes.wintypes import LPARAM
import numpy as np 
from math import pi, sin, cos, sqrt
import mpmath as mp

from scipy import special
from scipy.spatial import Delaunay
import scipy as sp


def gauss_quad(n, a=-1, b=1):
    #beta = .5./sqrt(1-(2*(1:n)).ˆ(-2))
    #alpha = mp.zeros(n, 1)
    #beta = mp.matrix([.5/sqrt(1-(2*(i))**(-2)) for i in range(1, n)], ctx=ctx)
    A = mp.zeros(n)
    for i in range(n-1):
        A[i, i+1] = A[i+1, i]= .5/sqrt(1-(2*(i+1))**(-2))
    x, v =  mp.eigsy(A)
    w = mp.matrix([(b - a)*v[0,i]**2 for i in range(v.cols)])

    return 0.5*(b - a)*x + 0.5*(b + a), w


def poly_area(X):
    x = X[:,0]
    y = X[:,1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def legendre_2d(i, j, x, y):
    return mp.legendre(i, x) * mp.legendre(j, y)


def vandermonde_2d(basis, X):
    V = mp.zeros(X.rows, basis.rows)
    for i in range(basis.rows):
        b = basis[i,:]
        for j in range(X.rows):
            V[j,i] = legendre_2d(b[0], b[1], X[j,0], X[j,1])
    return V


def interpolation_mat(basis, x1, x2):
    V1 = vandermonde_2d(basis, x1)
    V2 = vandermonde_2d(basis, x2)
    return V2 * V1**-1


def tensor_quad_duffy(nquad):
    xref, wref = gauss_quad(nquad, a=0, b=1)

    X = mp.zeros(nquad*nquad, 2)
    W = mp.zeros(nquad*nquad, 1)
    for j in range(nquad):
        for i in range(nquad):
            X[j*nquad + i, 0] = xref[i]
            X[j*nquad + i, 1] = xref[j]*xref[i]
            W[j*nquad + i, 0] = wref[i]*wref[j]
    return X, W


class Triangulation(object):
    def __init__(self, xn, ctx, deg):
        self.xn = xn
        self.ctx = ctx
        self.deg = deg
        self.xi, self.wi = self._int_points(xn, ctx, deg)

    def _int_points(self, xn, ctx, deg):
        A = poly_area(xn)
        basis = ctx.matrix([[0,0],[1,0],[0,1]])

        xt_ref = ctx.matrix([[0,0], [1,0], [1,1]])
        xi_ref, wi_ref = tensor_quad_duffy(deg)

        Mns = interpolation_mat(basis, xt_ref, xi_ref)

        xi_phy = Mns * xn
        wi_phy = wi_ref
        for i in range(wi_ref.rows):
            wi_phy[i, 0] *= 2 * A * xi_ref[i, 0]

        return xi_phy, wi_phy

    def integrate(self, f):
        I = 0
        for i in range(self.xi.rows):
            I += self.wi[i] * f(self.xi[i, 0], self.xi[i, 1])
        return I


class Polygon(object):
    def __init__(self, m, ctx, deg, xn=None):
        self.m = m
        self.ctx = ctx
        self.deg = deg
        xt = mp.zeros(3, 2)

        if xn is None:
            self.xn = mp.zeros(m, 2)
            xn_temp = np.zeros((m ,2))
            for i in range(m):
                self.xn[i,0] = ctx.cos(2*mp.pi*i/m)
                self.xn[i,1] = ctx.sin(2*mp.pi*i/m)
                xn_temp[i,:] = [cos(2*pi*i/m), sin(2*pi*i/m)]
        else:
            xn_temp = xn
            self.xn = mp.matrix(xn)

        self._tri = Delaunay(xn_temp)
        self.triang = []
        for s in self._tri.simplices:
            for i,n in enumerate(s):
                xt[i,0] = self.xn[n,0]
                xt[i,1] = self.xn[n,1]
            self.triang.append(Triangulation(xt, self.ctx, self.deg))

    def integrate(self, f):
        I = 0
        for t in self.triang:
            I += t.integrate(f)
        return I


def lp_basis(d, k_max, p):
    bt = np.zeros(((k_max+1)**d, d))
    if d == 2:
        for j in range(k_max + 1):
            for i in range(k_max + 1):
                bt[j*(k_max + 1) + i, 0] = i
                bt[j*(k_max + 1) + i, 1] = j
    
    return bt[np.linalg.norm(bt, ord=p, axis=1) <= k_max]


def mass_matrix(poly, basis):
    nbase = np.shape(basis)[0]
    M = mp.zeros(nbase, nbase)
    for i1, b1 in enumerate(basis):
        for i2, b2 in enumerate(basis):
            f = lambda x, y : legendre_2d(b1[0], b1[1], x, y)*legendre_2d(b2[0], b2[1], x, y)
            M[i1, i2] = poly.integrate(f)
    return M


def basis_int(poly, basis):
    nbase = np.shape(basis)[0]
    I = mp.zeros(nbase, 1)
    for i, b in enumerate(basis):
        f = lambda x, y : legendre_2d(b[0], b[1], x, y)
        I[i,0] = poly.integrate(f)
    return I



def norm_basis(M):
    L = mp.cholesky(M)
    I = L * M * L.T
    for i in range(L.rows):
        for j in range(L.cols):
            L[i,j] /= mp.sqrt(I[i, i])
    return L


def polyquad_basis(L, basis, tol=1e-16):
    src = ''

    src += 'std::vector<T> b_exact = {'

    for i in range(L.rows):
        for j in range(L.cols):
            src += f'T({mp.nstr(L[i,j], n=16)})'
            if i != L.rows - 1 or j != L.cols - 1:
                src += ',\n'
            else:
                src += '};'

    src += f'''
    
    '''

    return src

def orbits(m, ctx):
    dps = ctx.dps
    xn = mp.zeros(m, 2)
    for i in range(m):
        xn[i,0] = ctx.cos(2*mp.pi*i/m)
        xn[i,1] = ctx.sin(2*mp.pi*i/m)

    # Orbit 0
    src0 = 'pts.row(poff) = vector2d(T(0.), T(0.));'

    # Orbit 1
    src1 = ''
    for i in range(m):
        xc = 0.5*(xn[i,:] + xn[((i + 1) % m),:])
        src1 += f'pts.row(poff + {i}) = vector2d(a*T({mp.nstr(xc[0,0], n=dps)}), a*T({mp.nstr(xc[0,1], n=dps)}));\n'

    # Orbit 2
    src2 = ''
    for i in range(m):
        src2 += f'pts.row(poff + {i}) = vector2d(a*T({mp.nstr(xn[i,0], n=dps)}), a*T({mp.nstr(xn[i,1], n=dps)}));\n'

        
    # Orbit 3
    src3 = ''
    for i in range(m):
        xc1 = 0.5*(xn[i + 0,:] + xn[((i + 1) % m),:])
        xc2 = 0.5*(xn[i - 1,:] + xn[((i + 0) % m),:])

        src3 += f'pts.row(poff + {2*i}) = vector2d('
        src3 += f'b*a*T({mp.nstr(xc1[0,0], n=dps)}) + (1-b)*a*T({mp.nstr(xn[i,0], n=dps)}),'
        src3 += f'b*a*T({mp.nstr(xc1[0,1], n=dps)}) + (1-b)*a*T({mp.nstr(xn[i,1], n=dps)}));\n'
        src3 += f'pts.row(poff + {2*i + 1}) = vector2d('
        src3 += f'b*a*T({mp.nstr(xc2[0,0], n=dps)}) + (1-b)*a*T({mp.nstr(xn[i,0], n=dps)}),'
        src3 += f'b*a*T({mp.nstr(xc2[0,1], n=dps)}) + (1-b)*a*T({mp.nstr(xn[i,1], n=dps)}));\n'

    return src0, src1, src2, src3


def polyquad_class(m, k_max, I, ctx):
    o0, o1, o2, o3 = orbits(m, ctx)

    b_exact_src = ''
    for i in range(I.rows):
        b_exact_src += f'T({mp.nstr(I[i,0], n=ctx.dps)})'
        if i != I.rows-1:
            b_exact_src += ',\n'

    src = f'''
// Autogenerated polygon class for Polyquad
// sides = {m}
// k_max = {k_max}
// precision = {ctx.dps}

#ifndef POLYQUAD_SHAPES_POLY_HPP
#define POLYQUAD_SHAPES_POLY_HPP

#include "shapes/base.hpp"
#include "utils/ortho_poly.hpp"

#include <Eigen/Dense>

#include <array>
#include <cassert>
#include <iostream>
#include <vector>

namespace polyquad {{

template<typename T>
class PolyDomain : public BaseDomain<PolyDomain<T>, T, 2, 4>
{{
public:
    typedef BaseDomain<PolyDomain<T>, T, 2, 4> Base;
    typedef typename Base::MatrixXT MatrixXT;
    typedef typename Base::VectorXT VectorXT;
    typedef typename Base::ArrayXT ArrayXT;
    typedef typename Base::MatrixPtsT MatrixPtsT;
    typedef typename Base::VectorOrb VectorOrb;

    typedef Eigen::Matrix<T, 2, 1> Vector2T;

public:
    PolyDomain() : Base(T({mp.nstr(I[0,0], n=ctx.dps)}))
    {{}}

    void configure_d(int qdeg, bool poswts, const VectorOrb& orbits);

    static bool validate_orbit(const VectorOrb& orb)
    {{ return orb(0) <= 1; }}

private:
    friend class BaseDomain<PolyDomain<T>, T, 2, 4>;

    static constexpr int npts_for_orbit[] = {{1, {m}, {m}, {2*m}}};
    static constexpr int narg_for_orbit[] = {{0, 1, 1, 2}};
    static constexpr int nbfn_for_qdeg(int qdeg);

    void expand_orbit(int i, int aoff, int poff,
                      const VectorXT& args, MatrixPtsT& pts) const;

    void seed_orbit(int i, int aoff, VectorXT& args);

    template<typename D1, typename D2>
    void eval_orthob_block(const D1 pq, D2 out) const;

    template<typename ReplaceF>
    static void collapse_arg(int i, int aoff, const VectorXT& args, ReplaceF replace, const T& tol);

    static void clamp_arg(int i, int aoff, VectorXT& args);

    static void sort_arg(int i, int aoff, VectorXT& args);

private:
    static Vector2T vector2d(const T& p1, const T& p2)
    {{ return {{p1, p2}}; }}
}};

template<typename T>
inline void
PolyDomain<T>::configure_d(int qdeg,
        bool poswts,
        const VectorOrb& orbits)
{{
    std::vector<T> b_exact = {{
        {b_exact_src}
        }};

    for (int i=0, off = 0; i <= this->qdeg(); i += 1)
       {{
       for (int j=i; j <= this->qdeg()-i; j += 1, off++)
          {{
              this->b_(off) = b_exact[j*{k_max + 1} + i];
          }} 
       }}
}}

template<typename T>
inline constexpr int
PolyDomain<T>::nbfn_for_qdeg(int qdeg)
{{
    int n = 0;

    for (int i = 0; i <= qdeg; i += 1)
        for (int j = i; j <= qdeg - i; j += 1, ++n);

    return n;
}}

template<typename T>
void
PolyDomain<T>::expand_orbit(int i, int aoff, int poff,
                           const VectorXT& args, MatrixPtsT& pts) const
{{
    switch (i)
    {{
        case 0:
        {{
            {o0}
            break;
        }}
        case 1:
        {{
            const T& a = args(aoff + 0);
            {o1}
            break;
        }}
        case 2:
        {{
            const T& a = args(aoff + 0);
            {o2}
            break;
        }}
        case 3:
        {{
            const T& a = args(aoff + 0);
            const T& b = args(aoff + 1);
            {o3}
            break;
        }}
        default:
            assert(0 && "Bad orbit"), abort();
    }}
}}

template<typename T>
inline void
PolyDomain<T>::seed_orbit(int i, int aoff, VectorXT& args)
{{

    switch (i)
    {{
        case 0:
            break;
        case 1:
        case 2:
            args(aoff) = this->rand();
            break;
        case 3:
            args(aoff + 0) = this->rand();
            args(aoff + 1) = this->rand();
            break;
        default:
            assert(0 && "Bad orbit"), abort();
    }}
}}

template<typename T>
template<typename D1, typename D2>
inline void
PolyDomain<T>::eval_orthob_block(const D1 pq, D2 out) const
{{
    typedef Eigen::Array<T, D1::RowsAtCompileTime, 1> ArrayT;

    const auto& p = pq.col(0);
    const auto& q = pq.col(1);

    const ArrayT a = p;
    const ArrayT b = q;

    const T half = 0.5;

    JacobiP<ArrayT> jpa(0, 0, a);

    for (int i = 0, off = 0; i <= this->qdeg(); i += 1)
    {{
        JacobiP<ArrayT> jpb(0, 0, b);
        for (int j = i; j <= this->qdeg() - i; j += 1, ++off)
        {{
            out.row(off) = jpa(i)*jpb(j);
        }}
    }}
}}

template<typename T>
template<typename ReplaceF>
void inline
PolyDomain<T>::collapse_arg(int i, int aoff, const VectorXT& args,
                           ReplaceF replace, const T& tol)
{{
    
    if (i == 1 || i == 2 && abs(args(aoff)) < tol)
        replace(0);
    else if (i == 3)
    {{
        const T a = args(aoff + 0), b = args(aoff + 1);

        if (abs(a) < tol)
             replace(0);
        else if (abs(b) < tol)
            replace(2, a);
        else if (1 - abs(b) < tol)
            replace(1, a);
    }}
}}

template<typename T>
inline void
PolyDomain<T>::clamp_arg(int i, int aoff, VectorXT& args)
{{
    switch (i)
    {{
        case 0:
            break;
        case 1:
        case 2:
            args(aoff) = clamp(0, args(aoff), 1);
            break;
        case 3:
            args(aoff + 0) = clamp(0, args(aoff + 0), 1);
            args(aoff + 1) = clamp(0, args(aoff + 1), 1);
            break;
        default:
            assert(0 && "Bad orbit"), abort();
    }}
}}

template<typename T>
inline void
PolyDomain<T>::sort_arg(int i, int aoff, VectorXT& args)
{{
    
}}

}}   

#endif /* POLYQUAD_SHAPES_HEXA_HPP */
    '''
    return src


if __name__ == '__main__':
    sides = 6      # number of sides to polygon
    mp.mp.dps = 40 # number of decimal places
    nquad = 25     # number of points in 1d quadrature used in triangulation. 
    k_max = 20     # max order of polynomial basis
    lp = np.inf    # np.inf for infinity

    poly = Polygon(sides, mp.mp, nquad)
    basis = lp_basis(2, k_max, lp)

    I = basis_int(poly, basis)

    src = polyquad_class(sides, k_max, I, mp.mp)

    with open(f'polygon.hpp', 'w') as f:
        f.write(src)