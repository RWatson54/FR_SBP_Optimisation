/*
    This file is part of polyquad.
    Copyright (C) 2014  Freddie Witherden <freddie@witherden.org>

    Polyquad is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Polyquad is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with polyquad.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef POLYQUAD_SHAPES_BASE_HPP
#define POLYQUAD_SHAPES_BASE_HPP

#include <Eigen/Dense>
#include <Eigen/LevenbergMarquardt>

#include <algorithm>
#include <random>
#include <tuple>
#include <vector>

#include <stdio.h>

namespace polyquad {

template<typename T1, typename T2, typename T3>
T2 clamp(const T1& l, const T2& v, const T3& h)
{
    return (v < l) ? l : (v > h) ? h : v;
}

template<typename Derived, typename T, int Ndim, int Norbits>
class BaseDomain
{
public:
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;
    typedef Eigen::Array<T, Eigen::Dynamic, 1> ArrayXT;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Ndim> MatrixPtsT;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixObatT;
    typedef Eigen::Matrix<int, Norbits, 1> VectorOrb;
    typedef std::pair<VectorOrb, VectorXT> VectorOrbArgs;

    mutable VectorXT b_;

public:
    BaseDomain(const T& f0)
        : f0_(f0)
        , rand_eng_(std::random_device()())
    {}

    void configure(int qdeg, bool poswts, const VectorOrb& orbits);

    void seed();

    void seed(const VectorXT& args)
    { args_ = args; }

    VectorXT clamp_args(const VectorXT& args) const;

    void expand(const VectorXT& args, MatrixPtsT& pts) const;

    void eval_orthob(const MatrixPtsT& pts, MatrixObatT& out) const;

    void expand_wts(const VectorXT& wargs, VectorXT& wts) const;

    std::tuple<T, VectorXT> minimise(int maxfev);

    VectorXT wts(const VectorXT& args, VectorXT* resid=nullptr, VectorXT* cnorms=nullptr);

    int qdeg() const
    { return qdeg_; }

    int ndim() const
    { return Ndim; }

    int nbfn() const
    { return Derived::nbfn_for_qdeg(qdeg_); }

    const VectorOrb& orbits() const
    { return orbits_; }

    int npts() const
    { return npts(orbits_); }

    int nwts() const
    { return orbits_.sum(); }

    int ndof() const;

    std::vector<VectorOrbArgs> possible_reductions(const VectorXT& args,
                                                   const T& tol=1e-2);

protected:
    double rand(double a=0, double b=1)
    { return std::uniform_real_distribution<double>(a, b)(rand_eng_); }

    template<std::size_t N>
    double rand(double a, double b, const std::array<int, N>& wts);

    int arg_offset(int i, int j=0) const
    { return arg_offset(orbits_, i, j); }

public:
    static std::vector<VectorOrb> symm_decomps(int npts)
    { return symm_decomps(npts, VectorOrb::Constant(npts)); }

    static std::vector<VectorOrb> symm_decomps(int npts, const VectorOrb& limits);

private:
    static constexpr int npts(const VectorOrb& orb)
    {
        return std::inner_product(orb.data(), orb.data() + Norbits,
                                  Derived::npts_for_orbit, 0);
    }

    static constexpr int arg_offset(const VectorOrb& orb, int i, int j=0)
    {
        return std::inner_product(orb.data(), orb.data() + i,
                                  Derived::narg_for_orbit, 0)
             + Derived::narg_for_orbit[i]*j;
    }

    static VectorOrbArgs remove_orbit(const VectorOrb& orb, const VectorXT& args,
                                      int i, int j);

    static VectorOrbArgs add_orbit(const VectorOrb& orb, const VectorXT& args,
                                   int i, const std::initializer_list<T> aiargs);

    static void sort_args(const VectorOrb& orb, VectorXT& args);

    static void symm_decomps_recurse(VectorOrb coeffs,
                                     const VectorOrb& limits,
                                     int sum,
                                     VectorOrb& partsoln,
                                     std::vector<VectorOrb>& solns);

private:
    int qdeg_;
    bool poswts_;
    const T f0_;

    VectorOrb orbits_;
    VectorXT args_;

    std::mt19937 rand_eng_;

    mutable MatrixPtsT pts_;
    mutable MatrixObatT obat_;
    mutable MatrixXT A_;
    mutable VectorXT wts_;
};

template<typename Derived, typename T, int Ndim, int Norbits>
inline void
BaseDomain<Derived, T, Ndim, Norbits>::configure(
        int qdeg,
        bool poswts,
        const VectorOrb& orbits)
{
    Derived& derived = static_cast<Derived&>(*this);

    qdeg_ = qdeg;
    poswts_ = poswts;
    orbits_ = orbits;

    // Pre-allocate some scratch space
    pts_.resize(npts(), Ndim);
    obat_.resize(nbfn(), npts());
    A_.resize(nbfn(), nwts());
    b_.resize(nbfn());
    wts_.resize(nwts());

    b_.fill(0);
    b_(0) = f0_;

    // Run and domain specific configs
    derived.configure_d(qdeg, poswts, orbits);

}

template<typename Derived, typename T, int Ndim, int Norbits>
inline auto
BaseDomain<Derived, T, Ndim, Norbits>::symm_decomps(
    int npts,
    const VectorOrb& limits) -> std::vector<VectorOrb>
{
    VectorOrb coeffs(Derived::npts_for_orbit);

    std::vector<VectorOrb> solns;
    VectorOrb partsoln = VectorOrb::Zero();

    symm_decomps_recurse(coeffs, limits, npts, partsoln, solns);

    return solns;
}

template<typename Derived, typename T, int Ndim, int Norbits>
[[gnu::flatten]] inline void
BaseDomain<Derived, T, Ndim, Norbits>::seed()
{
    Derived& derived = static_cast<Derived&>(*this);
    args_ = VectorXT::Zero(ndof());

    for (int i = 0, aoff = 0; i < Norbits; ++i)
        for (int j = 0; j < orbits_(i); ++j, aoff += derived.narg_for_orbit[i])
            derived.seed_orbit(i, aoff, args_);
}

template<typename Derived, typename T, int Ndim, int Norbits>
[[gnu::flatten]] inline auto
BaseDomain<Derived, T, Ndim, Norbits>::clamp_args(const VectorXT& args) const -> VectorXT
{
    VectorXT nargs = args;

    for (int i = 0, aoff = 0; i < Norbits; ++i)
        for (int j = 0; j < orbits_(i); ++j, aoff += Derived::narg_for_orbit[i])
            Derived::clamp_arg(i, aoff, nargs);

    return nargs;
}

template<typename Derived, typename T, int Ndim, int Norbits>
[[gnu::flatten]] inline void
BaseDomain<Derived, T, Ndim, Norbits>::expand(
        const VectorXT& args,
        MatrixPtsT& pts) const
{
    const Derived& derived = static_cast<const Derived&>(*this);
    int aoff = 0, poff = 0;

    for (int i = 0; i < Norbits; ++i)
    {
        int ainc = derived.narg_for_orbit[i];
        int pinc = derived.npts_for_orbit[i];

        for (int j = 0; j < orbits_(i); ++j, aoff += ainc, poff += pinc)
            derived.expand_orbit(i, aoff, poff, args, pts);
    }
}

template<typename Derived, typename T, int Ndim, int Norbits>
[[gnu::flatten]] inline void
BaseDomain<Derived, T, Ndim, Norbits>::eval_orthob(
        const MatrixPtsT& pts,
        MatrixObatT& out) const
{
    const Derived& derived = static_cast<const Derived&>(*this);
    int n = out.cols();

    if (std::is_fundamental<T>::value)
    {
        int n16 = (n / 16)*16, n8 = (n / 8)*8, n4 = (n / 4)*4, n2 = (n / 2)*2;

        // Process points in blocks of sixteen
        for (int i = 0; i < n16; i += 16)
            derived.eval_orthob_block(pts.array().template middleRows<16>(i),
                                      out.array().template middleCols<16>(i));

        // Cleanup
        if (n16 < n8)
            derived.eval_orthob_block(pts.array().template middleRows<8>(n16),
                                      out.array().template middleCols<8>(n16));
        if (n8 < n4)
            derived.eval_orthob_block(pts.array().template middleRows<4>(n8),
                                      out.array().template middleCols<4>(n8));
        if (n4 < n2)
            derived.eval_orthob_block(pts.array().template middleRows<2>(n4),
                                      out.array().template middleCols<2>(n4));
        if (n2 < n)
            derived.eval_orthob_block(pts.array().template middleRows<1>(n2),
                                      out.array().template middleCols<1>(n2));
    }
    else
    {
        derived.eval_orthob_block(pts.array(), out.array());
    }
}

template<typename Derived, typename T, int Ndim, int Norbits>
inline void
BaseDomain<Derived, T, Ndim, Norbits>::expand_wts(
        const VectorXT& wargs,
        VectorXT& wtsout) const
{
    const Derived& derived = static_cast<const Derived&>(*this);

    for (int i = 0, poff = 0, woff = 0; i < Norbits; ++i)
    {
        int pinc = derived.npts_for_orbit[i];

        for (int j = 0; j < orbits_(i); ++j, poff += pinc, ++woff)
            wtsout.segment(poff, pinc).fill(wargs(woff));
    }
}


template<typename Derived, typename T, int Ndim, int Norbits>
inline auto
BaseDomain<Derived, T, Ndim, Norbits>::minimise(int maxfev) -> std::tuple<T, VectorXT>
{
    struct min_functor : Eigen::DenseFunctor<T>
    {
        min_functor(Derived& dom, bool poswts)
            : Eigen::DenseFunctor<T>(dom.ndof(), dom.nbfn() + 1 + poswts)
            , dom_(dom)
            , poswts_(poswts)
        {}

        int operator()(const VectorXT& x, VectorXT& f) const
        {
            // Compute the residual
            auto wts = dom_.wts(x, &f);

            // Handle the requirements for positive weights
            if (poswts_)
                f(f.size() - 2) = (wts.array() < 0).select(wts, 0).sum();

            // Account for invalid arguments
            f(f.size() - 1) = (x - dom_.clamp_args(x)).norm();

            return (f.cwiseAbs().maxCoeff() > 100) ? -1 : 0;
        }

        Derived& dom_;
        const bool poswts_;
    };

    Derived& derived = static_cast<Derived&>(*this);

    typedef Eigen::NumericalDiff<min_functor> min_functor_ndiff;
    min_functor f(derived, poswts_);
    min_functor_ndiff fd(f);

    // Perform the minimisation
    Eigen::LevenbergMarquardt<min_functor_ndiff> lm(fd);
    lm.setMaxfev(maxfev);
    lm.minimize(args_);

    // Clamp the arguments to ensure all points are inside the domain
    args_ = derived.clamp_args(args_);

    // Sort these clamped arguments into a canonical order
    derived.sort_args(orbits_, args_);

    // Compute the residual of these clamped points
    VectorXT resid(derived.nbfn());
    auto rwts = derived.wts(args_, &resid);

    // Compute the residual norm
    T rnorm = resid.norm();

    // If requested, penalize negative weights
    if (poswts_)
    {
        T pen = (rwts.array() < 0).select(rwts, 0).sum();
        rnorm = sqrt(rnorm*rnorm + pen*pen);
    }

    // Return
    return std::make_tuple(rnorm, args_);
}

template<typename Derived, typename T, int Ndim, int Norbits>
inline auto
BaseDomain<Derived, T, Ndim, Norbits>::wts(
        const VectorXT& args, VectorXT* resid, VectorXT* cnorms) -> VectorXT
{
    const Derived& derived = static_cast<const Derived&>(*this);

    derived.expand(args, pts_);
    derived.eval_orthob(pts_, obat_);

    for (int i = 0, poff = 0, coff = 0; i < Norbits; ++i)
    {
        int pinc = derived.npts_for_orbit[i];

        for (int j = 0; j < orbits_(i); ++j, ++coff, poff += pinc)
            for (int l = 0; l < A_.rows(); ++l)
                A_(l, coff) = obat_.block(l, poff, 1, pinc).sum();
    }

    // Compute the optimal set of weights
    wts_.noalias() = A_.householderQr().solve(b_);

    if (resid)
    {
        (*resid).head(b_.size()).noalias() = A_*wts_ - b_;
    }

    if (cnorms)
        *cnorms = A_.colwise().squaredNorm();

    return wts_;
}

template<typename Derived, typename T, int Ndim, int Norbits>
inline int
BaseDomain<Derived, T, Ndim, Norbits>::ndof() const
{
    const Derived& derived = static_cast<const Derived&>(*this);
    int s = 0;

    for (int i = 0; i < Norbits; ++i)
        s += orbits_(i)*derived.narg_for_orbit[i];

    return s;
}

template<typename Derived, typename T, int Ndim, int Norbits>
template<std::size_t N>
inline double
BaseDomain<Derived, T, Ndim, Norbits>::rand(
    double a,
    double b,
    const std::array<int, N>& wts)
{
    const int sum = std::accumulate(std::cbegin(wts), std::cend(wts), 0);
    const double step = (b - a) / N;

    int bin = std::uniform_int_distribution(0, sum - 1)(rand_eng_);
    for (int i = 0; i < N; bin -= wts[i++])
        if (bin < wts[i])
            return rand(a + i*step, a + (i + 1)*step);

    abort();
}

template<typename Derived, typename T, int Ndim, int Norbits>
inline void
BaseDomain<Derived, T, Ndim, Norbits>::sort_args(const VectorOrb& orb, VectorXT& args)
{
    for (int i = 0; i < Norbits; ++i)
    {
        int nobi = orb(i);
        int narg = Derived::narg_for_orbit[i];
        int aoff = arg_offset(orb, i);

        // Skip orbits which have no arguments/are not present
        if (!narg || !nobi)
            continue;

        // Sort the parameters inside of each orbit
        for (int j = 0; j < nobi; ++j)
            Derived::sort_arg(i, aoff + j*narg, args);

        // Map the arguments for the orbits of this type to a matrix
        Eigen::Map<MatrixXT> margs(args.data() + aoff, narg, nobi);

        // Now, reorder the orbits themselves
        for (int j = 0; j < nobi - 1; ++j)
        {
            int mix = j;

            for (int k = j + 1; k < nobi; ++k)
            {
                auto ck = margs.col(k).data(), cm = margs.col(mix).data();

                if (std::lexicographical_compare(ck, ck + narg, cm, cm + narg))
                    mix = k;
            }

            margs.col(mix).swap(margs.col(j));
        }
    }
}

template<typename Derived, typename T, int Ndim, int Norbits>
inline void
BaseDomain<Derived, T, Ndim, Norbits>::symm_decomps_recurse(
        VectorOrb coeffs,
        const VectorOrb& limits,
        int sum,
        VectorOrb& partsoln,
        std::vector<VectorOrb>& solns)
{
    int index;
    int mcoeff = coeffs.maxCoeff(&index);
    int range = std::min(sum / mcoeff, limits(index));

    if (range*mcoeff == sum)
    {
        partsoln(index) = range--;

        if (Derived::validate_orbit(partsoln))
            solns.push_back(partsoln);
    }

    if (coeffs.count() == 1)
        return;

    while (range >= 0)
    {
        int rem = sum - range*mcoeff;

        VectorOrb coeffCopy = coeffs;
        coeffCopy(index) = 0;

        partsoln(index) = range--;

        for (int i = 0; i < coeffCopy.rows(); ++i)
            if (coeffCopy(i) != 0)
                partsoln(i) = 0;

        symm_decomps_recurse(coeffCopy, limits, rem, partsoln, solns);
    }
}

template<typename Derived, typename T, int Ndim, int Norbits>
inline auto
BaseDomain<Derived, T, Ndim, Norbits>::remove_orbit(
    const VectorOrb& orb, const VectorXT& args,
    int i, int j) -> VectorOrbArgs
{
    VectorOrb norb = orb; --norb(i);

    int narg = Derived::narg_for_orbit[i];
    int nh = arg_offset(orb, i, j), nt = args.size() - nh - narg;
    VectorXT nargs(nh + nt);

    nargs.head(nh) = args.head(nh);
    nargs.tail(nt) = args.tail(nt);

    return {norb, nargs};
}

template<typename Derived, typename T, int Ndim, int Norbits>
inline auto
BaseDomain<Derived, T, Ndim, Norbits>::add_orbit(
    const VectorOrb& orb, const VectorXT& args,
    int i, const std::initializer_list<T> aiargs) -> VectorOrbArgs
{
    VectorOrb norb = orb; ++norb(i);

    int nh = arg_offset(orb, i), nt = args.size() - nh;
    VectorXT nargs(args.size() + aiargs.size());

    nargs.head(nh) = args.head(nh);
    std::copy(std::cbegin(aiargs), std::cend(aiargs), nargs.data() + nh);
    nargs.tail(nt) = args.tail(nt);

    sort_args(norb, nargs);
    return {norb, nargs};
}

template<typename Derived, typename T, int Ndim, int Norbits>
inline auto
BaseDomain<Derived, T, Ndim, Norbits>::possible_reductions(
    const VectorXT& args,
    const T& tol) -> std::vector<VectorOrbArgs>
{
    Derived& derived = static_cast<Derived&>(*this);
    int tnpts = npts();

    VectorXT cnorms(nwts());
    const auto wts = derived.wts(args, nullptr, &cnorms);

    for (int i = 0, woff = 0; i < Norbits; ++i)
        for (int j = 0; j < orbits_(i); ++j, ++woff)
            cnorms(woff) *= Derived::npts_for_orbit[i]*abs(wts(woff));

    int mcnormix;
    cnorms.minCoeff(&mcnormix);

    std::vector<VectorOrbArgs> reductions;

    // Loop over each orbit type
    for (int i = 0, aoff = 0, woff = 0; i < Norbits; ++i)
    {
        int nobi = orbits_(i);
        int narg = Derived::narg_for_orbit[i];

        // Check for small weights and unimportant orbits
        for (int j = 0; j < nobi; ++j, ++woff)
            if (woff == mcnormix)
                reductions.push_back(remove_orbit(orbits_, args, i, j));

        if (!nobi || !narg)
            continue;

        Eigen::Map<const MatrixXT> margs(args.data() + arg_offset(i), narg, nobi);

        for (int j = 0; j < nobi; ++j, aoff += narg)
        {
            // See if the orbit is similar to any other orbits
            for (int k = j + 1; k < nobi; ++k)
                if ((margs.col(j) - margs.col(k)).norm() < tol)
                    reductions.push_back(remove_orbit(orbits_, args, i, k));

            auto replace = [&](int p, auto... param)
            {
                auto [rorb, rargs] = remove_orbit(orbits_, args, i, j);
                auto [aorb, aargs] = add_orbit(rorb, rargs, p, {param...});

                if (Derived::validate_orbit(aorb))
                    reductions.push_back({aorb, aargs});
                else
                    reductions.push_back({rorb, rargs});
            };

            // See if the orbit exhibits degeneracy
            Derived::collapse_arg(i, aoff, args, replace, tol);
        }
    }

    std::sort(std::begin(reductions), std::end(reductions),
              [](const auto& p, const auto& q)
              { return npts(p.first) < npts(q.first); });
    return reductions;
}

}

#endif /* POLYQUAD_SHAPES_BASE_HPP */
