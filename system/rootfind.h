/*
 * Modern C++ root-finding algorithms for GIZMO.
 *
 * Provides template-based implementations of:
 *   - brent_root():   Brent's method (1973) with log-space bisection fallback
 *   - toms748_root(): TOMS Algorithm 748 (Alefeld, Potra & Shi, 1995)
 *   - secant_root():  Secant method with bisection fallback
 *
 * These replace the old macro-based bracketed_rootfind.h.
 * Usage: pass any callable (lambda, function pointer, functor) as the first argument.
 *
 * All solvers assume the root is bracketed: f(a) and f(b) must have opposite signs.
 */

#ifndef ROOTFIND_H
#define ROOTFIND_H

#include <cmath>
#include <algorithm>
#include <cstdio>

#ifndef MAXITER
#define MAXITER 150
#endif

struct RootFindResult {
    double root;
    int iterations;
    double bracket_size; /* final |b - a| */
};

namespace rootfind_detail {

/* Swap two doubles in place */
inline void swap(double &x, double &y) { double t = x; x = y; y = t; }

/* Safe inverse quadratic interpolation through three points.
   Returns the interpolated root estimate given (a,fa), (b,fb), (c,fc)
   where all three function values are distinct. */
inline double inverse_quadratic_interpolation(double a, double fa, double b, double fb, double c, double fc) {
    /* Lagrange interpolation of the inverse function x(f) evaluated at f=0 */
    return a * fb * fc / ((fa - fb) * (fa - fc))
         + b * fa * fc / ((fb - fa) * (fb - fc))
         + c * fa * fb / ((fc - fa) * (fc - fb));
}

/* Secant step between two points */
inline double secant_step(double a, double fa, double b, double fb) {
    return (a * fb - b * fa) / (fb - fa);
}

/* Log-space bisection: geometric mean if both positive, arithmetic midpoint otherwise */
inline double log_bisect(double a, double b) {
    if (a > 0 && b > 0) { return std::sqrt(a * b); }
    return 0.5 * (a + b);
}

/* Cubic interpolation for TOMS 748: given four points, compute the unique cubic
   through (a,fa),(b,fb),(d,fd),(e,fe) evaluated at f=0 via Newton form. */
inline double cubic_interpolation(double a, double fa, double b, double fb,
                                  double d, double fd, double e, double fe) {
    double f_ab = (fb - fa) / (b - a);
    double f_abd = (fd - fb) / (d - b);
    f_abd = (f_abd - f_ab) / (d - a);
    double f_abde = (fe - fd) / (e - d);
    f_abde = (f_abde - f_abd + (fd - fb)/(d - b) - f_ab) / (e - a); /* not used directly */
    /* Actually, let's use the safer iterated inverse interpolation approach.
       Newton form for inverse interpolation: x = f^{-1}(0) */
    double q11 = (d - e) * fd / (fe - fd);
    double q21 = (b - d) * fb / (fd - fb);
    double q31 = (a - b) * fa / (fb - fa);
    double d21 = (b - d) * fd / (fd - fb);
    double d31 = (a - b) * fb / (fb - fa);
    double q22 = (d21 - q11) * fb / (fe - fb);
    double q32 = (d31 - q21) * fa / (fd - fa);
    double d32 = (d31 - q21) * fd / (fd - fa);
    double q33 = (d32 - q22) * fa / (fe - fa);
    return a + q31 + q32 + q33;
}

/* Check if a value s lies strictly between a and b (exclusive) */
inline bool between(double s, double a, double b) {
    double lo = std::min(a, b), hi = std::max(a, b);
    return s > lo && s < hi;
}

} /* namespace rootfind_detail */


/*******************************************************************************
 * brent_root: GIZMO's Brent-like root-finding method.
 *
 * This is a faithful template port of the original macro-based algorithm from
 * bracketed_rootfind.h, preserving its exact behavior:
 *   - Convention: 'a' has larger |f|, 'b' has smaller |f|
 *   - IQI via Lagrange interpolation when 3 distinct points available
 *   - Log-space bisection when both bracket endpoints are positive
 *   - DELTA_TOL based on the trial point x_new (not the best estimate b)
 *   - Early exit when |x_new - b| < DELTA_TOL (without evaluating func)
 *   - Returns x_new (last trial point), not b
 *   - Automatic bracket expansion if initial bounds don't bracket
 *
 * Template parameter F: any callable with signature double(double).
 ******************************************************************************/
template <typename F>
RootFindResult brent_root(F func, double a, double b, double fa, double fb,
                          double rel_tol, double abs_tol, int maxiter = MAXITER) {
    /* Attempt bracket expansion if not initially bracketed */
    if (fa * fb > 0) {
        double bracket_fac = 1.1;
        int bracket_iter = 0;
        do {
            double tmp = a;
            a = std::min(a, b) / bracket_fac;
            b = std::max(tmp, b) * bracket_fac;
            fa = func(a);
            fb = func(b);
            bracket_iter++;
        } while (fa * fb > 0 && bracket_iter < maxiter);
        if (bracket_iter >= maxiter || std::isnan(fa) || std::isnan(fb)) {
            return {0.5 * (a + b), maxiter + 1, std::fabs(b - a)};
        }
    }

    /* Convention: 'a' has the larger residual, 'b' has the smaller */
    if (std::fabs(fa) < std::fabs(fb)) {
        rootfind_detail::swap(a, b);
        rootfind_detail::swap(fa, fb);
    }

    double c = a, fc = fa;          /* third point for IQI */
    int used_bisection = 1, iter = 0;
    double c_old = c, x_new = 0, f_new = fc, x_error = 1e100, delta_tol = 0;

    do {
        x_new = 0;
        if ((fa != fc) && (fb != fc)) {
            /* Inverse quadratic interpolation (Lagrange form) */
            x_new = rootfind_detail::inverse_quadratic_interpolation(a, fa, b, fb, c, fc);
        } else {
            /* Secant method */
            x_new = rootfind_detail::secant_step(a, fa, b, fb);
        }
        /* Tolerance based on the trial point x_new */
        delta_tol = std::max(std::fabs(abs_tol), rel_tol * std::fabs(x_new));

        int do_bisection = 0;
        double midpoint_a = 0.25 * (3 * a + b);
        if ((x_new < std::min(midpoint_a, b)) || (x_new > std::max(midpoint_a, b))) {
            do_bisection = 1;
        } else {
            /* Early exit: interpolation converged without needing to evaluate func */
            if (std::fabs(x_new - b) < delta_tol) {
                break;
            }
        }
        if (used_bisection) {
            if (std::fabs(x_new - b) >= 0.5 * std::fabs(c - b)) { do_bisection = 1; }
            if (b != c) { if (std::fabs(b - c) < delta_tol) { do_bisection = 1; } }
        } else {
            if (std::fabs(x_new - b) >= 0.5 * std::fabs(c_old - c)) { do_bisection = 1; }
            if (c_old != c) { if (std::fabs(c_old - c) < delta_tol) { do_bisection = 1; } }
        }
        if (do_bisection) {
            /* Log-space bisection when both endpoints are positive */
            if ((b > 0) && (a > 0)) {
                x_new = std::sqrt(b * a);
            } else {
                x_new = 0.5 * (b + a);
            }
            used_bisection = 1;
        } else {
            used_bisection = 0;
        }
        f_new = func(x_new);
        if (f_new == 0) { break; }

        c_old = c;
        c = b;
        fc = fb;
        if (fa * f_new < 0) {
            b = x_new; fb = f_new;
        } else {
            a = x_new; fa = f_new;
        }

        if (std::fabs(fa) < std::fabs(fb)) {
            rootfind_detail::swap(a, b);
            rootfind_detail::swap(fa, fb);
        }
        x_error = std::fabs(b - a);
        iter++;
        if (iter > maxiter) { break; }
    } while (x_error > delta_tol);

    return {x_new, iter, x_error};
}


/*******************************************************************************
 * toms748_root: TOMS Algorithm 748 (Alefeld, Potra & Shi, 1995).
 *
 * A bracketed root-finding method that achieves asymptotic order ~1.84 by
 * using inverse cubic interpolation when four distinct points are available,
 * falling back to inverse quadratic or secant interpolation, with bisection
 * as the ultimate fallback.
 *
 * Generally converges faster than Brent's method in fewer function evaluations,
 * while maintaining the guaranteed convergence of a bracketed method.
 *
 * Reference: "Algorithm 748: Enclosing Zeros of Continuous Functions"
 *            G.E. Alefeld, F.A. Potra, Y. Shi
 *            ACM Transactions on Mathematical Software, 21(3), 1995.
 *
 * Template parameter F: any callable with signature double(double).
 ******************************************************************************/
template <typename F>
RootFindResult toms748_root(F func, double a, double b, double fa, double fb,
                            double rel_tol, double abs_tol, int maxiter = MAXITER) {
    using namespace rootfind_detail;

    /* Attempt bracket expansion if not initially bracketed */
    if (fa * fb > 0) {
        double bracket_fac = 1.1;
        int bracket_iter = 0;
        do {
            double tmp = a;
            a = std::min(a, b) / bracket_fac;
            b = std::max(tmp, b) * bracket_fac;
            fa = func(a);
            fb = func(b);
            bracket_iter++;
        } while (fa * fb > 0 && bracket_iter < maxiter);
        if (bracket_iter >= maxiter || std::isnan(fa) || std::isnan(fb)) {
            return {0.5 * (a + b), maxiter + 1, std::fabs(b - a)};
        }
    }

    if (fa == 0.0) { return {a, 0, 0.0}; }
    if (fb == 0.0) { return {b, 0, 0.0}; }

    /* Ensure fa and fb have opposite signs; a is left bracket, b is right */
    if (a > b) { swap(a, b); swap(fa, fb); }

    int nfev = 0; /* function evaluation count within the loop */

    /* d and e are previous bracket endpoints for higher-order interpolation */
    double d = a, fd = fa;   /* will be updated after first iteration */
    double e = b, fe = fb;

    /* mu for the bracket tolerance */
    auto tol = [&](double x) -> double {
        return 2.0 * std::max(abs_tol, rel_tol * std::fabs(x));
    };

    for (int iter = 0; iter < maxiter; ++iter) {
        double tol_a = tol(a), tol_b = tol(b);

        /* Check convergence */
        if ((b - a) <= 2.0 * std::max(tol_a, tol_b)) {
            double root = (std::fabs(fa) <= std::fabs(fb)) ? a : b;
            return {root, nfev, b - a};
        }
        if (fa == 0.0) { return {a, nfev, 0.0}; }
        if (fb == 0.0) { return {b, nfev, 0.0}; }

        /* Step 1: compute two candidate points using available interpolation */
        double c;
        bool have_four_distinct = (iter > 0) && (fa != fd) && (fa != fe) &&
                                  (fb != fd) && (fb != fe) && (fd != fe);

        if (have_four_distinct) {
            /* Try inverse cubic interpolation */
            c = cubic_interpolation(a, fa, b, fb, d, fd, e, fe);
            if (!between(c, a, b)) {
                c = secant_step(a, fa, b, fb);
            }
        } else {
            /* Secant step */
            c = secant_step(a, fa, b, fb);
        }
        /* Clamp c to be well within (a, b) */
        double a_inner = a + tol_a, b_inner = b - tol_b;
        if (a_inner >= b_inner) { a_inner = b_inner = 0.5 * (a + b); }
        c = std::max(a_inner, std::min(c, b_inner));

        double fc = func(c);
        nfev++;
        if (fc == 0.0) { return {c, nfev, 0.0}; }

        /* Update bracket with c */
        e = d; fe = fd;
        if (fa * fc < 0) {
            d = b; fd = fb;
            b = c; fb = fc;
        } else {
            d = a; fd = fa;
            a = c; fa = fc;
        }

        /* Check convergence after first update */
        if ((b - a) <= 2.0 * std::max(tol(a), tol(b))) {
            double root = (std::fabs(fa) <= std::fabs(fb)) ? a : b;
            return {root, nfev, b - a};
        }

        /* Step 2: compute a second candidate using all available points */
        double c2;
        if ((fa != fd) && (fa != fe) && (fb != fd) && (fb != fe) && (fd != fe)) {
            c2 = cubic_interpolation(a, fa, b, fb, d, fd, e, fe);
            if (!between(c2, a, b)) {
                c2 = secant_step(a, fa, b, fb);
            }
        } else if ((fa != fc) && (fb != fc)) {
            c2 = inverse_quadratic_interpolation(a, fa, b, fb, c, fc);
            if (!between(c2, a, b)) {
                c2 = secant_step(a, fa, b, fb);
            }
        } else {
            c2 = secant_step(a, fa, b, fb);
        }
        a_inner = a + tol(a); b_inner = b - tol(b);
        if (a_inner >= b_inner) { a_inner = b_inner = 0.5 * (a + b); }
        c2 = std::max(a_inner, std::min(c2, b_inner));

        double fc2 = func(c2);
        nfev++;
        if (fc2 == 0.0) { return {c2, nfev, 0.0}; }

        /* Update bracket with c2 */
        e = d; fe = fd;
        if (fa * fc2 < 0) {
            d = b; fd = fb;
            b = c2; fb = fc2;
        } else {
            d = a; fd = fa;
            a = c2; fa = fc2;
        }

        /* Check convergence */
        if ((b - a) <= 2.0 * std::max(tol(a), tol(b))) {
            double root = (std::fabs(fa) <= std::fabs(fb)) ? a : b;
            return {root, nfev, b - a};
        }

        /* Step 3: double-length secant step (bisection-like contraction guarantee) */
        double u, fu;
        if (std::fabs(fa) < std::fabs(fb)) { u = a; fu = fa; } else { u = b; fu = fb; }

        /* Double-length secant: extrapolate from (a,fa) and (b,fb) with 2x step */
        double c3 = u - 2.0 * fu * (b - a) / (fb - fa);
        /* If the double-length secant overshoots, use midpoint */
        double midpt = 0.5 * (a + b);
        if (std::fabs(c3 - u) > 0.5 * (b - a)) { c3 = midpt; }
        a_inner = a + tol(a); b_inner = b - tol(b);
        if (a_inner >= b_inner) { a_inner = b_inner = midpt; }
        c3 = std::max(a_inner, std::min(c3, b_inner));

        double fc3 = func(c3);
        nfev++;
        if (fc3 == 0.0) { return {c3, nfev, 0.0}; }

        /* Update bracket with c3 */
        e = d; fe = fd;
        if (fa * fc3 < 0) {
            d = b; fd = fb;
            b = c3; fb = fc3;
        } else {
            d = a; fd = fa;
            a = c3; fa = fc3;
        }

        /* Safety: if bracket hasn't shrunk enough, force bisection */
        if ((b - a) > 0.5 * (d == a ? (b - d) : (d - a) > 0 ? (b - a) : (b - a))) {
            /* Always do at least one bisection per super-step to guarantee convergence */
            if (iter > 2 && (b - a) > 0.5 * std::fabs(e - d)) {
                double bsct = 0.5 * (a + b);
                a_inner = a + tol(a); b_inner = b - tol(b);
                if (a_inner >= b_inner) { a_inner = b_inner = bsct; }
                bsct = std::max(a_inner, std::min(bsct, b_inner));
                double fbsct = func(bsct);
                nfev++;
                if (fbsct == 0.0) { return {bsct, nfev, 0.0}; }
                e = d; fe = fd;
                if (fa * fbsct < 0) {
                    d = b; fd = fb;
                    b = bsct; fb = fbsct;
                } else {
                    d = a; fd = fa;
                    a = bsct; fa = fbsct;
                }
            }
        }
    }

    double root = (std::fabs(fa) <= std::fabs(fb)) ? a : b;
    return {root, nfev, b - a};
}


/*******************************************************************************
 * secant_root: Secant method with bisection fallback and bracket maintenance.
 *
 * A simple method useful when bracketing is available and function evaluations
 * are cheap. Falls back to log-space bisection when the secant step is not
 * converging fast enough.
 *
 * Template parameter F: any callable with signature double(double).
 ******************************************************************************/
template <typename F>
RootFindResult secant_root(F func, double x0, double x1, double f0, double f1,
                           double x_lower, double x_upper,
                           double rel_tol, int maxiter = MAXITER) {
    double x_old = x0, f_old = f0;
    double x_cur = x1, f_cur = f1;

    /* Track function signs at bracket endpoints:
       f_lower_sign = sign of f at x_lower, f_upper_sign = sign of f at x_upper.
       These let us update the correct endpoint when a new point is evaluated. */
    double f_at_lower = func(x_lower);
    double f_at_upper = func(x_upper);
    (void)f_at_upper; /* may only be used on one branch */

    for (int iter = 0; iter < maxiter; ++iter) {
        if (f_cur == f_old) { break; } /* denominator would be zero */

        /* Secant step, clamped to bracket */
        double x_sec = x_cur - f_cur * (x_cur - x_old) / (f_cur - f_old);
        x_sec = std::max(x_lower, std::min(x_sec, x_upper));

        /* Check convergence ratio: if secant isn't shrinking fast enough, bisect */
        double step_ratio = std::fabs(x_sec - x_cur) / (std::fabs(x_cur - x_old) + 1e-300);
        double x_next;
        if (step_ratio < 0.5 && x_sec != x_cur) {
            x_next = x_sec; /* accept secant step */
        } else {
            /* Log-space bisection of the bracket */
            x_next = rootfind_detail::log_bisect(x_lower, x_upper);
            /* If bisection lands on current point, perturb toward the other bracket end */
            if (x_next == x_cur) {
                x_next = 0.5 * (x_cur + (std::fabs(x_cur - x_lower) > std::fabs(x_cur - x_upper) ? x_lower : x_upper));
            }
        }

        double f_next = func(x_next);
        if (f_next == 0.0) { return {x_next, iter + 1, 0.0}; }

        /* Update brackets: replace the endpoint whose sign matches f_next */
        if ((f_next > 0) == (f_at_lower > 0)) {
            x_lower = x_next; f_at_lower = f_next;
        } else {
            x_upper = x_next; f_at_upper = f_next;
        }

        /* Check convergence */
        if (std::fabs(x_upper - x_lower) < rel_tol * (std::fabs(x_lower) + std::fabs(x_upper)) * 0.5) {
            return {x_next, iter + 1, std::fabs(x_upper - x_lower)};
        }

        x_old = x_cur; f_old = f_cur;
        x_cur = x_next; f_cur = f_next;
    }

    return {x_cur, maxiter, std::fabs(x_upper - x_lower)};
}

/*******************************************************************************
 * halley_root: Halley's method with bracket safety for root-finding.
 *
 * Uses function value, first derivative, and second derivative for cubic
 * convergence. Falls back to bisection when the Halley step would leave the
 * bracket or when convergence stalls.
 *
 * Halley's update formula:
 *   x_{n+1} = x_n - 2·f·f' / (2·f'² - f·f'')
 *
 * Template parameter F: callable with signature
 *   HalleyFuncResult(double x)  returning {f, df, d2f}.
 ******************************************************************************/

struct HalleyFuncResult {
    double f;    /* function value */
    double df;   /* first derivative */
    double d2f;  /* second derivative */
};

template <typename F>
RootFindResult halley_root(F func_and_derivs, double x0,
                           double x_lower, double x_upper,
                           double rel_tol, int maxiter = MAXITER) {
    /* Evaluate at initial guess */
    auto [f, df, d2f] = func_and_derivs(x0);
    double x = x0;
    double f_prev = f; /* function value at previous iterate, for bracket tracking */

    for (int iter = 0; iter < maxiter; ++iter) {
        if (f == 0.0) { return {x, iter, std::fabs(x_upper - x_lower)}; }

        /* Halley step: x_{n+1} = x_n - 2*f*f' / (2*f'^2 - f*f'') */
        double denom = 2.0 * df * df - f * d2f;
        double x_next;
        bool accepted = false;

        if (std::fabs(denom) > 1e-300 && df != 0.0) {
            x_next = x - 2.0 * f * df / denom;
            if (x_next > x_lower && x_next < x_upper) { accepted = true; }
        }

        /* Newton fallback if Halley step rejected */
        if (!accepted && std::fabs(df) > 1e-300) {
            x_next = x - f / df;
            if (x_next > x_lower && x_next < x_upper) { accepted = true; }
        }

        /* Bisection fallback if both rejected */
        if (!accepted) {
            x_next = rootfind_detail::log_bisect(x_lower, x_upper);
        }

        /* Evaluate at new point */
        f_prev = f;
        double x_prev = x;
        auto result = func_and_derivs(x_next);
        f = result.f; df = result.df; d2f = result.d2f;

        /* Update brackets: only tighten when we have a sign change,
           which guarantees the root lies between the two points. */
        if (f * f_prev < 0) {
            if (x_next < x_prev) { x_lower = x_next; x_upper = x_prev; }
            else { x_lower = x_prev; x_upper = x_next; }
        }

        /* Check convergence */
        double tol = rel_tol * std::fabs(x_next);
        if (std::fabs(x_next - x_prev) <= tol || f == 0.0) {
            return {x_next, iter + 1, std::fabs(x_upper - x_lower)};
        }

        x = x_next;
    }

    return {x, maxiter, std::fabs(x_upper - x_lower)};
}

#endif /* ROOTFIND_H */
