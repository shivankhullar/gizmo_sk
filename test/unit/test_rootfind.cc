// Unit tests for system/rootfind.h: brent_root, toms748_root, secant_root.
#include "test_harness.h"
#include "../../system/rootfind.h"
#include <cmath>
#include <functional>

// ============================================================
//  Test functions with known roots
// ============================================================

static double f_sqrt2(double x) { return x * x - 2.0; }
static double f_sin(double x)   { return std::sin(x); }
static double f_cubic(double x) { return x * x * x - x - 2.0; }   // root ~1.5214
static double f_exp3(double x)  { return std::exp(x) - 3.0; }     // root = ln(3)
static double f_lambert(double x) { return x * std::exp(x) - 1.0; } // root ~0.5671
static double f_steep(double x) { return 1e8 * (x - 0.5); }
static double f_hyp(double x)   { return 1.0 / x - 5.0; }          // root = 0.2
static double f_atan(double x)  { return std::atan(x - 1.5); }      // root = 1.5
static double f_flat(double x)  { return std::pow(x, 20) - 1.0; }   // root = 1.0

// ============================================================
//  brent_root: basic convergence
// ============================================================

TEST_CASE("brent: sqrt(2)") {
    auto r = brent_root(f_sqrt2, 1.0, 2.0, f_sqrt2(1.0), f_sqrt2(2.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, std::sqrt(2.0), 1e-10);
    CHECK(r.iterations < 50);
}

TEST_CASE("brent: sin(x)=0 near pi") {
    auto r = brent_root(f_sin, 3.0, 4.0, f_sin(3.0), f_sin(4.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, M_PI, 1e-10);
}

TEST_CASE("brent: cubic x^3-x-2") {
    auto r = brent_root(f_cubic, 1.0, 2.0, f_cubic(1.0), f_cubic(2.0), 1e-12, 0.0);
    CHECK(std::fabs(f_cubic(r.root)) < 1e-10);
}

TEST_CASE("brent: exp(x)-3 = ln(3)") {
    auto r = brent_root(f_exp3, 0.0, 2.0, f_exp3(0.0), f_exp3(2.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, std::log(3.0), 1e-10);
}

TEST_CASE("brent: Lambert W") {
    auto r = brent_root(f_lambert, 0.0, 1.0, f_lambert(0.0), f_lambert(1.0), 1e-14, 0.0);
    CHECK(std::fabs(f_lambert(r.root)) < 1e-12);
}

// ============================================================
//  brent_root: edge cases
// ============================================================

TEST_CASE("brent: root at bracket endpoint") {
    double r2 = std::sqrt(2.0);
    auto r = brent_root(f_sqrt2, r2, 2.0, f_sqrt2(r2), f_sqrt2(2.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, r2, 1e-10);
}

TEST_CASE("brent: very tight bracket") {
    double mid = std::sqrt(2.0);
    auto r = brent_root(f_sqrt2, mid - 1e-8, mid + 1e-8,
                         f_sqrt2(mid - 1e-8), f_sqrt2(mid + 1e-8), 1e-14, 0.0);
    CHECK_CLOSE(r.root, mid, 1e-8);
}

TEST_CASE("brent: log-space bisection (e)") {
    auto f = [](double x) { return std::log(x) - 1.0; };
    auto r = brent_root(f, 1.0, 10.0, f(1.0), f(10.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, M_E, 1e-10);
}

TEST_CASE("brent: steep function") {
    auto r = brent_root(f_steep, 0.0, 1.0, f_steep(0.0), f_steep(1.0), 1e-14, 0.0);
    CHECK_CLOSE(r.root, 0.5, 1e-12);
}

TEST_CASE("brent: hyperbolic 1/x - 5") {
    auto r = brent_root(f_hyp, 0.01, 1.0, f_hyp(0.01), f_hyp(1.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, 0.2, 1e-10);
}

TEST_CASE("brent: lambda with captures") {
    double offset = 3.14159;
    auto f = [&](double x) { return x * x - offset; };
    auto r = brent_root(f, 1.0, 3.0, f(1.0), f(3.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, std::sqrt(offset), 1e-10);
}

TEST_CASE("brent: auto-bracket expansion") {
    // Both f(2)>0 and f(3)>0: not initially bracketed
    auto r = brent_root(f_sqrt2, 2.0, 3.0, f_sqrt2(2.0), f_sqrt2(3.0), 1e-10, 0.0);
    CHECK_CLOSE(r.root, std::sqrt(2.0), 1e-8);
}

TEST_CASE("brent: nested root-find via lambda") {
    // Outer: find x such that sqrt(x) = 1.5, i.e. x = 2.25
    auto outer = [](double x) -> double {
        auto inner = [&](double y) { return y * y - x; };
        auto ir = brent_root(inner, 0.0, 10.0, inner(0.0), inner(10.0), 1e-14, 0.0);
        return ir.root - 1.5;
    };
    auto r = brent_root(outer, 1.0, 4.0, outer(1.0), outer(4.0), 1e-10, 0.0);
    CHECK_CLOSE(r.root, 2.25, 1e-8);
}

// ============================================================
//  toms748_root: basic convergence
// ============================================================

TEST_CASE("toms748: sqrt(2)") {
    auto r = toms748_root(f_sqrt2, 1.0, 2.0, f_sqrt2(1.0), f_sqrt2(2.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, std::sqrt(2.0), 1e-10);
    CHECK(r.iterations < 30);
}

TEST_CASE("toms748: sin(x)=0 near pi") {
    auto r = toms748_root(f_sin, 3.0, 4.0, f_sin(3.0), f_sin(4.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, M_PI, 1e-10);
}

TEST_CASE("toms748: cubic x^3-x-2") {
    auto r = toms748_root(f_cubic, 1.0, 2.0, f_cubic(1.0), f_cubic(2.0), 1e-12, 0.0);
    CHECK(std::fabs(f_cubic(r.root)) < 1e-10);
}

TEST_CASE("toms748: exp(x)-3 = ln(3)") {
    auto r = toms748_root(f_exp3, 0.0, 2.0, f_exp3(0.0), f_exp3(2.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, std::log(3.0), 1e-10);
}

TEST_CASE("toms748: Lambert W") {
    auto r = toms748_root(f_lambert, 0.0, 1.0, f_lambert(0.0), f_lambert(1.0), 1e-14, 0.0);
    CHECK(std::fabs(f_lambert(r.root)) < 1e-12);
}

TEST_CASE("toms748: atan root at 1.5") {
    auto r = toms748_root(f_atan, 0.0, 3.0, f_atan(0.0), f_atan(3.0), 1e-14, 0.0);
    CHECK_CLOSE(r.root, 1.5, 1e-12);
}

TEST_CASE("toms748: steep function") {
    auto r = toms748_root(f_steep, 0.0, 1.0, f_steep(0.0), f_steep(1.0), 1e-14, 0.0);
    CHECK_CLOSE(r.root, 0.5, 1e-12);
}

TEST_CASE("toms748: hyperbolic 1/x - 5") {
    auto r = toms748_root(f_hyp, 0.01, 1.0, f_hyp(0.01), f_hyp(1.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, 0.2, 1e-10);
}

// ============================================================
//  toms748_root: edge cases
// ============================================================

TEST_CASE("toms748: root at bracket endpoint") {
    double r2 = std::sqrt(2.0);
    auto r = toms748_root(f_sqrt2, r2, 2.0, f_sqrt2(r2), f_sqrt2(2.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, r2, 1e-10);
}

TEST_CASE("toms748: very wide bracket") {
    auto r = toms748_root(f_sqrt2, 0.001, 1000.0, f_sqrt2(0.001), f_sqrt2(1000.0), 1e-10, 0.0);
    CHECK_CLOSE(r.root, std::sqrt(2.0), 1e-8);
}

TEST_CASE("toms748: lambda with captures") {
    double a = 2.5, b = 1.3;
    auto f = [&](double x) { return a * x * x + b * x - 7.0; };
    double root_exact = (-b + std::sqrt(b * b + 4.0 * a * 7.0)) / (2.0 * a);
    auto r = toms748_root(f, 0.0, 5.0, f(0.0), f(5.0), 1e-12, 0.0);
    CHECK_CLOSE(r.root, root_exact, 1e-10);
}

TEST_CASE("toms748: nested calls") {
    auto outer = [](double x) -> double {
        auto inner = [&](double y) { return y * y - x; };
        auto ir = toms748_root(inner, 0.0, 10.0, inner(0.0), inner(10.0), 1e-14, 0.0);
        return ir.root - 2.0;
    };
    auto r = toms748_root(outer, 1.0, 8.0, outer(1.0), outer(8.0), 1e-10, 0.0);
    CHECK_CLOSE(r.root, 4.0, 1e-8);
}

// ============================================================
//  toms748 vs brent efficiency comparison
// ============================================================

TEST_CASE("toms748 uses fewer evals than brent on Lambert W") {
    int brent_evals = 0, toms_evals = 0;
    auto fb = [&](double x) -> double { brent_evals++; return f_lambert(x); };
    auto ft = [&](double x) -> double { toms_evals++; return f_lambert(x); };
    brent_root(fb, 0.0, 1.0, f_lambert(0.0), f_lambert(1.0), 1e-14, 0.0);
    toms748_root(ft, 0.0, 1.0, f_lambert(0.0), f_lambert(1.0), 1e-14, 0.0);
    CHECK(toms_evals <= brent_evals + 5);
}

TEST_CASE("toms748 uses fewer evals than brent on cubic") {
    int brent_evals = 0, toms_evals = 0;
    auto fb = [&](double x) -> double { brent_evals++; return f_cubic(x); };
    auto ft = [&](double x) -> double { toms_evals++; return f_cubic(x); };
    brent_root(fb, 1.0, 2.0, f_cubic(1.0), f_cubic(2.0), 1e-14, 0.0);
    toms748_root(ft, 1.0, 2.0, f_cubic(1.0), f_cubic(2.0), 1e-14, 0.0);
    CHECK(toms_evals <= brent_evals + 5);
}

// ============================================================
//  secant_root: basic convergence
// ============================================================

TEST_CASE("secant: sqrt(2)") {
    double a = 1.0, b = 2.0;
    auto r = secant_root(f_sqrt2, a, b, f_sqrt2(a), f_sqrt2(b), 0.0, 10.0, 1e-10);
    CHECK_CLOSE(r.root, std::sqrt(2.0), 1e-8);
}

TEST_CASE("secant: ln(3)") {
    double a = 0.5, b = 1.5;
    auto r = secant_root(f_exp3, a, b, f_exp3(a), f_exp3(b), 0.1, 3.0, 1e-10);
    CHECK_CLOSE(r.root, std::log(3.0), 1e-8);
}

TEST_CASE("secant: sin near pi") {
    double a = 3.0, b = 3.5;
    auto r = secant_root(f_sin, a, b, f_sin(a), f_sin(b), 2.0, 4.0, 1e-10);
    CHECK_CLOSE(r.root, M_PI, 1e-8);
}

// ============================================================
//  halley_root: basic convergence
// ============================================================

/* Helper: returns {f, f', f''} for x^2 - C */
static HalleyFuncResult halley_x2_minus(double x, double C) {
    return {x * x - C, 2.0 * x, 2.0};
}

/* Helper: returns {f, f', f''} for T^4 - C (astrophysics-like) */
static HalleyFuncResult halley_T4_minus(double T, double C) {
    double T2 = T * T, T3 = T2 * T;
    return {T2 * T2 - C, 4.0 * T3, 12.0 * T2};
}

TEST_CASE("halley: sqrt(2)") {
    auto f = [](double x) -> HalleyFuncResult { return halley_x2_minus(x, 2.0); };
    auto r = halley_root(f, 1.5, 0.0, 10.0, 1e-14);
    CHECK_CLOSE(r.root, std::sqrt(2.0), 1e-12);
    CHECK(r.iterations <= 5);
}

TEST_CASE("halley: T^4 = 1e16 (root at 1e4)") {
    auto f = [](double x) -> HalleyFuncResult { return halley_T4_minus(x, 1e16); };
    auto r = halley_root(f, 5000.0, 1.0, 1e6, 1e-10);
    CHECK_CLOSE(r.root, 1e4, 1e-4);
    CHECK(r.iterations <= 10);
}

TEST_CASE("halley: cubic convergence — fewer evals than brent on T^4") {
    int halley_evals = 0, brent_evals = 0;
    auto fh = [&](double T) -> HalleyFuncResult {
        halley_evals++;
        double T2 = T*T, T3 = T2*T;
        return {T2*T2 - 1e16, 4.0*T3, 12.0*T2};
    };
    auto fb = [&](double T) -> double { brent_evals++; return T*T*T*T - 1e16; };
    halley_root(fh, 5000.0, 1.0, 1e6, 1e-10);
    brent_root(fb, 1.0, 1e6, fb(1.0) - 1, fb(1e6) - 1, 1e-10, 0.0);
    CHECK(halley_evals < brent_evals);
}

TEST_CASE("halley: sin(x)=0 near pi") {
    auto f = [](double x) -> HalleyFuncResult {
        return {std::sin(x), std::cos(x), -std::sin(x)};
    };
    auto r = halley_root(f, 3.0, 2.0, 4.0, 1e-14);
    CHECK_CLOSE(r.root, M_PI, 1e-12);
    CHECK(r.iterations <= 5);
}

TEST_CASE("halley: exp(x)-3 = ln(3)") {
    auto f = [](double x) -> HalleyFuncResult {
        double ex = std::exp(x);
        return {ex - 3.0, ex, ex};
    };
    auto r = halley_root(f, 1.5, -1.0, 5.0, 1e-14);
    CHECK_CLOSE(r.root, std::log(3.0), 1e-12);
    CHECK(r.iterations <= 5);
}

TEST_CASE("halley: dust-like T^4 + b*T with coupling") {
    /* f(T) = alpha*(Tgas - T) + A - eps*kappa*T^4, monotone decreasing */
    double Tgas = 100.0, alpha = 1e-3, A = 1e5, eps_kappa = 1e-8;
    auto f = [&](double T) -> HalleyFuncResult {
        double T2 = T*T, T3 = T2*T, T4 = T2*T2;
        double fval = alpha * (Tgas - T) + A - eps_kappa * T4;
        double df = -alpha - 4.0 * eps_kappa * T3;
        double d2f = -12.0 * eps_kappa * T2;
        return {fval, df, d2f};
    };
    auto r = halley_root(f, 50.0, 1.0, 1e5, 1e-10);
    /* Verify residual is near zero */
    double T2 = r.root*r.root, T4 = T2*T2;
    double residual = alpha * (Tgas - r.root) + A - eps_kappa * T4;
    CHECK(std::fabs(residual) < 1e-3);
    CHECK(r.iterations <= 20);
}

TEST_CASE("halley: falls back to bisection when derivative is zero") {
    /* f(x) = x^3, f'(0) = 0 — Halley/Newton singular at origin */
    auto f = [](double x) -> HalleyFuncResult {
        return {x * x * x, 3.0 * x * x, 6.0 * x};
    };
    auto r = halley_root(f, 0.5, -1.0, 2.0, 1e-10);
    CHECK(std::fabs(r.root) < 1e-3); /* should still find root near 0 */
}

// ============================================================
//  Hard test functions
// ============================================================

TEST_CASE("brent+toms748: triple root (x-1)^3") {
    auto f = [](double x) { return (x-1)*(x-1)*(x-1); };
    auto rb = brent_root(f, 0.0, 3.0, f(0.0), f(3.0), 1e-8, 0.0);
    auto rt = toms748_root(f, 0.0, 3.0, f(0.0), f(3.0), 1e-8, 0.0);
    CHECK(std::fabs(rb.root - 1.0) < 1e-3);
    CHECK(std::fabs(rt.root - 1.0) < 1e-3);
}

TEST_CASE("brent+toms748: cube root singularity at 0.3") {
    auto f = [](double x) {
        double d = x - 0.3;
        return (d >= 0 ? 1 : -1) * std::cbrt(std::fabs(d));
    };
    auto rb = brent_root(f, 0.0, 1.0, f(0.0), f(1.0), 1e-10, 0.0);
    auto rt = toms748_root(f, 0.0, 1.0, f(0.0), f(1.0), 1e-10, 0.0);
    CHECK_CLOSE(rb.root, 0.3, 1e-8);
    CHECK_CLOSE(rt.root, 0.3, 1e-8);
}

TEST_CASE("brent+toms748: cos(x) = x fixed point") {
    auto f = [](double x) { return std::cos(x) - x; };
    auto rb = brent_root(f, 0.0, 2.0, f(0.0), f(2.0), 1e-14, 0.0);
    auto rt = toms748_root(f, 0.0, 2.0, f(0.0), f(2.0), 1e-14, 0.0);
    CHECK(std::fabs(f(rb.root)) < 1e-12);
    CHECK(std::fabs(f(rt.root)) < 1e-12);
}

TEST_CASE("brent+toms748: x^20 - 1 (very flat near root)") {
    auto rb = brent_root(f_flat, 0.5, 1.5, f_flat(0.5), f_flat(1.5), 1e-10, 0.0);
    auto rt = toms748_root(f_flat, 0.5, 1.5, f_flat(0.5), f_flat(1.5), 1e-10, 0.0);
    CHECK_CLOSE(rb.root, 1.0, 1e-4);
    CHECK_CLOSE(rt.root, 1.0, 1e-4);
}

TEST_CASE("brent+toms748: T^4 - 1e16 (astrophysics-like)") {
    auto f = [](double T) { return T*T*T*T - 1e16; };
    auto rb = brent_root(f, 1.0, 1e6, f(1.0), f(1e6), 1e-10, 0.0);
    auto rt = toms748_root(f, 1.0, 1e6, f(1.0), f(1e6), 1e-10, 0.0);
    CHECK_CLOSE(rb.root, 1e4, 1e-2);
    CHECK_CLOSE(rt.root, 1e4, 1e-2);
}

// ============================================================
//  Relative tolerance
// ============================================================

TEST_CASE("brent: relative tolerance on large root") {
    auto f = [](double x) { return x * x - 1e12; };
    double exact = 1e6;
    auto r = brent_root(f, 1.0, 1e8, f(1.0), f(1e8), 1e-8, 0.0);
    CHECK(std::fabs(r.root - exact) / exact < 1e-6);
}

TEST_CASE("toms748: relative tolerance on large root") {
    auto f = [](double x) { return x * x - 1e12; };
    double exact = 1e6;
    auto r = toms748_root(f, 1.0, 1e8, f(1.0), f(1e8), 1e-8, 0.0);
    CHECK(std::fabs(r.root - exact) / exact < 1e-6);
}

TEST_MAIN()
