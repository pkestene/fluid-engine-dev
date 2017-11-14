// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include <pch.h>

#include <jet/exp_grid_fluid_solver2.h>

using namespace jet;
using namespace experimental;

#define IX(i, j) ((i) + (N + 2) * (j))
#define FOR_EACH_CELL          \
    for (i = 1; i <= N; i++) { \
        for (j = 1; j <= N; j++) {
#define END_FOR \
    }           \
    }

namespace {

void addSource(size_t N, float* x, float* s, float dt) {
    size_t i, size = (N + 2) * (N + 2);
    for (i = 0; i < size; i++) {
        x[i] += dt * s[i];
    }
}

void setBoundaryCondition(size_t N, size_t b, float* x) {
    size_t i;

    for (i = 1; i <= N; i++) {
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void solveLinearSystem(size_t N, size_t b, float* x, float* x0, float a,
                       float c) {
    size_t i, j, k;

    for (k = 0; k < 20; k++) {
        FOR_EACH_CELL
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                           x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                      c;
        END_FOR
        setBoundaryCondition(N, b, x);
    }
}

void diffuse(size_t N, size_t b, float* x, float* x0, float diff, float dt) {
    float a = dt * diff * N * N;
    solveLinearSystem(N, b, x, x0, a, 1 + 4 * a);
}

void advect(size_t N, size_t b, float* d, float* d0, float* u, float* v,
            float dt) {
    size_t i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * N;
    FOR_EACH_CELL
    x = i - dt0 * u[IX(i, j)];
    y = j - dt0 * v[IX(i, j)];
    if (x < 0.5f) {
        x = 0.5f;
    }
    if (x > N + 0.5f) {
        x = N + 0.5f;
    }
    i0 = (size_t)x;
    i1 = i0 + 1;
    if (y < 0.5f) {
        y = 0.5f;
    }
    if (y > N + 0.5f) {
        y = N + 0.5f;
    }
    j0 = (size_t)y;
    j1 = j0 + 1;
    s1 = x - i0;
    s0 = 1 - s1;
    t1 = y - j0;
    t0 = 1 - t1;
    d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                  s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
    END_FOR
    setBoundaryCondition(N, b, d);
}

void project(size_t N, float* u, float* v, float* p, float* div) {
    size_t i, j;

    FOR_EACH_CELL
    div[IX(i, j)] = -0.5f *
                    (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] -
                     v[IX(i, j - 1)]) /
                    N;
    p[IX(i, j)] = 0;
    END_FOR
    setBoundaryCondition(N, 0, div);
    setBoundaryCondition(N, 0, p);

    solveLinearSystem(N, 0, p, div, 1, 4);

    FOR_EACH_CELL
    u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
    v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
    END_FOR
    setBoundaryCondition(N, 1, u);
    setBoundaryCondition(N, 2, v);
}

void densityStep(size_t N, float* x, float* x0, float* u, float* v, float diff,
                 float dt) {
    addSource(N, x, x0, dt);
    std::swap(x0, x);
    diffuse(N, 0, x, x0, diff, dt);
    std::swap(x0, x);
    advect(N, 0, x, x0, u, v, dt);
}

void velocityStep(size_t N, float* u, float* v, float* u0, float* v0,
                  float visc, float dt) {
    addSource(N, u, u0, dt);
    addSource(N, v, v0, dt);
    std::swap(u0, u);
    diffuse(N, 1, u, u0, visc, dt);
    std::swap(v0, v);
    diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0);
    std::swap(u0, u);
    std::swap(v0, v);
    advect(N, 1, u, u0, u0, v0, dt);
    advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
}

}  // namespace

GridFluidSolver2::GridFluidSolver2() {}

GridFluidSolver2::~GridFluidSolver2() {}

void GridFluidSolver2::resizeGrid(const Size2& newSize,
                                     const Vector2D& newGridSpacing,
                                     const Vector2D& newGridOrigin) {
    UNUSED_VARIABLE(newGridSpacing);
    UNUSED_VARIABLE(newGridOrigin);

    _u.resize(newSize);
    _uTemp.resize(newSize);
    _v.resize(newSize);
    _vTemp.resize(newSize);
    _den.resize(newSize);
    _denTemp.resize(newSize);
}

ConstArrayAccessor2<float> GridFluidSolver2::density() const {
    return _den.constAccessor();
}

ArrayAccessor2<float> GridFluidSolver2::density() {
    return _den.accessor();
}

void GridFluidSolver2::onAdvanceTimeStep(double timeIntervalInSeconds) {
    size_t N = _u.size().x - 2;
    float visc = 0.0f;
    float diff = 0.0f;
    float dt = (float)timeIntervalInSeconds;

    _uTemp.set(0.0f);
    _vTemp.set(0.0f);
    _denTemp.set(0.0f);

    velocityStep(N, _u.data(), _v.data(), _uTemp.data(), _vTemp.data(), visc,
                 dt);
    densityStep(N, _den.data(), _denTemp.data(), _u.data(), _v.data(), diff,
                dt);
}
