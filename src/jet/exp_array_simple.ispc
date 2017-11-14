// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

export void ispcAxpy(uniform float a, uniform float x[], uniform float y[],
                     uniform float result[], uniform int count) {
    foreach (index = 0 ... count) { result[index] = a * x[index] + y[index]; }
}
