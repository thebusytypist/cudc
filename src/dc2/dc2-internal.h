#ifndef DC2_INTERNAL_H
#define DC2_INTERNAL_H

#include "dc2.h"

template <FunctionType>
void Sample2(
    const Function& f,
    const float* x, const float* y, float* s, int n);

template <FunctionType>
void SampleGradient2(
    const Function& f,
    const float* x, const float* y,
    float* ds0, float* ds1,
    int n);

template <FunctionType>
void CollectIntersectionEdges(
    const Function& f,
    const float* x0, const float* y0,
    const float* x1, const float* y1,
    const float* v0, const float* v1, int n,
    float* xlow, float* ylow,
    float* xhigh, float* yhigh,
    int* ens,
    int* en);

template <FunctionType>
void SolveIntersection(
    const Function& f,
    const float* xlow, const float* ylow,
    const float* xhigh, const float* yhigh,
    float* x,
    float* y,
    int n);

#endif
