#ifndef DC2_INTERNAL_H
#define DC2_INTERNAL_H

#include "dc2.h"

/**
  \brief Sample the 2D mesh grid.

  @param[in] x0 starting x coordinates.
  @param[in] y0 starting y coordinates.
  @param[in] x1 ending x coordinates.
  @param[in] y1 ending y coordinates.
  @param[out] s sampled values.
  @param n number of mesh grid vertices.
*/
template <FunctionType>
void Sample2(
    const Function& f,
    float x0, float x1,
    float y0, float y1,
    float* s, int n);

/**
  \brief Sample the 2D mesh grid for the gradients.

  @param[in] x x coordinates.
  @param[in] y y coordinates.
  @param[out] ds0 sampled gradient along the x direction.
  @param[out] ds1 sampled gradient along the y direction.
  @param n number of mesh grid vertices.
*/
template <FunctionType>
void SampleGradient2(
    const Function& f,
    const float* x, const float* y,
    float* ds0, float* ds1,
    int n);

/**
  \brief Collect the intersection edges for each grid cell.

  @param[in] x0 x coordinates of the first row of cells.
  @param[in] y0 y coordinates of the first row of cells.
  @param[in] x1 x coordinates of the second row of cells.
  @param[in] y1 y coordinates of the second row of cells.
  @param[in] v0 sampled values of the first row of cells.
  @param[in] v1 sampled values of the second row of cells.
  @param n number of mesh grid vertices.
  @param[out] xlow x coordinates of vertices who has negative sampled values.
  @param[out] ylow y coordinates of vertices who has negative sampled values.
  @param[out] xhigh x coordinates of vertices
              who has non-negative sampled values.
  @param[out] yhigh y coordinates of vertices
              who has non-negative sampled values.
  @param[out] ens a bit field indicating how many edges are intersected
              in a grid cell.
  @param[out] en total number of intersected edges.
*/
void CollectIntersectionEdges2(
    const float* x0, const float* y0,
    const float* x1, const float* y1,
    const float* v0, const float* v1, int n,
    float* xlow, float* ylow,
    float* xhigh, float* yhigh,
    int* ens,
    int* en);

/**
  \brief Solve for the intersection points.

  @param[in] xlow x coordinates of vertices who has negative sampled values.
  @param[in] ylow y coordinates of vertices who has negative sampled values.
  @param[in] xhigh x coordinates of vertices
              who has non-negative sampled values.
  @param[in] yhigh y coordinates of vertices
              who has non-negative sampled values.
  @param[out] x x coordinates of the intersection points.
  @param[out] y y coordinates of the intersection points.
  @param n total number of intersected edges.
*/
template <FunctionType>
void SolveIntersection2(
    const Function& f,
    const float* xlow, const float* ylow,
    const float* xhigh, const float* yhigh,
    float* x,
    float* y,
    int n);

/**
  \brief Construct QEF using two rows intersection points.

  @param[in] ix0 x coordinates of the intersection points of the first row.
  @param[in] iy0 y coordinates of the intersection points of the first row.
  @param[in] ix1 x coordinates of the intersection points of the second row.
  @param[in] iy1 y coordinates of the intersection points of the second row.
  @param[in] nx0 gradient along the x direction of the first row.
  @param[in] ny0 gradient along the y direction of the first row.
  @param[in] nx1 gradient along the x direction of the second row.
  @param[in] ny1 gradient along the y direction of the second row.
  @param[in] ens0 a bit field indicating how many edges are intersected
             in a grid cell of the first row.
  @param[in] ens1 a bit field indicating how many edges are intersected
             in a grid cell of the second row.
  @param n number of cells in a row.
  @param[out] f an array of QEF coefficients.
              Each QEF has 7 coefficients.
  @param[out] h an array of size n - 1 indicating
              if the cell contributes a mesh vertex.
  @param[out] m total number of QEFs.
*/
void ConstructQEF2(
    const float* ix0, const float* iy0,
    const float* ix1, const float* iy1,
    const float* nx0, const float* ny0,
    const float* nx1, const float* ny1,
    const int* ens0, const int* ens1, int n,
    float* f, bool* h, int* m);

/**
  \brief Solve the QEFs for the mesh vertices.

  @param[in] f an array of QEF coefficients.
             Each QEF has 7 coefficients.
  @param[out] p the solution array. Each element is a 2D vertex.
  @param n total number of QEFs.
*/
void SolveQEF2(const float* f, float* p, int n);

#endif
