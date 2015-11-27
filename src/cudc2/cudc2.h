#ifndef CUDC2_H
#define CUDC2_H

enum FunctionType {
    FT_UNIT_SPHERE = 0,
    cFunctionTypeCount
};

/**
  \brief Dual contour a 2D function.

  @param xs the starting x coordinates.
  @param xt the ending x coordinates.
  @param ys the starting y coordinates.
  @param yt the ending y coordinates.

  @param[out] p the generated mesh vertices stored as array of float2.
  @param pcap the capacity of the vertices.
  @param[out] pcnt the actual count of generated vertices.

  @param[out] edges the generated edges stored as array int2.
  @param ecap the capacity of the edges.
  @param[out] ecnt the actual count of generated edges.
*/
bool DualContour2(
    FunctionType ft,
    float xs, float xt,
    float ys, float yt, int n,
    float* p, int pcap, int* pcnt,
    int* edges, int ecap, int* ecnt);
    
#endif