#ifndef DC2_H
#define DC2_H

enum FunctionType {
    FT_UNIT_SPHERE = 0,
    cFunctionTypeCount
};

bool DualContour2(
    FunctionType ft,
    float xs, float xt,
    float ys, float yt, int n,
    float* p, int pcap, int* pcnt,
    int* edge, int ecap, int* ecnt);

#endif
