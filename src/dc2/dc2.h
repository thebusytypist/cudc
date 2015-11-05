#ifndef DC2_H
#define DC2_H

enum FunctionType {
    FT_UNIT_SPHERE = 0,
    cFunctionTypeCount
};

template <FunctionType>
bool DualContour2(
    float xs, float xt,
    float ys, float yt, int n,
    float* p, int pcap,
    int* edge, int ecap);

#endif
