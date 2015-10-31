#ifndef PINV_H
#define PINV_H

bool svd2(
    const float* M,
    float* U,
    float* S,
    float* VT);

bool pinv2(
    const float* M,
    float* pinv);

#endif
