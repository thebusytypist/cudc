#ifndef PINV_H
#define PINV_H

bool SVD2(
    const float* M,
    float* U,
    float* S,
    float* VT);

bool PseudoInverse2(
    const float* M,
    float* pinv);

#endif
