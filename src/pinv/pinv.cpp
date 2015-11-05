#include "pinv.h"
#include <cmath>
using namespace std;

bool SVD2(
    const float* M,
    float* U,
    float* S,
    float* VT) {
    const float e = 0.5f * (M[0] + M[3]);
    const float f = 0.5f * (M[0] - M[3]);
    const float g = 0.5f * (M[1] + M[2]);
    const float h = 0.5f * (M[1] - M[2]);

    // one half of w1 + w2
    const float hw1pw2 = sqrt(e * e + h * h);
    const float hw1mw2 = sqrt(f * f + g * g);

    // cos(y - b)
    const float cymb = hw1mw2 != 0.0f ? f / hw1mw2 : 1.0f;

    // cos(y + b)
    const float cypb = hw1pw2 != 0.0f ? e / hw1pw2 : 1.0f;

    // | cos(y) + cos(b) |
    const float cc = sqrt((1 + cymb) * (1 + cypb));
    // | cos(y) - cos(b) |
    const float ss = sqrt((1 - cymb) * (1 - cypb));
    // | sin(y) - sin(b) |
    const float cs = sqrt((1 + cymb) * (1 - cypb));
    // | sin(y) + sin(b) |
    const float sc = sqrt((1 - cymb) * (1 + cypb));

    const float cb = (cc - ss) * 0.5f;
    const float sb = (sc + cs) * 0.5f;

    U[0] = cb;
    U[1] = -sb;
    U[2] = sb;
    U[3] = cb;

    S[0] = hw1pw2 + hw1mw2;
    S[1] = hw1pw2 - hw1mw2;

    const float invS0 = 1.0f / S[0];
    const float invS1 = 1.0f / S[1];

    const float v0 = invS0 * U[0];
    const float v1 = invS0 * U[2];
    const float v2 = invS1 * U[1];
    const float v3 = invS1 * U[3];

    VT[0] = v0 * M[0] + v1 * M[2];
    VT[1] = v0 * M[1] + v1 * M[3];
    VT[2] = v2 * M[0] + v3 * M[2];
    VT[3] = v2 * M[1] + v3 * M[3];

    return true;
}

bool PseudoInverse2(
    const float* M,
    float* pinv) {
    const float EPS = 1e-12f;

    float U[4], S[2], VT[4];
    SVD2(M, U, S, VT);

    S[0] = abs(S[0]) < EPS ? 0.0f : 1.0f / S[0];
    S[1] = abs(S[1]) < EPS ? 0.0f : 1.0f / S[1];

    const float v0 = S[0] * VT[0];
    const float v1 = S[1] * VT[2];
    const float v2 = S[0] * VT[1];
    const float v3 = S[1] * VT[3];

    pinv[0] = v0 * U[0] + v1 * U[1];
    pinv[1] = v0 * U[2] + v1 * U[3];
    pinv[2] = v2 * U[0] + v3 * U[1];
    pinv[3] = v2 * U[2] + v3 * U[3];

    return true;
}