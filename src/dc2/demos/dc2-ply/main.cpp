#include <dc2.h>
#include <cstdio>
using namespace std;

int main(int argc, char* argv[]) {
    const int pcap = 1024;
    const int ecap = 1024;

    float p[pcap * 2];
    int edges[ecap * 2];

    int pcnt, ecnt;

    bool r = DualContour2(FT_UNIT_SPHERE,
        -2.0f, 2.0f,
        -2.0f, 2.0f,
        4,
        p, pcap, &pcnt,
        edges, ecap, &ecnt);

    printf("pcnt: %d\necnt: %d\n", pcnt, ecnt);

    return 0;
}