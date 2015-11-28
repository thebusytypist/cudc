#include <cudc2.h>
#include <cassert>
#include <cstdio>
#include <cstring>
using namespace std;

int main(int argc, char* argv[]) {
    const int pcap = 1024;
    const int ecap = 1024;

    float start, end;
    int ticks;
    FunctionType ft;
    assert(argc == 5);

    if (strncmp(argv[1], "UNIT_SPHERE", 11) == 0)
        ft = FT_UNIT_SPHERE;
    else if (strncmp(argv[1], "HEART", 5) == 0)
        ft = FT_HEART;
    else
        assert(false);

    sscanf(argv[2], "%f", &start);
    sscanf(argv[3], "%f", &end);
    sscanf(argv[4], "%d", &ticks);

    float p[pcap * 2];
    int edges[ecap * 2];

    int pcnt, ecnt;

    bool r = DualContour2(ft,
        start, end,
        start, end,
        ticks,
        p, pcap, &pcnt,
        edges, ecap, &ecnt);
    assert(r);

    printf("%d %d\n", pcnt, ecnt);

    for (int i = 0; i < pcnt; ++i) {
        printf("%f %f\n", p[i * 2], p[2 * i + 1]);
    }

    for (int i = 0; i < ecnt; ++i) {
        printf("%d %d\n", edges[i * 2], edges[i * 2 + 1]);
    }

    return 0;
}