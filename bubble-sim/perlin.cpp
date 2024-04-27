#include "perlin.h"
#include <algorithm>
#include <cmath>
#include <random>

// based on http://flafla2.github.io/2014/08/09/perlinnoise.html

// 6t^5 -15t^4 + 10t^3
static double fade(double t) {
    return t * t * t * (t * (t * 6 - 15) + 10);
}

static double grad(int hash, double x, double y, double z)
{
    switch(hash & 0xF)
    {
        case 0x0: return  x + y;
        case 0x1: return -x + y;
        case 0x2: return  x - y;
        case 0x3: return -x - y;
        case 0x4: return  x + z;
        case 0x5: return -x + z;
        case 0x6: return  x - z;
        case 0x7: return -x - z;
        case 0x8: return  y + z;
        case 0x9: return -y + z;
        case 0xA: return  y - z;
        case 0xB: return -y - z;
        case 0xC: return  y + x;
        case 0xD: return -y + z;
        case 0xE: return  y - x;
        case 0xF: return -y - z;
        default: return 0; // never happens
    }
}

static double lerp(double a, double b, double x) {
    return a + x * (b - a);
}

// output ranges from [-1.0, 1.0]
double perlin(double x, double y, double z, const PerlinPermutation &p) {
    int xi = (int)x & 255;
    int yi = (int)y & 255;
    int zi = (int)z & 255;
    double xf = x - std::floor(x);
    double yf = y - std::floor(y);
    double zf = z - std::floor(z);

    double u = fade(xf);
    double v = fade(yf);
    double w = fade(zf);

    int aaa, aba, aab, abb, baa, bba, bab, bbb;
    aaa = p[p[p[xi  ]+yi  ]+zi  ];
    aba = p[p[p[xi  ]+yi+1]+zi  ];
    aab = p[p[p[xi  ]+yi  ]+zi+1];
    abb = p[p[p[xi  ]+yi+1]+zi+1];
    baa = p[p[p[xi+1]+yi  ]+zi  ];
    bba = p[p[p[xi+1]+yi+1]+zi  ];
    bab = p[p[p[xi+1]+yi  ]+zi+1];
    bbb = p[p[p[xi+1]+yi+1]+zi+1];

    double x1, x2, y1, y2;
    x1 = lerp(grad(aaa, xf, yf, zf), grad(baa, xf-1, yf, zf), u);
    x2 = lerp(grad(aba, xf, yf-1, zf), grad(bba, xf-1, yf-1, zf), u);
    y1 = lerp(x1, x2, v);

    x1 = lerp(grad(aab, xf, yf, zf-1), grad(bab, xf-1, yf, zf-1), u);
    x2 = lerp(grad(abb, xf, yf-1, zf-1), grad(bbb, xf-1, yf-1, zf-1), u);
    y2 = lerp(x1, x2, v);

    return lerp(y1, y2, w);
}

PerlinPermutation generate_permutation(int seed=-1) {
    std::mt19937 rng;
    if (seed != -1) {
        rng.seed(seed);
    }
    else {
        std::random_device rd;
        rng.seed(rd());
    }

    PerlinPermutation p;
    for (int i = 0; i < p.size(); i += 1) {
        p[i] = i;
    }
    
    std::shuffle(p.begin(), p.end(), rng);
    return p;
}