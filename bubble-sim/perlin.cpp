#include "perlin.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

#include "vmath.hpp"

// based on https://github.com/rajabala/CurlNoise/blob/master/CurlNoise/Noise.cpp

// 6t^5 -15t^4 + 10t^3
using vmath::Vector3;

static const unsigned int s_kHashIndexMask = 255;
static const unsigned int s_kGradientIndexMask = 15;

static const Vector3 s_Gradients3D[]  =
{
	// center of cube to center of each of the 12 edges of a cube
	Vector3(1.0f, 1.0f, 0.0f),
	Vector3(-1.0f, 1.0f, 0.0f),
	Vector3(1.0f,-1.0f, 0.0f),
	Vector3(-1.0f,-1.0f, 0.0f),
	Vector3(1.0f, 0.0f, 1.0f),
	Vector3(-1.0f, 0.0f, 1.0f),
	Vector3(1.0f, 0.0f,-1.0f),
	Vector3(-1.0f, 0.0f,-1.0f),
	Vector3(0.0f, 1.0f, 1.0f),
	Vector3(0.0f,-1.0f, 1.0f),
	Vector3(0.0f, 1.0f,-1.0f),
	Vector3(0.0f,-1.0f,-1.0f),

	// Repeat some to skew distribution
	Vector3(1.0f, 1.0f, 0.0f),
	Vector3(-1.0f, 1.0f, 0.0f),
	Vector3(0.0f,-1.0f, 1.0f),
	Vector3(0.0f,-1.0f,-1.0f)
};

inline static float Dot(Vector3 g, float x, float y, float z) {
    return g.getX() * x + g.getY() * y + g.getZ() * z;
}

inline static float Smooth(float t) {
    // 6t^5 - 15t^4 + 10t^3
    return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

inline static float SmoothDerivative(float t) {
    // 30t^4 - 60t^3 + 30t^2
    return 30.0f * t * t * (t * (t - 2.0f) + 1.0f);
}

//-----------------------------------------------------------------------------------------------
// The PerlinNoise3 implementation calculates a noise value in the range [-1f, 1f] given a 3D position & frequency.
// It also calculates the analytical derivative (change wrt x, y & z) which can be used to quickly calculate the curl,
// when the potential field doesn't need to be modified.
NoiseSample perlin(const Vector3 &point, float frequency, const PerlinPermutation &permutation) {
  Vector3 p = point * frequency;
  int ix0 = static_cast<int>(floorf(p.getX()));
  int iy0 = static_cast<int>(floorf(p.getY()));
  int iz0 = static_cast<int>(floorf(p.getZ()));
  float tx0 = p.getX() - ix0;
  float ty0 = p.getY() - iy0;
  float tz0 = p.getZ() - iz0;
  float tx1 = tx0 - 1.0f;
  float ty1 = ty0 - 1.0f;
  float tz1 = tz0 - 1.0f;
  ix0 &= s_kHashIndexMask;
  iy0 &= s_kHashIndexMask;
  iz0 &= s_kHashIndexMask;
  int ix1 = ix0 + 1;
  int iy1 = iy0 + 1;
  int iz1 = iz0 + 1;

  int h0 = permutation[ix0];
  int h1 = permutation[ix1];
  int h00 = permutation[h0 + iy0];
  int h10 = permutation[h1 + iy0];
  int h01 = permutation[h0 + iy1];
  int h11 = permutation[h1 + iy1];

  // gradients at each of the 8 lattice points
  Vector3 g000 = s_Gradients3D[permutation[h00 + iz0] & 15];
  Vector3 g100 = s_Gradients3D[permutation[h10 + iz0] & 15];
  Vector3 g010 = s_Gradients3D[permutation[h01 + iz0] & 15];
  Vector3 g110 = s_Gradients3D[permutation[h11 + iz0] & 15];
  Vector3 g001 = s_Gradients3D[permutation[h00 + iz1] & 15];
  Vector3 g101 = s_Gradients3D[permutation[h10 + iz1] & 15];
  Vector3 g011 = s_Gradients3D[permutation[h01 + iz1] & 15];
  Vector3 g111 = s_Gradients3D[permutation[h11 + iz1] & 15];

  float v000 = Dot(g000, tx0, ty0, tz0);
  float v100 = Dot(g100, tx1, ty0, tz0);
  float v010 = Dot(g010, tx0, ty1, tz0);
  float v110 = Dot(g110, tx1, ty1, tz0);
  float v001 = Dot(g001, tx0, ty0, tz1);
  float v101 = Dot(g101, tx1, ty0, tz1);
  float v011 = Dot(g011, tx0, ty1, tz1);
  float v111 = Dot(g111, tx1, ty1, tz1);

  float dtx = SmoothDerivative(tx0);
  float dty = SmoothDerivative(ty0);
  float dtz = SmoothDerivative(tz0);
  float tx = Smooth(tx0);
  float ty = Smooth(ty0);
  float tz = Smooth(tz0);

  float a = v000;
  float b = v100 - v000;
  float c = v010 - v000;
  float d = v001 - v000;
  float e = v110 - v010 - v100 + v000;
  float f = v101 - v001 - v100 + v000;
  float g = v011 - v001 - v010 + v000;
  float h = v111 - v011 - v101 + v001 - v110 + v010 + v100 - v000;

  Vector3 da = g000;
  Vector3 db = g100 - g000;
  Vector3 dc = g010 - g000;
  Vector3 dd = g001 - g000;
  Vector3 de = g110 - g010 - g100 + g000;
  Vector3 df = g101 - g001 - g100 + g000;
  Vector3 dg = g011 - g001 - g010 + g000;
  Vector3 dh = g111 - g011 - g101 + g001 - g110 + g010 + g100 - g000;

  float value = a + b * tx + (c + e * tx) * ty + (d + f * tx + (g + h * tx) * ty) * tz;

  // calculate derivative analytically
  Vector3 derivative = da + db * tx + (dc + de * tx) * ty +
                      (dd + df * tx + (dg + dh * tx) * ty) * tz;

  Vector3 delta = Vector3((b + e * ty + (f + h * ty) * tz) * dtx,
                          (c + e * tx + (g + h * tx) * tz) * dty,
                          (d + f * tx + (g + h * tx) * ty) * dtz);

  derivative += delta;
  derivative *= frequency;
  
  
  return NoiseSample {.value = value, .derivative = derivative};
}

NoiseSample perlin(float x, float y, float z, float frequency, const PerlinPermutation& permutation) {
    return perlin(vmath::Vector3(x, y, z), frequency, permutation);
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
    for (int i = 0; i < 256; i += 1) {
        p[i] = i;
    }
    
    std::shuffle(p.begin(), p.begin() + 256, rng);
    for (int i = 256; i < 512; i += 1) {
        p[i] = p[i-256];
    }
    return p;
}

void print_permutation(const PerlinPermutation &p) {
    for (int i = 0; i < p.size(); i += 1) {
        std::cout << p[i] << " ";
    }
    std::cout << "\n";
    std::cout << p[0] << " " << p[256] << " \n"; 
}

// incompressible perlin-based curl noise
void perlin_curl(
    const Eigen::MatrixXd &vertices, 
    float frequency,
    Eigen::MatrixXd &curl_noise, 
    const PerlinPermutation &x, 
    const PerlinPermutation &y, 
    const PerlinPermutation &z
) {
    static constexpr double DELTA = 1e-8;

    // print_permutation(x);
    // print_permutation(y);
    // print_permutation(z);

    curl_noise.resize(vertices.rows(), vertices.cols());

    for (int i = 0; i < vertices.rows(); i += 1) {
        NoiseSample psi_x = perlin(vertices(i, 0), vertices(i, 1), vertices(i, 2), frequency, x);
        NoiseSample psi_y = perlin(vertices(i, 0), vertices(i, 1), vertices(i, 2), frequency, y);
        NoiseSample psi_z = perlin(vertices(i, 0), vertices(i, 1), vertices(i, 2), frequency, z);
        
        double curl_x = psi_z.derivative.getY() - psi_y.derivative.getZ();
        double curl_y = psi_x.derivative.getZ() - psi_z.derivative.getX();
        double curl_z = psi_y.derivative.getX() - psi_x.derivative.getY();
        
        curl_noise(i, 0) = curl_x;
        curl_noise(i, 1) = curl_y;
        curl_noise(i, 2) = curl_z;
    }
}