#pragma once
#include <array>
#include <Eigen/Eigen>

#include "vmath.hpp"

struct NoiseSample {
    float value;
    vmath::Vector3 derivative;
};

typedef std::array<int, 512> PerlinPermutation;

PerlinPermutation generate_permutation(int seed);
NoiseSample perlin(const vmath::Vector3 &point, float frequency, const PerlinPermutation &permutation);
NoiseSample perlin(float x, float y, float z, float frequency, const PerlinPermutation &permutation);
void perlin_curl(
    const Eigen::MatrixXd &vertices, 
    float frequency,
    Eigen::MatrixXd &curl_noise,
    const PerlinPermutation &x,
    const PerlinPermutation &y,
    const PerlinPermutation &z
);