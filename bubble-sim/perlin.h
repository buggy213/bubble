#pragma once
#include <array>

typedef std::array<int, 512> PerlinPermutation;

PerlinPermutation generate_permutation(int seed);
double perlin(double x, double y, double z, const PerlinPermutation &permutation);