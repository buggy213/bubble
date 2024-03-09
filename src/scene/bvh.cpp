#include "bvh.h"

#include "CGL/CGL.h"
#include "misc.h"
#include "scene/primitive.h"
#include "triangle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

enum SplitMethod {
  Midpoint,
  Median,
  SAH
};

std::vector<Primitive*>::iterator split_primitives(
  std::vector<Primitive *>::iterator start,
  std::vector<Primitive *>::iterator end,
  int dim,
  BBox& centroid_bounds,
  SplitMethod method
) {
  switch (method) {
  case Midpoint:
  {
    double mid;
    mid = centroid_bounds.centroid()[dim];
    auto split_point = std::partition(start, end, [dim, mid](Primitive* p){
      double v = p->get_bbox().centroid()[dim];
      return v < mid;
    });

    if (split_point != start && split_point != end) {
      return split_point;
    }
    else {
      return split_primitives(start, end, dim, centroid_bounds, SplitMethod::Median);
    }
  }
  case Median:
  {
    int n = std::distance(start, end);
    auto mid = start + (n / 2);
    std::nth_element(start, mid, end, [dim](Primitive* a, Primitive* b) {
      double a_v = a->get_bbox().centroid()[dim];
      double b_v = b->get_bbox().centroid()[dim];
      return a_v < b_v;
    });
    
    return mid;
  }
  case SAH:
    printf("SAH not implemented yet");
    fflush(stdout);
    exit(-1);
    return start;
  }

  return start;
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox bbox;

  for (auto p = start; p != end; p++) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
  }

  BVHNode *node = new BVHNode(bbox);
  if (std::distance(start, end) <= max_leaf_size) {
    // leaf node
    node->start = start;
    node->end = end;

    return node;
  }
  else {
    // basically lifted from PBR book
    // 1. choose axis based on extent of projection of centroids
    BBox centroid_bbox;
    for (auto p = start; p != end; p++) {
      Vector3D centroid = (*p)->get_bbox().centroid();
      centroid_bbox.expand(centroid);
    }

    double extent_x = centroid_bbox.extent.x;
    double extent_y = centroid_bbox.extent.y;
    double extent_z = centroid_bbox.extent.z;

    if (centroid_bbox.empty()) {
      // give up
      node->start = start;
      node->end = end;

      return node;
    }
    int split_dim;
    if (extent_x > extent_y && extent_x > extent_z) {
      // partition along x
      split_dim = 0;
    }
    else if (extent_y > extent_x && extent_y > extent_z) {
      // partition along y
      split_dim = 1;
    }
    else {
      // partition along z
      split_dim = 2;
    }

    auto split_point = split_primitives(
      start, 
      end, 
      split_dim, 
      bbox, 
      SplitMethod::Midpoint
    );

    auto left = construct_bvh(start, split_point, max_leaf_size);
    auto right = construct_bvh(split_point, end, max_leaf_size);

    node->l = left;
    node->r = right;

    return node;
  }
}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.



  for (auto p : primitives) {
    total_isects++;
    if (p->has_intersection(ray))
      return true;
  }
  return false;


}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.



  bool hit = false;
  for (auto p : primitives) {
    total_isects++;
    hit = p->intersect(ray, i) || hit;
  }
  return hit;


}

} // namespace SceneObjects
} // namespace CGL
