#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"
#include "misc.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.

  // moller-trumbore algorithm mostly just transcribed from wikipedia
  Vector3D e1 = p2 - p1;
  Vector3D e2 = p3 - p1;
  Vector3D ray_cross_e2 = cross(r.d, e2);
  double det = dot(e1, ray_cross_e2);
  if (abs(det) < EPS_D) { return false; }

  double inv_det = 1.0 / det;
  Vector3D s = r.o - p1;
  double u = inv_det * dot(s, ray_cross_e2);
  if (u < 0 || u > 1) { return false; }

  Vector3D s_cross_e1 = cross(s, e1);
  double v = inv_det * dot(r.d, s_cross_e1);
  if (v < 0 || u + v > 1) { return false; }

  double ray_t = inv_det * dot(e2, s_cross_e1);
  if (ray_t < r.min_t || ray_t > r.max_t) {
    return false;
  }

  // ray intersected triangle, and it was earliest such intersection so far
  r.max_t = ray_t;

  return true;
}

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly

  // moller-trumbore algorithm mostly just transcribed from wikipedia
  Vector3D e1 = p2 - p1;
  Vector3D e2 = p3 - p1;
  Vector3D ray_cross_e2 = cross(r.d, e2);
  double det = dot(e1, ray_cross_e2);
  if (abs(det) < EPS_D) { return false; }

  double inv_det = 1.0 / det;
  Vector3D s = r.o - p1;
  double u = inv_det * dot(s, ray_cross_e2);
  if (u < 0 || u > 1) { return false; }

  Vector3D s_cross_e1 = cross(s, e1);
  double v = inv_det * dot(r.d, s_cross_e1);
  if (v < 0 || u + v > 1) { return false; }

  double ray_t = inv_det * dot(e2, s_cross_e1);
  if (ray_t < r.min_t || ray_t > r.max_t) {
    return false;
  }

  // ray intersected triangle, and it was earliest such intersection so far
  r.max_t = ray_t;
  
  double tri_t = 1.0 - u - v;
  Vector3D interp_n = n1 * tri_t + n2 * u + n3 * v;
  interp_n = interp_n.unit();
  
  isect->t = ray_t;
  isect->n = interp_n;
  isect->primitive = this;
  isect->bsdf = get_bsdf();

  return true;
}

void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
