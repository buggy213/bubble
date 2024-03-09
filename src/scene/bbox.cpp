#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  double x_slab_1 = (min.x - r.o.x) / r.d.x;
  double x_slab_2 = (max.x - r.o.x) / r.d.x;
  double x_slab_min = std::min(x_slab_1, x_slab_2);
  double x_slab_max = std::max(x_slab_1, x_slab_2);

  double y_slab_1 = (min.y - r.o.y) / r.d.y;
  double y_slab_2 = (max.y - r.o.y) / r.d.y;
  double y_slab_min = std::min(y_slab_1, y_slab_2);
  double y_slab_max = std::max(y_slab_1, y_slab_2);

  double z_slab_1 = (min.z - r.o.z) / r.d.z;
  double z_slab_2 = (max.z - r.o.z) / r.d.z;
  double z_slab_min = std::min(z_slab_1, z_slab_2);
  double z_slab_max = std::max(z_slab_1, z_slab_2);

  double t_min = std::max(x_slab_min, std::max(y_slab_min, z_slab_min));
  double t_max = std::min(x_slab_max, std::min(y_slab_max, z_slab_max));

  if (t_min >= t_max) { return false; }
  if (t_min > t1 || t_max < t0) { return false; } 

  t0 = t_min;
  t1 = t_max;
  return true;
  
}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
