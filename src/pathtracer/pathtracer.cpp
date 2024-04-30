#include "pathtracer.h"

#include "misc.h"
#include "pathtracer/intersection.h"
#include "scene/bvh.h"
#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"
#include "util/random_util.h"
#include "vector3D.h"

using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  albedoBuffer.resize(width, height);
  normalBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);

  albedoBuffer.clear();
  normalBuffer.clear();
  albedoBuffer.resize(0, 0);
  normalBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

void PathTracer::write_to_hdr_buffer(HDRImageBuffer &hdr_buffer, size_t x0, size_t y0,
                                      size_t x1, size_t y1) {
  sampleBuffer.blit(hdr_buffer, x0, y0, x1, y1);
}

void PathTracer::write_to_albedo_buffer(HDRImageBuffer &albedo_buffer, size_t x0, size_t y0,
                                          size_t x1, size_t y1) {
  albedoBuffer.blit(albedo_buffer, x0, y0, x1, y1);
}


void PathTracer::write_to_normal_buffer(HDRImageBuffer &normal_buffer, size_t x0, size_t y0,
                                          size_t x1, size_t y1) {
  normalBuffer.blit(normal_buffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling,
  // you may find the "glow" around the light source is gone. This is totally
  // fine: the area lights in importance sampling has directionality, however in
  // hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out;

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead
  // of normal shading
  for (int sample = 0; sample < num_samples; sample += 1) {
    Vector3D w_in = hemisphereSampler->get_sample();
    double pdf = 1.0 / (2.0 * PI); // uniform over 2pi steradians of hemisphere

    Ray shadow_ray{
        hit_p,
        o2w * w_in}; // transform back to world coordinates before casting
    shadow_ray.min_t = EPS_F; // avoid numerical precision issues
    Intersection shadow_ray_isect;
    bool hit = bvh->intersect(shadow_ray, &shadow_ray_isect);
    if (hit) {
      Vector3D Li = zero_bounce_radiance(
          shadow_ray, shadow_ray_isect); // radiance from light, if any
      // abs here is unnecessary since sample will always be from upper
      // hemisphere, but this is still ok
      double cos_theta = std::abs(
          w_in.z); // since z aligned with normal in local coordinate system,
                   // dot product w/ (0, 0, 1) is just z-coordinate
      L_out +=
          (isect.bsdf->f(w_out, w_in) * Li * cos_theta / pdf) / num_samples;
    }
  }

  return L_out;
}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;

  for (auto light : scene->lights) {
    int num_samples = light->is_delta_light() ? 1 : ns_area_light;
    for (int s = 0; s < num_samples; s += 1) {
      Vector3D w_in;
      double light_distance;
      double light_pdf;
      Vector3D Li = light->sample_L(hit_p, &w_in, &light_distance, &light_pdf);

      if (dot(w_in, isect.n) < 0.0) {
        continue; // light is behind surface, cannot contribute
      }

      Ray shadow_ray{hit_p, w_in}; // w_in already in world coordinates
      shadow_ray.min_t = EPS_F;    // avoid numerical precision issues
      shadow_ray.max_t =
          light_distance - EPS_F; // don't want to intersect light
      bool hit = bvh->has_intersection(shadow_ray);
      if (hit) {
        continue; // light is occluded
      }

      // point lights have pdf of 1 to ensure that this works
      Vector3D w_in_local = w2o * w_in;
      // abs here is unnecessary due to earlier check, but this should still be
      // ok
      double cos_theta = std::abs(w_in_local.z);
      L_out += (isect.bsdf->f(w_out, w_in_local) * Li * cos_theta / light_pdf) /
               num_samples;
    }
  }

  return L_out;
}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light

  // not exactly sure what units are supposed to be here,
  // TODO: think about this maybe
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D emission = isect.bsdf->get_emission(w_out);
  return emission;
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

  if (direct_hemisphere_sample) {
    return estimate_direct_lighting_hemisphere(r, isect);
  } else {
    return estimate_direct_lighting_importance(r, isect);
  }
}

struct Bounce {
  Vector3D dir; // local coordinates
  double pdf;
  Vector3D f;

  Ray ray; // world coordinates
  Intersection isect;
  bool hit;
};

static Bounce calculate_bounce(const Ray &r, const Intersection &isect, const BVHAccel *bvh) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D bounce_dir;
  double bounce_pdf;
  Vector3D f = isect.bsdf->sample_f(w_out, &bounce_dir, &bounce_pdf);
  Ray bounce_ray{hit_p, o2w * bounce_dir,
                 (int)(r.depth - 1)}; // convert to world coords before casting
  bounce_ray.min_t = EPS_F;
  Intersection bounce_isect;
  bool bounce_hit = bvh->intersect(bounce_ray, &bounce_isect);

  return Bounce {
    .dir = bounce_dir,
    .pdf = bounce_pdf,
    .f = f,
    .ray = bounce_ray,
    .isect = bounce_isect,
    .hit = bounce_hit
  };
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Vector3D L_out(0, 0, 0);

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this
  // point. Should be called recursively to simulate extra bounces.
  Vector3D L_direct;
  Bounce bounce;
  if (isect.bsdf->is_delta()) {
    bounce = calculate_bounce(r, isect, bvh);
    if (bounce.hit) {
      L_direct = zero_bounce_radiance(bounce.ray, bounce.isect);
    }
    else {
      if (envLight != nullptr) {
        L_direct = envLight->sample_dir(bounce.ray);
      }
    }
  } else {
    L_direct = one_bounce_radiance(r, isect);
  }

  if (russian_roulette) {
    bool roulette_terminate = !coin_flip(continuation_probability);
    if (roulette_terminate) {
      return L_direct;
    }
  }

  if (r.depth <= 1) { // out of bounces
    return L_direct;
  }

  Vector3D L_indirect{0.0, 0.0, 0.0};
  // don't recalculate bounce if we were forced to compute it for direct lighting (delta BSDFs)
  if (!isect.bsdf->is_delta()) {
    bounce = calculate_bounce(r, isect, bvh);
  }
  if (bounce.hit) {
    Vector3D Li_indirect =
        at_least_one_bounce_radiance(bounce.ray, bounce.isect);
    // bounce_dir.z might be negative (e.g. transmission)
    double cos_theta = std::abs(bounce.dir.z);
    L_indirect = bounce.f * Li_indirect * cos_theta / bounce.pdf;
    if (russian_roulette) {
      L_indirect /= continuation_probability;
    }
  }

  if (isAccumBounces) {
    if (!indirect_only || r.depth != max_ray_depth) {
      L_out = L_direct + L_indirect;
    } else {
      L_out = L_indirect;
    }
  } else {
    L_out = L_indirect;
  }

  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Vector3D _albedo;
  Vector3D _normal;
  return est_radiance_global_illumination(r, _albedo, _normal);
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r, Vector3D &albedo, Vector3D &normal) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  // L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

  if (!bvh->intersect(r, &isect)) {
    Vector3D environment_L = envLight ? envLight->sample_dir(r) : L_out;

    // not entirely clear what to return in this case
    albedo = environment_L; 
    normal = Vector3D();
    return environment_L;
  }

  // otherwise, isect is valid
  // OIDN documentation says albedo of 1 is ok for dielectric
  if (isect.bsdf->is_delta()) {
    albedo = Vector3D(1.0, 1.0, 1.0);
  }
  else {
    // otherwise, ask BSDF for advice
    albedo = isect.bsdf->albedo();
  }
  normal = isect.n;

  

  // TODO (Part 3): Return the direct illumination.

  // only situation to not include direct emission are
  // 1. not accumulating bounces and depth > 0
  // 2. indirect only
  if (!(r.depth > 0 && !isAccumBounces) && !indirect_only) {
    L_out = zero_bounce_radiance(r, isect);
  }
  if (r.depth > 0) {
    L_out += at_least_one_bounce_radiance(r, isect);
  }
  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct

  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"
  int num_samples = ns_aa;          // total samples to evaluate
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel

  Vector3D radiance;
  Vector3D albedo;
  Vector3D normal;

  int sample_count = 0;

  double s1 = 0.0, s2 = 0.0;

  for (sample_count = 0; sample_count < num_samples;) {
    for (int i = 0; i < samplesPerBatch; i += 1) {
      Vector2D random_offset = gridSampler->get_sample();
      double normalized_x =
          ((double)x + random_offset.x) / (double)sampleBuffer.w;
      double normalized_y =
          ((double)y + random_offset.y) / (double)sampleBuffer.h;
      Ray ray = camera->generate_ray(normalized_x, normalized_y);
      ray.depth = max_ray_depth;
      Vector3D ray_normal, ray_albedo;
      Vector3D ray_radiance = est_radiance_global_illumination(ray, ray_albedo, ray_normal);
      radiance += ray_radiance;
      albedo += ray_albedo;
      normal += ray_normal;
      s1 += ray_radiance.illum();
      s2 += ray_radiance.illum() * ray_radiance.illum();
    }

    sample_count += samplesPerBatch;
    if (adaptive_sampling) {
      double mean = s1 / sample_count;
      double variance = (s2 - ((s1 * s1) / sample_count)) / (sample_count - 1);

      double I = 1.96 * sqrt(variance) / sqrt(sample_count);
      if (I <= maxTolerance * mean) {
        break;
      }
    }
  }

  radiance /= sample_count;
  normal /= sample_count;
  albedo /= sample_count;
  sampleBuffer.update_pixel(radiance, x, y);
  normalBuffer.update_pixel(normal, x, y);
  albedoBuffer.update_pixel(albedo, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = sample_count;
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
