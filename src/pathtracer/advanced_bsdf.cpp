#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 1
  // Implement MirrorBSDF
  *pdf = 1;
  reflect(wo, wi);
  return reflectance/abs_cos_theta(*wi);
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
      DragDouble3("Reflectance", &reflectance[0], 0.005);
      ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Project 3-2: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  double tan_2 = (sin_theta(h) / cos_theta(h)) * (sin_theta(h) / cos_theta(h));
  double e_val = exp(-1 * tan_2/(alpha * alpha));
  double denom = PI * alpha * alpha * pow(cos_theta(h), 4);
  
  return e_val / denom;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.
  
  Vector3D R_S = eta * eta + k * k - 2 * eta * cos_theta(wi) + cos_theta(wi) * cos_theta(wi);
  R_S = R_S / (eta * eta + k * k + 2 * eta * cos_theta(wi) + cos_theta(wi) * cos_theta(wi));
  
  Vector3D R_P = (eta * eta + k * k) * (cos_theta(wi) * cos_theta(wi)) - 2 * eta * cos_theta(wi) + 1;
  R_P = R_P / ((eta * eta + k * k) * (cos_theta(wi) * cos_theta(wi)) + 2 * eta * cos_theta(wi) + 1);
  return (R_S + R_P) / 2;
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Implement microfacet model here.
  
  Vector3D h = (wo + wi) * 0.5;
  h.normalize();
  Vector3D F_term = F(wi); //fresnel
  double G_term = G(wo, wi); //shadow masking
  double D_term = D(h);
  return (F_term * G_term * D_term)/ (4.0 * wo.z * wi.z);
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
    
  Vector2D uniform_two = sampler.get_sample();

  double theta_h = atan(sqrt(-1 * alpha * alpha * log(1 - uniform_two.x)));
  double phi_h = 2 * PI * uniform_two.y;
  //from previous section
  double half_x = cos(phi_h) * sin(theta_h);
  double half_y = sin(phi_h) * sin(theta_h);
  double half_z = cos(theta_h); // positive because face up on norm
  Vector3D h = Vector3D(half_x, half_y, half_z);
  
  // reflect out angle around the half
  *wi = 2 * h - wo;
  
  // calculate the pdfs
  double pdf_theta = (2.0 * sin(theta_h) * exp(-1 * tan(theta_h) * tan(theta_h) / (alpha * alpha))) / (alpha * alpha * pow(cos(theta_h), 3));
  float pdf_phi = 0.5 * PI;
  float pdf_fin_h = (pdf_theta * pdf_phi) / sin(theta_h);

  float pdf_fin_w = pdf_fin_h /( 4.0 * dot(*wi, h));
  *pdf = pdf_fin_w;
  Vector3D norm(0, 0, 1);
  if (dot(*wi, norm) <= 0 || dot(wo, norm) <= 0) {
      return Vector3D();
  }
  return MicrofacetBSDF::f(wo, *wi);
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 1
  // Implement RefractionBSDF
  if (!refract(wo, wi, ior)) {
      return Vector3D();
  }
  double eta = wo.z > 0 ? 1.0 / ior : ior;
  return transmittance / abs_cos_theta(*wi) / (eta * eta);
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO:
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305
  if (!refract(wo, wi, ior)) { // exists total internal reflection
    reflect(wo, wi);
    *pdf = 1.0;
    return reflectance / abs_cos_theta(*wi);
  }
  else {
    // use slick approximation
    double R_O = ((1.0 - ior) / (1.0 + ior)) * ((1.0 - ior) / (1.0 + ior));
    double R_val = R_O + (1.0 - R_O) * pow((1.0 - abs_cos_theta(wo)), 5);
    // coinflip
    if (coin_flip(R_val)) {
      reflect(wo, wi);
      *pdf = R_val;
      return R_val * reflectance / abs_cos_theta(*wi);
    }
    // false then refract
    else {
      refract(wo, wi, ior);
      *pdf = 1.0 - R_val;
      double eta = wo.z > 0 ? 1.0/ior : ior;
      return (1.0 - R_val) * transmittance / abs_cos_theta(*wi) / (eta * eta);
    }
  }
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {
  // TODO:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {
  // TODO Project 3-2: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double eta;
  double sign;
  if (wo.z < 0) { // incoming ray exits the refractive surface
      eta = ior;
      sign = 1;
  } else {
      eta = 1 / ior;
      sign = -1;
  }
  if ((1 - eta * eta * (1 - wo.z * wo.z) < 0)) {
      return false; //total internal refraction
  }
  // else can set wi
  wi->x = -1 * eta * wo.x;
  wi->y = -1* eta * wo.y;
  wi->z = sign * sqrt((1 - eta * eta * (1 - wo.z * wo.z)));
  return true;
}

} // namespace CGL
