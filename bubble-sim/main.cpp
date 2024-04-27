#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <igl/upsample.h>

#include "imgui.h"
#include "meshio.h"
#include "simulator.h"
// Inline mesh of a cube
  
const Eigen::MatrixXd cube_V = (Eigen::MatrixXd(8,3)<<
  0.0,0.0,0.0,
  0.0,0.0,1.0,
  0.0,1.0,0.0,
  0.0,1.0,1.0,
  1.0,0.0,0.0,
  1.0,0.0,1.0,
  1.0,1.0,0.0,
  1.0,1.0,1.0).finished();

const Eigen::MatrixXi cube_F = (Eigen::MatrixXi(12,3)<<
  0,6,4,
  0,2,6,
  0,3,2,
  0,1,3,
  2,7,6,
  2,3,7,
  4,6,7,
  4,7,5,
  0,4,5,
  0,5,1,
  1,5,7,
  1,7,3).finished();

int main(int argc, char *argv[])
{
  Eigen::MatrixXd subdivided_V;
  Eigen::MatrixXi subdivided_F;

  igl::upsample(cube_V, cube_F, subdivided_V, subdivided_F, 4);

  int fps = 30;
  int steps_per_frame = 5;
  double beta = 0.05;
  double delta_t = 1.0 / (fps * steps_per_frame);

  Simulator sim {
    std::move(subdivided_V),
    std::move(subdivided_F),
    beta,
    delta_t
  };
  
  igl::opengl::glfw::Viewer viewer;

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);

  // Customize the menu

  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();

    // Add new group
    if (ImGui::CollapsingHeader("Simulation Parameters", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // Expose variable directly ...
      bool any_changed = false;
      any_changed |= ImGui::InputDouble("beta", &beta, 0, 0, "%.4f");
      any_changed |= ImGui::InputInt("FPS", &fps, 0, 0);
      any_changed |= ImGui::InputInt("steps per frame", &steps_per_frame, 0, 0);
      sim.display_stats();

      if (any_changed) {
        delta_t = 1.0 / (fps * steps_per_frame);
        sim.set_params(SimParameters(beta, delta_t));
      }
    }
  };


  // Plot the mesh
  viewer.data().set_mesh(sim.get_verts(), sim.get_faces());
  viewer.data().set_face_based(true);

  bool running = false;
  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifier) -> bool {
    if (key == '.') {
      sim.step();
      viewer.data().clear();
      viewer.data().set_mesh(sim.get_verts(), sim.get_faces());
    }

    if (key == 'r') {
      running = !running;
      viewer.core().is_animating = running;
    }

    if (key == 's') {
      save_mesh_as_obj("test.obj", sim.get_verts(), sim.get_faces());
    }
    
    return false;
  };

  viewer.callback_post_draw = [&](igl::opengl::glfw::Viewer& viewer) -> bool {
    if (running) {
      sim.step();
      
      viewer.data().clear();
      viewer.data().set_mesh(sim.get_verts(), sim.get_faces());
    }

    return false;
  };

  viewer.launch();
}
