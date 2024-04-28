#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <igl/upsample.h>
#include <optional>
#include <getopt.h>

#include "igl/opengl/ViewerData.h"
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



struct ApplicationArgs {
  std::optional<std::string> input_mesh;
  std::optional<std::string> output_folder;
};

ApplicationArgs parse_args(int argc, char *argv[]) {
  ApplicationArgs args;

  int c;
  while (1) {
    static struct option long_options[] = {
      {"mesh", required_argument, 0, 'm'},
      {"output", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    int option_index = 0;
    c = getopt_long(argc, argv, "m:o:", long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch(c) {
      case 'm':
        args.input_mesh = std::make_optional(optarg);
        break;
      case 'o':
        args.output_folder = std::make_optional(optarg);
        break;
    }
  }

  return args;
}

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  ApplicationArgs args = parse_args(argc, argv);

  if (auto &mesh_filename = args.input_mesh) {
    if (!read_obj_as_mesh(*mesh_filename, V, F)) {
      std::cerr << "unable to parse obj file" << std::endl;
      return 1;
    }
  }
  else {
    // use default cube if none given
    igl::upsample(cube_V, cube_F, V, F, 4);
  }

  int fps = 30;
  int steps_per_frame = 5;
  double beta = 1.0;
  double delta_t = 1.0 / (fps * steps_per_frame);
  double gravity = 0.0;
  bool visualize_wind = false;

  SimParameters sim_params{
    beta,
    delta_t,
    gravity
  };
  Simulator sim {
    std::move(V),
    std::move(F),
    sim_params
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
    // (don't) Draw parent menu content
    // menu.draw_viewer_menu();

    // Add new group
    if (ImGui::CollapsingHeader("Simulation Parameters", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // Expose variable directly ...
      bool any_changed = false;
      any_changed |= ImGui::InputDouble("beta", &beta, 0, 0, "%.4f");
      any_changed |= ImGui::InputInt("FPS", &fps, 0, 0);
      any_changed |= ImGui::InputInt("steps per frame", &steps_per_frame, 0, 0);
      any_changed |= ImGui::InputDouble("gravity", &gravity, 0, 0, "%.4f");
      any_changed |= ImGui::Checkbox("visualize wind", &visualize_wind);
      sim.display_stats();

      if (any_changed) {
        if (visualize_wind) {
          viewer.data().clear_edges();
          sim.visualize_wind([&](const Eigen::MatrixXd &P1, const Eigen::MatrixXd &P2, const Eigen::MatrixXd &C){
            viewer.data().add_edges(P1, P2, C);
            // viewer.data().add_points(P1, C);
          });
        }
        else {
          viewer.data().clear_edges();
        }

        delta_t = 1.0 / (fps * steps_per_frame);
        sim_params = {beta, delta_t, gravity};
        sim.set_params(sim_params);
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
