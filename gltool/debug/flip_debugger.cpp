#include <ui.h>
//#include "cell_centered_graph2_model.h"
//#include "particle_system_model.h"

class GUI : public nanogui::Screen {
 public:
  GUI(/*furoo::AdaptiveCenteredPic2 &pic,
      ParticleSystemModel2 &pmodel,
      CellCenteredGraph2Model &gmodel*/)
  /*: nanogui::Screen(), _pic(pic), _pmodel(pmodel), _gmodel(gmodel)*/ {
    parser.addIntArgument("minLevel", 5);
    parser.addIntArgument("maxLevel", 7);
    parser.addStringArgument("config", std::string());
    parser.addStringArgument("particlesOutput", std::string("."));
    parser.addStringArgument("cellGraphOutput", std::string("."));
    parser.addStringArgument("substepsOutput", std::string("."));
    parser.addDoubleArgument("particleSpacing", 0.00390625);
    parser.addIntArgument("frames", 60);
    nanogui::Window *window = new nanogui::Window(this, "Options");
    window->setPosition(Eigen::Vector2i(15, 15));
    window->setLayout(new nanogui::GroupLayout());
    auto b = new nanogui::Button(window, "Load Config File");
    b->setCallback([&] {
      parser.parse(nanogui::file_dialog(
          {{"png", "Portable Network Graphics"}, {"txt", "Text file"}},
          false).c_str());
    });
    auto *paths = new nanogui::Widget(window);
    paths->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
                                            nanogui::Alignment::Middle,
                                            0, 2));
    b = new nanogui::Button(paths, "Particles Path", ENTYPO_ICON_ARCHIVE);
    b->setCallback([&]() {
      auto f = nanogui::folder_dialog();
      if (f.size()) {
        particlePath->setTooltip(f[0]);
        particlePath->setValue(f[0].substr(f[0].size() - 8, 8));
      }
    });
    particlePath = new nanogui::TextBox(paths);
    paths = new nanogui::Widget(window);
    paths->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
                                            nanogui::Alignment::Middle,
                                            0, 2));
    b = new nanogui::Button(paths, "Structure Path", ENTYPO_ICON_ARCHIVE);
    b->setCallback([&]() {
      auto f = nanogui::folder_dialog();
      if (f.size()) {
        structurePath->setTooltip(f[0]);
        structurePath->setValue(f[0].substr(f[0].size() - 8, 8));
      }
    });
    structurePath = new nanogui::TextBox(paths);

    curFrame = curSubFrame = 0;
    window = new nanogui::Window(this, "Frame Control");
    window->setPosition(Eigen::Vector2i(35, 15));
    window->setLayout(new nanogui::GroupLayout());
    paths = new nanogui::Widget(window);
    paths->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
                                            nanogui::Alignment::Middle,
                                            0, 2));
    b = new nanogui::Button(paths, "", ENTYPO_ICON_CONTROLLER_JUMP_TO_START);
    b->setCallback([&]() {
      curFrame = std::max(0, curFrame - 1);
      curSubFrame = 0;
      loadFrame(curFrame, curSubFrame);
    });
    b = new nanogui::Button(paths, "", ENTYPO_ICON_CHEVRON_THIN_LEFT);
    b->setCallback([&]() {
      curSubFrame = std::max(0, curSubFrame - 1);
      loadFrame(curFrame, curSubFrame);
    });
    curFrameLabel = new nanogui::TextBox(paths);
    curFrameLabel->setFixedWidth(80);
    curFrameLabel->setValue("");
    b = new nanogui::Button(paths, "", ENTYPO_ICON_CHEVRON_THIN_RIGHT);
    b->setCallback([&]() {
      curSubFrame = std::max(0, curSubFrame + 1);
      loadFrame(curFrame, curSubFrame);
    });
    b = new nanogui::Button(paths, "", ENTYPO_ICON_CONTROLLER_NEXT);
    b->setCallback([&]() {
      curFrame = std::max(0, curFrame + 1);
      curSubFrame = 0;
      loadFrame(curFrame, curSubFrame);
    });
  }
 private:
  void loadFrame(int frame, int subFrame) {
    curFrame = frame;
    curSubFrame = subFrame;
    std::stringstream ss;
    ss << frame << "." << subFrame;
    curFrameLabel->setValue(ss.str());
  }
  // simulation & rendering
//  furoo::AdaptiveCenteredPic2 &_pic;
//  ParticleSystemModel2 &_pmodel;
//  CellCenteredGraph2Model &_gmodel;

  int curSubFrame;
  furoo::ArgsParser parser;
  // strings
  nanogui::TextBox *particlePath, *structurePath, *curFrameLabel;
  int curFrame;
};

int main() {
  aergia::SceneApp<> app(800, 800, "FLIP Debugger");
  app.init();
//  furoo::AdaptiveCenteredPic2 pic(furoo::BBox2d::squareBox(1), 4, 6);
//  ParticleSystemModel2 particles(pic.particles(), 0.001);
//  CellCenteredGraph2Model
//      model(dynamic_cast<furoo::CellGraph2 *>(pic.structure()));
//  GUI gui(pic, particles, model);
  try {
    GUI gui;
    gui.initialize(aergia::GraphicsDisplay::instance().getGLFWwindow(), true);
    gui.setVisible(true);
    gui.performLayout();
    app.mouseCallback =
        [&](double x, double y) { gui.cursorPosCallbackEvent(x, y); };
    app.buttonCallback = [&](int button, int action, int mods) {
      gui.mouseButtonCallbackEvent(button, action, mods);
    };
    app.keyCallback = [&](int key, int scancode, int action, int mods) {
      gui.keyCallbackEvent(key,
                           scancode,
                           action,
                           mods);
    };
    app.charCallback = [&](unsigned int pointcode) {
      gui.charCallbackEvent(pointcode);
    };
    app.dropCallback = [&](int count, const char **filenames) {
      gui.dropCallbackEvent(count,
                            filenames);
    };
    app.renderCallback = [&]() {
      gui.drawContents();
      gui.drawWidgets();
    };
    app.run();
  } catch (std::string e) {
    std::cerr << e;
  }
  return 0;
}