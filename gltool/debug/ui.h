#ifndef FUROO_UI_H
#define FUROO_UI_H

#include "utils.h"
#include <cell_graph2_model.h>
#include <fstream>
#include <furoo.h>
#include <iomanip>
#include <nanogui_utils.h>
#include <particle_system_model.h>
#include <solid_model.h>

template <int D> class UI : public aergia::SceneObject {
public:
  explicit UI();
  void init();
  void log(std::string line);
  // override
  void draw(const aergia::CameraInterface *camera, ponos::Transform t) override;
  bool intersect(const ponos::Ray3 &r, float *t) override;
  void dropPath(std::string path);
  void movePlane(int axis, int direction);
  // callbacks
  std::function<void(std::string)> saveFrameImageCallback;
  void loadConfigFile(std::string filename, std::string frameDir, int frameId);

private:
  void setupOptions();
  void saveFrameImage();
  void updateAll();
  // Simulation
  furoo::CustomSimulation<D> _sim;
  std::shared_ptr<SimulationModel<D>> _model;
  std::vector<std::shared_ptr<SolidModel<D>>> _solids;
  // GUI
  std::shared_ptr<aergia::NanoGUIScreen> _screen;
  std::shared_ptr<IOMenu<D>> _ioMenu;
  std::shared_ptr<SimulationPanel<D>> _simPanel;
  std::shared_ptr<SceneMenu<D>> _sceneMenu;
  std::shared_ptr<VisPanel<D>> _visPanel;
  std::shared_ptr<StepControl<D>> _stepControl;
  std::shared_ptr<SimulationStatistics<D>> _simStats;
  std::shared_ptr<TestFunctionPanel<D>> _testFunctionPanel;
  std::string _inputDirectory = "";
  // Statistics
  // Log system
  std::vector<nanogui::Label *> _logLines;
};

#include "ui.inl"

#endif // FUROO_UI_H
