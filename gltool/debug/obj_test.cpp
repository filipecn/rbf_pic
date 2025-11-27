#include <ponos_wrapper.h>

#include <mesh_model.h>
#include <ui.h>

using namespace furoo;

template <int D> class OBJDebugger : public aergia::SceneObject {
public:
  OBJDebugger() {
    _screen.reset(new aergia::NanoGUIScreen());
    auto *_gui = new nanogui::FormHelper(_screen.get());
    nanogui::ref<nanogui::Window> win =
        _gui->addWindow(Eigen::Vector2i(10, 10), "Settings");
    win->setLayout(new nanogui::GroupLayout());
    auto *b = new nanogui::Button(win, "load");
    b->setCallback([&]() {
      auto file =
          nanogui::file_dialog({{"obj", "wavefront object file"}}, false);
      if (!file.empty()) {
        furoo::IO::loadOBJ<D>(file, &_mesh);
        _solidModel.reset(new SolidObject<D>(&_mesh));
        _meshModel.reset(new MeshModel<D>(_mesh));
        _bvh.reset(new BVH(_mesh));
        setGraph();
        static auto cube = ponos::RawMeshes::cube();
        _pointSet = std::make_shared<GraphicElementSet>(cube.get());
        updatePoints();
      }
    });
    _screen->setVisible(true);
    _screen->performLayout();
  }

  void setGraph();

  void updatePoints() {
    if (!_bvh)
      return;
    size_t N = 10;
    double step = 1. / N;
    _points.clear();
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) {
        if (D == 2)
          _points.emplace_back(i * step, j * step);
        else
          for (int k = 0; k < N; k++)
            _points.emplace_back(i * step, j * step, k * step);
      }
    _pointSet->resize(_points.size());
    // update selected field min and max values
    {
      // unactivate dead cells, unfortunatelly we have to kill them all first
      // and make active ones live again
      auto a = _pointSet->activeBuffer();
      memset(a, 0, sizeof(unsigned int) * _pointSet->count());
    }
    // update positions
    size_t bufferElementIndex = 0;
    for (auto point : _points) {
      auto pos = point;
      auto m = _pointSet->transformBuffer(bufferElementIndex);
      float size = 0.005;
      float t[16];
      (ponos::translate(ponos::vec3(pos[0], pos[1], (D == 3) ? pos[2] : 0)) *
       ponos::scale(size, size, size))
          .matrix()
          .column_major(t);
      for (size_t k = 0; k < 16; k++)
        m[k] = t[k];

      aergia::Color color = aergia::COLOR_BLACK;
      std::cerr << "test " << point << std::endl;
      if (_bvh->isInside(point))
        color = aergia::COLOR_RED;
      auto c = _pointSet->colorBuffer(bufferElementIndex);
      c[0] = color.r;
      c[1] = color.g;
      c[2] = color.b;
      c[3] = color.a;
      auto a = _pointSet->activeBuffer(bufferElementIndex);
      a[0] = 1;
      bufferElementIndex++;
    }
  }

  void draw(const aergia::CameraInterface *camera,
            ponos::Transform t) override {
    if (_meshModel) {
      _meshModel->draw(camera, t);
      _pointSet->draw(camera, t);
    }
    if (_solidModel)
      _solidModel->draw(camera, t);
    _graphModel.draw(camera, t);
  }

  std::vector<Point<double, D>> _points;
  std::shared_ptr<GraphicElementSet> _pointSet;
  CellGraphModel<D> _graphModel;
  std::shared_ptr<StructureInterface<D>> _graph;
  std::shared_ptr<BVH<D>> _bvh;
  std::shared_ptr<aergia::NanoGUIScreen> _screen;
  std::shared_ptr<MeshModel<D>> _meshModel;
  std::shared_ptr<SolidObject<D>> _solidModel;
  furoo::Mesh<D> _mesh;
};

template <> inline void OBJDebugger<3>::setGraph() {
  _graph.reset(new CellGraph3(furoo::BBox3d::squareBox(1.f)));
  _graph->addCellField<int>(0);
  auto graph = dynamic_cast<CellGraph3 *>(_graph.get());
  std::queue<size_t> q;
  graph->refine(1, [&](const CellGraph3::Node node) {
    if (node.level() > 5)
      return false;
    std::cerr << "testing point " << node.region().center() << std::endl;
    if (_bvh->isInside(node.region().center()))
      return true;
    return false;
  });
  _graphModel.set(_graph.get());
}

int main() {
  aergia::SceneApp<> app(800, 800, "OBJ");
  OBJDebugger<3> obj;
  app.scene.add(&obj);
  app.scene.add(new aergia::CartesianGrid(5));
  app.run();
  return 0;
}