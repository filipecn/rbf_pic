#include "cell_centered_graph2_model.h"
#include "ponos_wrapper.h"

CellCenteredGraph2Model::CellCenteredGraph2Model(furoo::CellGraph2 *graph)
    : _graph(graph) {
  dt = 0.05;
  debugMode = false;
  nodeColor = aergia::Color(0, 0, 0, 0.4);
  regionColor = aergia::Color(1, 0, 0, 0.2);
}

void CellCenteredGraph2Model::draw() {
  auto cellMaterialField = *_graph->field<furoo::Definitions::Material>(2);
  auto cellPressureField = *_graph->field<double>(4);
  auto cellVelocityXField = *_graph->field<double>(8);
  auto cellVelocityYField = *_graph->field<double>(9);
  auto cellDivergenceField = *_graph->field<double>(5);
  auto cellGradientXField = *_graph->field<double>(6);
  auto cellGradientYField = *_graph->field<double>(7);
  auto faceVelocityField = *_graph->field<double>(10);
  auto facePressureField = *_graph->field<double>(11);
  auto surfaceField = *_graph->field<unsigned char>(3);
  auto faceMaterialField = *_graph->field<furoo::Definitions::Material>(1);
  aergia::ColorPalette pressurePalette = aergia::HEAT_MATLAB_PALETTE;
  aergia::TextRenderer text;
  _graph->iterateCells([&](size_t cellId) {
    auto region = _graph->cellRegion(cellId);
    auto cellColor = aergia::COLOR_TRANSPARENT;
    if (cellMaterialField[cellId] == furoo::Definitions::Material::FLUID)
      cellColor = aergia::Color(0, 0, 1, 0.3);
    if (surfaceField[cellId]) {
      cellColor = aergia::Color(0.75, 0, 1, 0.3);
    }
//    cellColor =
//        pressurePalette(1.f - ponos::linearStep(cellPressureField[cellId],
//                                                _minPressure,
//                                                _maxPressure));
//    cellColor.a = 0.6;
    furoo::Ponos::drawRegion(region, cellColor);
    if (debugMode) {
      std::ostringstream stringStream2;
      stringStream2 << cellId;
//          << "vx:" << cellVelocityXField[cellId]
//          << "vy: " << cellVelocityYField[cellId];
//          << "d: " << cellDivergenceField[cellId]
//          << "p: " << cellPressureField[cellId];
//          << "gx: " << cellGradientXField[cellId];
//          << "gy: " << cellGradientYField[cellId];

      auto p = _graph->cellRegion(cellId).center();
      auto labelPosition =
          aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
      text.render(stringStream2.str(),
                  labelPosition,
                  0.3f,
                  aergia::COLOR_BLACK);
    }
  });
  for (auto sc : selectedCells) {
    auto region = _graph->cellRegion(sc);
    furoo::Ponos::drawRegion(region, aergia::Color(1, 1, 0, 0.4));
  }
  {
//    _graph->iterateFaces([&](size_t id) {
//      auto p = _graph->faceCenterPosition(id);
//      glPointSize(5.);
//      glBegin(GL_POINTS);
//      glColor3f(1, 0, 0);
////      if (id == 9)
////        aergia::glVertex(furoo::Ponos::point2(p));
//      glColor4f(0, 0, 0, 0.2);
//      if (faceMaterialField[id] == furoo::Definitions::Material::FLUID) {
//        glColor3f(0, 0, 1);
//        aergia::glVertex(furoo::Ponos::point2(p));
//      }
//      glEnd();
//      glPointSize(1.);
//      std::ostringstream stringStream2;
//      stringStream2 << id;
//      glLineWidth(2.);
//
//      aergia::glColor(aergia::COLOR_RED);
////      if (_graph->faceOrientation(id)
////          == furoo::Definitions::Orientation::HORIZONTAL)
////        aergia::draw_vector(furoo::Ponos::point2(p),
////                            ponos::Vector2(0., faceVelocityField[id]) * dt);
////      else
////        aergia::draw_vector(furoo::Ponos::point2(p),
////                            ponos::Vector2(faceVelocityField[id], 0) * dt);
//      glLineWidth(1.);
////      if (debugMode) {
//      if (id == 9) {
//        auto labelPosition =
//            aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
//        _text->render(stringStream2.str(),
//                      labelPosition,
//                      0.6f,
//                      aergia::COLOR_RED);
////      }
//      }
//    });
  }
//  return;
  if (debugMode)
    _graph->iterateFaces([&](size_t id) {
      std::ostringstream stringStream2;
      stringStream2 << faceVelocityField[id];
      auto p = _graph->faceCenterPosition(id);
      auto labelPosition =
          aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
      text.render(stringStream2.str(),
                  labelPosition,
                  0.3f,
                  aergia::COLOR_BLACK);

    });
//  _graph->iterateCells([&](size_t id) {
//    std::ostringstream stringStream2;
//    stringStream2 << cellDivergenceField[id];
////      stringStream2 << cellVelocityXField[id] << ", g:"
////                    << cellVelocityYField[id];
//    auto p = _graph->cellCenterPosition(id);
////      aergia::draw_vector(furoo::Ponos::point2(p),
////                          ponos::Vector2(cellVelocityXField[id],
////                                         cellVelocityYField[id]) * dt);
//    auto labelPosition =
//        aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
//    _text->render(stringStream2.str(),
//                  labelPosition,
//                  0.3f,
//                  aergia::COLOR_WHITE);
//  });
}

void CellCenteredGraph2Model::updatePressure() {
  auto cellPressureField = *_graph->field<double>(4);
  _minPressure = INFINITY;
  _maxPressure = -INFINITY;
  _graph->iterateCells([&](size_t id) {
    _minPressure = std::min(_minPressure, cellPressureField[id]);
    _maxPressure = std::max(_maxPressure, cellPressureField[id]);
  });
  std::cerr << "pressure updated " << _minPressure << " " << _maxPressure
            << std::endl;
}