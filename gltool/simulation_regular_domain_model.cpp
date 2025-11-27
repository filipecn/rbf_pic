#include "simulation_regular_domain_model.h"
#include <sstream>
#include "ponos_wrapper.h"

namespace furoo {

SimulationRegularDomain2Model::SimulationRegularDomain2Model(
    SimRegularDomain2 *d)
    : _domain(d) {
  _selectedCell = _selectedFace = -1;
  drawVelocities = false;
  drawPressure = false;
  drawDivergence = false;
  pressurePalette = aergia::HEAT_GREEN_COLOR_PALETTE;
  _text.reset(new aergia::Text(FONT_PATH));
  divergence.resize(_domain->cellCount(), 0.);
  _faces.reset(new PointZGrid2(16u, BBox2D::make_unit_bbox()));
  for (unsigned i = 0; i < _domain->faceCount(); i++)
    _faces->add(_domain->faceCenterPosition(i));
}

void SimulationRegularDomain2Model::draw() const {
  std::ostringstream stringStream;
  stringStream << "Selected cell: " << _selectedCell;
  ponos::Point3 labelPosition = aergia::glGetMVPTransform()(ponos::Point3(0.f, 0.f, 0.f));
  _text->render(stringStream.str(), labelPosition, 0.4f, aergia::COLOR_BLACK);
  if (_selectedCell >= 0) {
    std::ostringstream stringStream2;
    stringStream2 << "p: " << _domain->scalarAtCell(0, _selectedCell);
    auto p = _domain->cellCenterPosition(_selectedCell);
    labelPosition = aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
    _text->render(stringStream2.str(), labelPosition, 0.4f, aergia::COLOR_BLACK);
  }
  if(_selectedFace >= 0){
    std::ostringstream stringStream2;
    stringStream2 << "g: " << _domain->cellGradientAtFaceCenter(0, _selectedFace);
    auto p = _domain->faceCenterPosition(_selectedFace);
    labelPosition = aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
    _text->render(stringStream2.str(), labelPosition, 0.4f, aergia::COLOR_BLACK);
  }
  GL_DRAW_POINTS(1.f,
                 for (unsigned i = 0; i < _domain->cellCount(); ++i) {

                   aergia::glColor(cellCenterColor);
                   if (static_cast<int>(i) == _selectedCell)
                     aergia::glColor(aergia::COLOR_PURPLE);
                   auto p = _domain->cellCenterPosition(i);
                   aergia::glVertex(Ponos::point2(p));
                 });
  if (_selectedCell >= 0) {
    auto r = _domain->cellRegion(static_cast<unsigned int>(_selectedCell));
    aergia::draw_bbox(Ponos::bbox2D(r), aergia::COLOR_PURPLE,
                      aergia::COLOR_TRANSPARENT);
  }
  if (drawVelocities) {
    aergia::glColor(faceCenterColor);
    _domain->iterateDxScalar(0, [](double &v, Point2d p) {
      aergia::draw_vector(Ponos::point2(p), ponos::Vector2(v, 0) * 0.05);
    });
    _domain->iterateDyScalar(0, [](double &v, Point2d p) {
      aergia::draw_vector(Ponos::point2(p), ponos::Vector2(0, v) * 0.05);
    });
  }
  if (drawPressure) {
    _domain->iterateCellScalar(0, [&](size_t id, double &v) {
      double t = 0.;
      if (maxPressure - minPressure > 0)
        t = (v - minPressure) / (maxPressure - minPressure);
      aergia::draw_bbox(Ponos::bbox2D(_domain->cellRegion(id)),
                        aergia::COLOR_TRANSPARENT,
                        pressurePalette(static_cast<float>(t), 0.6f));
    });
  }
  if (drawDivergence) {
    for (unsigned id = 0; id < _domain->cellCount(); ++id) {
      double v = divergence[id];
      double t = 0.;
      if (maxDivergence - minDivergence > 0)
        t = (v - minDivergence) / (maxDivergence - minDivergence);
      aergia::draw_bbox(Ponos::bbox2D(_domain->cellRegion(id)),
                        aergia::COLOR_TRANSPARENT,
                        pressurePalette(static_cast<float>(t), 0.6f));
    }
  }
  for (unsigned i = 0; i < _domain->faceCount(); i++) {
    auto c = aergia::Color(0.f, 0.f, 0.f, 0.1);
    if (static_cast<int>(i) == _selectedFace)
      c = aergia::Color(1.f, 0.f, 0.f, 0.3);
    aergia::draw_bbox(Ponos::bbox2D(BBox2D(_domain->faceCenterPosition(i), 0.0005)),
                      aergia::COLOR_TRANSPARENT, c);
  }
}

bool SimulationRegularDomain2Model::intersect(const ponos::Ray3 &r, float *t) {
  UNUSED_VARIABLE(t);
  Point2d p(r.o.x, r.o.y);
  _selectedCell = -1;
  _selectedFace = -1;
  if (_domain->region().inside(p))
    _selectedCell = _domain->cellId(p);
  _faces->search(BBox2D(p, 0.0005), [&](size_t i) {
    _selectedFace = static_cast<int>(i);
  });
  return _selectedCell >= 0;
}

void SimulationRegularDomain2Model::updatePressure() {
  minPressure = INFINITY;
  maxPressure = -(INFINITY - 1);
  _domain->iterateCellScalar(0, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    minPressure = std::min(minPressure, v);
    maxPressure = std::max(maxPressure, v);
  });
}

void SimulationRegularDomain2Model::updateDivergence() {
  minDivergence = INFINITY;
  maxDivergence = -(INFINITY - 1);

  for (unsigned i = 0; i < _domain->cellCount(); i++) {
    divergence[i] = _domain->faceDivergentAtCellCenter(0, i);
    minDivergence = std::min(minDivergence, divergence[i]);
    maxDivergence = std::max(maxDivergence, divergence[i]);
  }
}

}  // namespace furoo
