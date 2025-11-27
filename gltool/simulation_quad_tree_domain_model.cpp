#include "simulation_quad_tree_domain_model.h"
#include "ponos_wrapper.h"
#include <sstream>

namespace furoo {

SimulationQuadTreeDomainModel::SimulationQuadTreeDomainModel(SimQuadTree *d,
                                                             DomainBoundaryHandler2 *bh)
    : _domain(d), _bh(bh) {
  pressurePalette = aergia::HEAT_GREEN_COLOR_PALETTE;
  divergence.resize(_domain->cellCount(), 0.);
  _selectedCell = _selectedFace = -1;
  _faces.reset(new PointZGrid2(16u, BBox2D::make_unit_bbox()));
  for (unsigned i = 0; i < _domain->faceCount(); i++)
    _faces->add(_domain->faceCenterPosition(i));
  for (unsigned i = 0; i < _domain->cellCount(); i++)
    boxes.emplace_back(_domain->cellRegion(i));
  _drawMaterial = false;
}

void SimulationQuadTreeDomainModel::draw() {
  aergia::TextRenderer text;
  for (unsigned i = 0; i < _domain->cellCount(); i++) {
    if (static_cast<int>(i) == _selectedCell && _drawPressure) {
      std::ostringstream stringStream2;
      stringStream2 << "p: " << _domain->scalarAtCell(0, i);
      auto p = _domain->cellCenterPosition(i);
      auto labelPosition =
          aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
      text.render(stringStream2.str(),
                  labelPosition,
                  0.4f,
                  aergia::COLOR_BLACK);
    }
    if (_bh->cellType(i) == DomainBoundaryHandler2::MaterialType::FLUID) {
      aergia::draw_bbox(Ponos::bbox2D(boxes[i]), aergia::COLOR_BLACK,
                        (_drawMaterial) ? aergia::Color(0, 0, 1, 0.5)
                                        : aergia::COLOR_TRANSPARENT);
    } else {
      aergia::draw_bbox(Ponos::bbox2D(boxes[i]), aergia::COLOR_BLACK,
                        _selectedCell != static_cast<int>(i)
                        ? aergia::COLOR_TRANSPARENT : aergia::Color(1,
                                                                    0,
                                                                    0,
                                                                    0.1));
    }
  }
  for (unsigned i = 0; i < _domain->faceCount(); i++) {
    auto c = aergia::COLOR_TRANSPARENT;
    if (_bh->faceBoundaryType(i)
        == DomainBoundaryHandler2::BoundaryType::NEUMANN)
      c = aergia::Color(0.f, 1.f, 0.f, 1.f);
    else if (_bh->faceBoundaryType(i)
        == DomainBoundaryHandler2::BoundaryType::DIRICHLET)
      c = aergia::Color(1.f, 0.f, 0.f, 1.f);
    auto f = aergia::Color(0, 0, 0, 0.2);
    if (_bh->faceType(i) == DomainBoundaryHandler2::MaterialType::SOLID)
      f = aergia::Color(0, 0, 0, 0.9);
    else if (_bh->faceType(i) == DomainBoundaryHandler2::MaterialType::FLUID)
      f = aergia::Color(0, 1, 1, 0.9);
    if (static_cast<int>(i) == _selectedFace) {
      f = aergia::Color(1.f, 0.f, 0.f, 0.5);
      std::ostringstream stringStream2;
      stringStream2 << "v: " << _domain->scalarAtFace(0, _selectedFace)
                    << " g: "
                    << _domain->cellGradientAtFaceCenter(0, _selectedFace);
      auto p = _domain->faceCenterPosition(_selectedFace);
      auto labelPosition =
          aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
      text.render(stringStream2.str(),
                  labelPosition,
                  0.4f,
                  aergia::COLOR_BLACK);
    }
    glLineWidth(4.);
    aergia::draw_bbox(Ponos::bbox2D(BBox2D(_domain->faceCenterPosition(i),
                                           0.0005)),
                      c, f);
    glLineWidth(1.);
  }
  if (_drawPressure) {
    _domain->iterateCellScalar(0, [&](size_t id, double &v) {
      double t = 0.;
      if (maxPressure - minPressure > 0)
        t = (v - minPressure) / (maxPressure - minPressure);
      aergia::draw_bbox(Ponos::bbox2D(_domain->cellRegion(id)),
                        aergia::COLOR_TRANSPARENT,
                        pressurePalette(static_cast<float>(t), 0.6f));
    });
  }/*
  aergia::glColor(aergia::COLOR_BLACK);
  _domain->iterateDxScalar(0, [](double &v, Point2d p) {
    aergia::draw_vector(Ponos::point2(p), ponos::Vector2(v, 0) * 0.05);
  });
  _domain->iterateDyScalar(0, [](double &v, Point2d p) {
    aergia::draw_vector(Ponos::point2(p), ponos::Vector2(0, v) * 0.05);
  });*/
  if (_drawDivergence) {
    for (unsigned id = 0; id < _domain->cellCount(); ++id) {
      if (_bh->cellType(id) != DomainBoundaryHandler2::MaterialType::FLUID)
        continue;
      double v = divergence[id];
      double t = 0.;
      if (maxDivergence - minDivergence > 0)
        t = (v - minDivergence) / (maxDivergence - minDivergence);
      aergia::draw_bbox(Ponos::bbox2D(_domain->cellRegion(id)),
                        aergia::COLOR_TRANSPARENT,
                        pressurePalette(static_cast<float>(t), 0.6f));
      if (_selectedCell >= 0 && _selectedCell == static_cast<int>(id)
          && _selectedFace < 0) {
        std::ostringstream stringStream2;
        stringStream2 << "d: " << divergence[id] << " p: "
                      << _domain->scalarAtCell(0, id);
        auto p = _domain->cellCenterPosition(_selectedCell);
        auto labelPosition =
            aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
        text.render(stringStream2.str(),
                    labelPosition,
                    0.4f,
                    aergia::COLOR_BLACK);
      }
    }
  }
}

bool SimulationQuadTreeDomainModel::intersect(const ponos::Ray3 &r, float *t) {
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

void SimulationQuadTreeDomainModel::updatePressure() {
  minPressure = INFINITY;
  maxPressure = -(INFINITY - 1);
  _domain->iterateCellScalar(0, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    minPressure = std::min(minPressure, v);
    maxPressure = std::max(maxPressure, v);
  });
}

void SimulationQuadTreeDomainModel::updateDivergence() {
  minDivergence = INFINITY;
  maxDivergence = -(INFINITY - 1);

  for (unsigned i = 0; i < _domain->cellCount(); i++) {
    if (_bh->cellType(i) != DomainBoundaryHandler2::MaterialType::FLUID)
      continue;
    divergence[i] = _domain->faceDivergentAtCellCenter(0, i);
    minDivergence = std::min(minDivergence, divergence[i]);
    maxDivergence = std::max(maxDivergence, divergence[i]);
  }
  // std::cerr << minDivergence <<
  //           " " << maxDivergence << std::endl;
}
void SimulationQuadTreeDomainModel::set(SimQuadTree *d,
                                        DomainBoundaryHandler2 *bh) {
  _domain = d;
  _bh = bh;
  divergence.resize(_domain->cellCount(), 0.);
  _selectedCell = _selectedFace = -1;
  _faces.reset(new PointZGrid2(16u, BBox2D::make_unit_bbox()));
  boxes.clear();
  for (unsigned i = 0; i < _domain->faceCount(); i++)
    _faces->add(_domain->faceCenterPosition(i));
  for (unsigned i = 0; i < _domain->cellCount(); i++)
    boxes.emplace_back(_domain->cellRegion(i));
}
void SimulationQuadTreeDomainModel::drawDivergence() {
  _drawMaterial = false;
  _drawVelocities = false;
  _drawPressure = false;
  _drawDivergence = true;
}
void SimulationQuadTreeDomainModel::drawPressure() {
  _drawMaterial = false;
  _drawVelocities = false;
  _drawPressure = true;
  _drawDivergence = false;
}
void SimulationQuadTreeDomainModel::drawMaterial() {
  _drawMaterial = true;
  _drawVelocities = false;
  _drawPressure = false;
  _drawDivergence = false;
}

}
