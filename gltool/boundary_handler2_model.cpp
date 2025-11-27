#include "boundary_handler2_model.h"
#include "ponos_wrapper.h"

namespace furoo {

DomainBoundaryHandler2Model::DomainBoundaryHandler2Model(SimRegularDomain2 *d,
                                                         DomainBoundaryHandler2 *b)
    : _dm(d),
      _bh(b) {
  drawCells = true;
}

void DomainBoundaryHandler2Model::draw() {
  if (drawCells)
    for (unsigned i = 0; i < _dm->cellCount(); ++i) {
      auto r = _dm->cellRegion(i);
      aergia::Color c = aergia::COLOR_WHITE;
      if (_bh->cellType(i) == DomainBoundaryHandler2::MaterialType::FLUID)
        c = aergia::COLOR_BLUE;
      else if (_bh->cellType(i) == DomainBoundaryHandler2::MaterialType::SOLID)
        c = aergia::COLOR_BLACK;
      c.a = 0.4;
      aergia::draw_bbox(Ponos::bbox2D(r), aergia::COLOR_TRANSPARENT, c);
    }
  auto h = _dm->grid().spacing() / 2.;
  auto uCount = _dm->uFacesCount();
  GL_DRAW_LINES(1.f,
                for (unsigned i = 0; i < _dm->faceCount(); ++i) {
                  auto p = _dm->faceCenterPosition(i);
                  glColor4f(0, 0, 0, 0.2);
                  if (_bh->faceBoundaryType(i)
                      == DomainBoundaryHandler2::BoundaryType::DIRICHLET)
                    glColor4f(1, 0, 0, 0.5);
                  else if (_bh->faceBoundaryType(i)
                      == DomainBoundaryHandler2::BoundaryType::NEUMANN)
                    glColor4f(0, 1, 0, 0.5);
                  aergia::glVertex(ponos::Point2(
                      p.x() - ((i >= uCount) ? h.x() : 0.),
                      p.y() - ((i < uCount) ? h.y() : 0.)
                  ));
                  aergia::glVertex(ponos::Point2(
                      p.x() + ((i >= uCount) ? h.x() : 0.),
                      p.y() + ((i < uCount) ? h.y() : 0.)
                  ));
                }

  );
}

} // furoo namespace
