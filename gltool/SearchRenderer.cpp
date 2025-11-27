//
// Created by filipecn on 11/28/17.
//

#include "SearchRenderer.h"
#include "ponos_wrapper.h"

namespace furoo {

SearchRenderer::SearchRenderer(ParticleSystem2 *p, SimQuadTree *q)
    : _ps(p), _domain(q) {
  _selectedParticle = _selectedCell = _selectedFace = -1;
  _faces.reset(new PointZGrid2(16u, BBox2D::make_unit_bbox()));
  for (unsigned i = 0; i < _domain->faceCount(); i++)
    _faces->add(_domain->faceCenterPosition(i));
}

void SearchRenderer::draw() {
  if (_selectedParticle >= 0) {
    _domain->iterateVFaces((*_ps)[_selectedParticle], /*0.01,*/ [&](size_t i) {

      auto c = aergia::Color(0.9f, 0.f, 0.6f, 0.4);
      aergia::draw_bbox(Ponos::bbox2D(BBox2D(_domain->faceCenterPosition(i),
                                             0.0010)),
                        aergia::COLOR_TRANSPARENT, c);
    });
    _domain->iterateUFaces((*_ps)[_selectedParticle], /*0.01,*/ [&](size_t i) {

      auto c = aergia::Color(0.4f, 0.6f, 0.9f, 0.4);
      aergia::draw_bbox(Ponos::bbox2D(BBox2D(_domain->faceCenterPosition(i),
                                             0.0010)),
                        aergia::COLOR_TRANSPARENT, c);
    }, 2);
  }
  if (_selectedCell >= 0 && _selectedParticle < 0 && _selectedFace < 0) {
    auto p = _domain->cellCenterPosition(_selectedCell);
    _domain->iterateVFaces(p, 0.01, [&](size_t i) {
      auto c = aergia::Color(0.9f, 0.f, 0.6f, 0.4);
      aergia::draw_bbox(Ponos::bbox2D(BBox2D(_domain->faceCenterPosition(i),
                                             0.0010)),
                        aergia::COLOR_TRANSPARENT, c);
    });
    _domain->iterateUFaces(p, 0.01, [&](size_t i) {
      auto c = aergia::Color(0.9f, 0.f, 0.6f, 0.4);
      aergia::draw_bbox(Ponos::bbox2D(BBox2D(_domain->faceCenterPosition(i),
                                             0.0010)),
                        aergia::COLOR_TRANSPARENT, c);
    });
  }
  if (_selectedFace >= 0) {
    auto p = _domain->faceCenterPosition(_selectedFace);
    aergia::glColor(aergia::COLOR_BLACK);
    GL_DRAW_LINES(3.,
                  _domain->iterateCells(p, 0.005, [&](size_t id) {
                    aergia::glVertex(Ponos::point2(p));
                    aergia::glVertex(Ponos::point2(_domain->cellCenterPosition(
                        id)));
                  });
    );
  }
  if (_domain->region().inside(mouse)) {
    glColor4f(0, 1, 0, 1);
    auto p = mouse;
    aergia::draw_bbox(ponos::BBox2D(ponos::Point2(p.x() - 0.001, p.y() - 0.001),
                                    ponos::Point2(p.x() + 0.001,
                                                  p.y() + 0.001)),
                      aergia::COLOR_BLACK,
                      aergia::COLOR_GREEN);
    glPointSize(5.);
    GL_DRAW_POINTS(3.,
                   _ps->iterateClosestParticles(p,
                                                10,
                                                [](size_t i, Point2d pa) {
                                                  UNUSED_VARIABLE(i);
                                                  aergia::glVertex(Ponos::point2(
                                                      pa));
                                                });
    );
  }

}

bool SearchRenderer::intersect(const ponos::Ray3 &r, float *t) {
  UNUSED_VARIABLE(t);
  Point2d c(r.o.x, r.o.y);
  mouse = c;
  _selectedParticle = -1;
  _selectedCell = -1;
  if (_domain->region().inside(c))
    _selectedCell = _domain->cellId(c);
  double particleRadius = 0.0005;
  _ps->iterateParticles(BBox2D(c, 5 * particleRadius),
                        [&](unsigned int id, Point2d p) {
                          if (p.distance(c) <= particleRadius * 1.5)
                            _selectedParticle = id;
                        });
  _selectedFace = -1;
  if (_domain->region().inside(c))
    _selectedCell = _domain->cellId(c);
  _faces->search(BBox2D(c, 0.005), [&](size_t i) {
    _selectedFace = static_cast<int>(i);
  });
  return _selectedParticle >= 0;
}
void SearchRenderer::set(ParticleSystem2 *p, SimQuadTree *q) {
  _ps = p;
  _domain = q;
  _selectedParticle = _selectedCell = _selectedFace = -1;
  _faces.reset(new PointZGrid2(16u, BBox2D::make_unit_bbox()));
  for (unsigned i = 0; i < _domain->faceCount(); i++)
    _faces->add(_domain->faceCenterPosition(i));
}

}
