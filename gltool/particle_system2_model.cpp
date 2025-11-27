#include "particle_system2_model.h"
#include <sstream>
#include "ponos_wrapper.h"

namespace furoo {

ParticleSystem2Model::ParticleSystem2Model(ParticleSystem2 *p, double r)
    : _ps(p) {
  dt = 0.05;
  particleRadius = r;
  _selectedParticle = -1;
  particleColor = aergia::COLOR_BLUE;
  pallete = aergia::HEAT_MATLAB_PALETTE;
  _minValue = _maxValue = 0;
}

void ParticleSystem2Model::draw() {
  std::ostringstream stringStream;
  stringStream << "Selected particle: " << _selectedParticle;
  aergia::TextRenderer text;
  text.render(stringStream.str(), -0.9, -0.9, 0.4f, aergia::COLOR_BLACK);
  _ps->iterateParticles([&](unsigned int id, Point2d p) {
                          double radius = particleRadius;
                          if (static_cast<int>(id) == _selectedParticle) {
                            aergia::glColor(aergia::COLOR_PURPLE);
                            radius = particleRadius * 2.;
                            glLineWidth(3.);
                            aergia::draw_vector(Ponos::point2(p),
                                                ponos::Vector2(_ps->getScalarProperty(0, id),
                                                               _ps->getScalarProperty(1, id)) * dt);
                            glLineWidth(1.);
                          }
                          ponos::Circle c(Ponos::point2(p), static_cast<float>(radius));
                          aergia::glColor(particleColor);
                          { // particle color
                            Vector2d v(_ps->getScalarProperty(0, id), _ps->getScalarProperty(1, id));
                            double m = v.length2();
                            aergia::glColor(pallete(1.f - ponos::smoothStep(m, _minValue, _maxValue),
                                                    1.f));
                          }
                          aergia::draw_circle(c);
                          if (static_cast<int>(id) == _selectedParticle) {
                            std::ostringstream stringStream2;
                            stringStream2 << "v: " << _ps->getScalarProperty(0, _selectedParticle)
                                          << ", "
                                          << _ps->getScalarProperty(1, _selectedParticle);
                            auto p = (*_ps)[_selectedParticle];
                            auto labelPosition =
                                aergia::glGetMVPTransform()(ponos::Point3(p.x(), p.y(), 0.f));
                            text.render(stringStream2.str(),
                                        labelPosition,
                                        0.4f,
                                        aergia::COLOR_WHITE);
                          }
                        }
  );
}

bool ParticleSystem2Model::intersect(const ponos::Ray3 &r, float *t) {
  UNUSED_VARIABLE(t);
  Point2d c(r.o.x, r.o.y);
  _selectedParticle = -1;
  _ps->iterateParticles(BBox2D(c, 5 * particleRadius),
                        [&](unsigned int id, Point2d p) {
                          if (p.distance(c) <= particleRadius * 1.5)
                            _selectedParticle = id;
                        });
  return _selectedParticle >= 0;
}
void ParticleSystem2Model::update() {
  _ps->iterateParticles([&](unsigned int id, Point2d p) {
    UNUSED_VARIABLE(p);
    Vector2d v(_ps->getScalarProperty(0, id), _ps->getScalarProperty(1, id));
    double m = v.length2();
    _maxValue = std::max(_maxValue, m);
  });
}

}  // namespace furoo
