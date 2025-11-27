#include <furoo.h>

class Injector {
 public:
  explicit Injector(furoo::ParticleSystem2 *ps, double s = 0.03)
      : _spacing(s), _particleSystem(ps) {}
  void setSpacing(double s);
  void setupBoxShape(furoo::Point2d pMin, furoo::Point2d pMax);
  void setupCircleShape(furoo::Point2d center, double radius);
  void loadFromFile(const char *filename);

 private:
  double _spacing;
  furoo::ParticleSystem2 *_particleSystem;
};

class PICSolverDebuggerClass : public furoo::PICSolver2 {
 public:
  PICSolverDebuggerClass() = default;
  PICSolverDebuggerClass(furoo::DomainBoundaryHandler2 *b,
                         furoo::SimulationDomain2 *d,
                         furoo::ParticleSystem2 *ps)
      : PICSolver2(b, d, ps) {}

  void runOnInitialize() { this->onInitialize(); }

  void runAdvanceFrame() { this->advanceFrame(); }

  void runApplyBoundaryCondition() { this->applyBoundaryCondition(); }

  void runBeginAdvanceTimeStep(double dt) { this->beginAdvanceTimeStep(dt); }

  void runAdvanceTimeStep(double dt) { this->onAdvanceTimeStep(dt); }

  void runComputeExternalForces(double dt) { this->computeExternalForces(dt); }

  void runComputeViscosity(double dt) { this->computeViscosity(dt); }

  void runComputePressure(double dt) { this->computePressure(dt); }

  void runCorrectVelocities(double dt) { this->correctVelocities(dt); }

  void runTransferFromParticlesToGrid() { this->transferFromParticlesToGrid(); }

  void runTransferFromGridToParticles() { this->transferFromGridToParticles(); }

  void runComputeAdvection(double dt) { this->computeAdvection(dt); }

  void runEndAdvanceTimeStep(double dt) { this->endAdvanceTimeStep(dt); }

  void runExtrapolateToDomain() { this->propagateVelocitiesToDomain(); }

  void setParticleSystemPtr(furoo::ParticleSystem2 *p) { this->_particles = p; }

  int getNumberOfSubTimeSteps(double dt) {
    return this->numberOfSubTimeSteps(dt);
  }
};

void saveFramebuffer(const char *filename, unsigned int WIDTH,
                     unsigned int HEIGHT);

furoo::Point2d
advectExact(furoo::Point2d p,
            double dt,
            std::function<furoo::Vector2d(furoo::Point2d)> velocity,
            unsigned int subStepsCount = 2);
