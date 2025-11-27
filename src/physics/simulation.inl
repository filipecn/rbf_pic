template<int D>
Simulation<D>::Simulation(StructureInterface <D> *structure,
                          ParticleSystem<double, D> *particles,
                          TemporalIntegratorInterface <D> *integrator,
                          Solver <D> *solver)
    : _structure(structure),
      _particles(particles),
      _integrator(integrator),
      _solver(solver) {

}

template<int D>
Simulation<D>::~Simulation() {

}

template<int D>
ParticleSystem<double, D> *Simulation<D>::particles() {
  return _particles.get();
}

template<int D>
void Simulation<D>::setRegion(BBox<double, D> region) {
  if (_particles)
    _particles->setDomainRegion(region);
}

template<int D>
StructureInterface <D> *Simulation<D>::structure() {
  return _structure.get();
}

template<int D>
void Simulation<D>::setIntegrator(TemporalIntegratorInterface <D> *integrator) {
  _integrator = integrator;
}
template<int D>
void Simulation<D>::setSolver(Solver <D> *solver) {
  _solver = solver;
}
