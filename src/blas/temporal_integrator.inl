template<size_t D>
Point<double, D>
EulerTemporalIntegrator<D>::integrate(Point<double, D> curr,
                                      Vector<double, D> direction,
                                      double dt) const {
  return curr + dt * direction;
}

template<size_t D>
Point<double, D> RungeKutta2TemporalIntegrator<D>::integrate(
    Point<double, D> curr, Vector<double, D> direction, double dt) const {
  // TODO: [filipecn] is it right?! I think we are creating bigger velocity vectors then we should...
  Vector<double, D> k1 = dt * direction;
  Vector<double, D> k2 = dt * direction + .5 * dt * k1;
  Vector<double, D> k3 = dt * direction + .75 * dt * k2;
  return curr + (2. / 9.) * dt * k1 + (3. / 9.) * dt * k2 + (4. / 9.) * dt * k3;
}
