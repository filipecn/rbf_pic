#include <furoo.h>
#include <gtest/gtest.h>

using namespace furoo;

TEST(ParticleSystem, IterateAll) {
  {
    auto domainRegion = BBox2d::squareBox();
    AdaptiveCenteredPic2 simulation(domainRegion, 6, 6);
    auto particles = simulation.particles();
    {
      Injector2 injector(particles, 0.007);
      injector.setupBoxShape(BBox2d(Point2d(0.3, 0.7), Point2d(0.7, 0.8)));
      size_t total = 0;
      std::set<size_t> ids;
      particles->iterateParticles([&](size_t id, Point2d pos) {
        UNUSED_VARIABLE(pos);
        total++;
        if (ids.find(id) != ids.end())
          std::cerr << "aosifhasfh " << id << std::endl;
        ids.insert(id);
        particles->setScalarProperty(1, id, -2 * 0.4);
      });
      std::cerr << ids.size() << "total " << total << std::endl;
    }
  }
  {
    auto region = BBox2d::squareBox();
    ParticleSystem2d ps(new PointZGrid2d(64));
    ps.addScalarProperty(0.);
    ps.addScalarProperty(0.);
    ps.setDomainRegion(region);
    Injector2 injector(&ps, 0.007);
    size_t total =
        injector.setupBoxShape(BBox2d(Point2d(0.3, 0.7), Point2d(0.7, 0.8)));
    std::cerr << total << std::endl;
    std::set<size_t> ids;
    ps.iterateParticles([&](size_t id, Point2d p) {
      UNUSED_VARIABLE(p);
      ids.insert(id);
    });
    EXPECT_EQ(total, ids.size());
  }
}

TEST(ParticleSystem, Shepard) {
  {
    auto region = BBox2d::squareBox();
    ParticleSystem2d ps(new PointZGrid2d(64)/*, new MLS<2>()*/);
    ps.setDomainRegion(region);
    size_t vx = ps.addScalarProperty(0.);
    size_t vy = ps.addScalarProperty(0.);
    BoxSampler sampler;
    for (size_t i = 0; i < 10; i++) {
      size_t id = ps.addParticle(sampler.sample(region));
      EXPECT_EQ(id, i);
      ps.setScalarProperty(vx, id, 1. * id);
      ps.setScalarProperty(vy, id, -1. * id);
    }
    for (size_t i = 0; i < 10; i++) {
      ASSERT_NEAR(ps.getScalarProperty(vx, i), 1. * i, 1e-8);
      ASSERT_NEAR(ps.getScalarProperty(vy, i), -1. * i, 1e-8);
      ps.setScalarProperty(vx, i, std::sin(ps[i].x()));
      ps.setScalarProperty(vy, i, std::cos(ps[i].y()));
    }
    for (size_t i = 0; i < 1000; i++) {
      auto p = sampler.sample(region);
      ASSERT_NEAR(ps.sampleScalarProperty(vx, 0.05, p),
                  std::sin(p.x()),
                  2e-2);
      ASSERT_NEAR(ps.sampleScalarProperty(vy, 0.05, p),
                  std::cos(p.y()),
                  2e-2);
    }
  }
}

TEST(ParticleSystem2, Constructors) {
  {
    auto region = BBox2D::make_unit_bbox();
    ParticleSystem2 ps(new PointZGrid2(64, 1. / 64.),
                       new RBFInterpolant<Point2d>());
    ps.setDomainRegion(region);
    size_t vx = ps.addScalarProperty(0.);
    size_t vy = ps.addScalarProperty(0.);
    BoxSampler sampler;
    for (size_t i = 0; i < 10000; i++) {
      size_t id = ps.addParticle(sampler.sample(region));
      EXPECT_EQ(id, i);
      ps.setScalarProperty(vx, id, 1. * id);
      ps.setScalarProperty(vy, id, -1. * id);
    }
    for (size_t i = 0; i < 10000; i++) {
      ASSERT_NEAR(ps.getScalarProperty(vx, i), 1. * i, 1e-8);
      ASSERT_NEAR(ps.getScalarProperty(vy, i), -1. * i, 1e-8);
      ps.setScalarProperty(vx, i, std::sin(ps[i].x()));
      ps.setScalarProperty(vy, i, std::cos(ps[i].y()));
    }
    for (size_t i = 0; i < 1000; i++) {
      auto p = sampler.sample(region);
      ASSERT_NEAR(ps.sampleScalarProperty(vx, 0.05, p),
                  std::sin(p.x()),
                  1e-3);
      ASSERT_NEAR(ps.sampleScalarProperty(vy, 0.05, p),
                  std::cos(p.y()),
                  1e-3);
    }
  }
}

TEST(ParticleSystem2, Shepard) {
  // return;
  {
    auto region = BBox2D::make_unit_bbox();
    ParticleSystem2 ps(new PointZGrid2(64, 1. / 64.),
                       new ShepardInterpolant2());
    ps.setDomainRegion(region);
    size_t vx = ps.addScalarProperty(0.);
    size_t vy = ps.addScalarProperty(0.);
    BoxSampler sampler;
    for (size_t i = 0; i < 10; i++) {
      size_t id = ps.addParticle(sampler.sample(region));
      EXPECT_EQ(id, i);
      ps.setScalarProperty(vx, id, 1. * id);
      ps.setScalarProperty(vy, id, -1. * id);
    }
    for (size_t i = 0; i < 10; i++) {
      ASSERT_NEAR(ps.getScalarProperty(vx, i), 1. * i, 1e-8);
      ASSERT_NEAR(ps.getScalarProperty(vy, i), -1. * i, 1e-8);
      ps.setScalarProperty(vx, i, std::sin(ps[i].x()));
      ps.setScalarProperty(vy, i, std::cos(ps[i].y()));
    }
    for (size_t i = 0; i < 1000; i++) {
      auto p = sampler.sample(region);
      ASSERT_NEAR(ps.sampleScalarProperty(vx, 0.05, p),
                  std::sin(p.x()),
                  2e-2);
      ASSERT_NEAR(ps.sampleScalarProperty(vy, 0.05, p),
                  std::cos(p.y()),
                  2e-2);
    }
  }
}

TEST(ParticleSystem2, TrilinearHatKernel) {
  return;
  /* {
     auto region = BBox2D::make_unit_bbox();
     ParticleSystem2 ps(
         new PointZGrid2(64, 1. / 64.),
         new KernelInterpolant<Point2d>(Point2d(1 / 64., 1 / 64.)));
     ps.setDomainRegion(region);
     size_t vx = ps.addScalarProperty(0.);
     size_t vy = ps.addScalarProperty(0.);
     BoxSampler sampler;
     for (size_t i = 0; i < 10000; i++) {
       size_t id = ps.addParticle(sampler.sample(region));
       ps.setScalarProperty(vx, id, 1. * id);
       ps.setScalarProperty(vy, id, -1. * id);
     }
     for (size_t i = 0; i < 10000; i++) {
       ps.setScalarProperty(vx, i, std::sin(ps[i].x()));
       ps.setScalarProperty(vy, i, std::cos(ps[i].y()));
     }
     for (size_t i = 0; i < 1000; i++)
       std::cerr << particles->getScalarProperty(0, id) << " " <<
                 enright(p).x() << " | ";
     {
       std::cout << i << '\n';
       auto p = sampler.sample(region);
       ASSERT_NEAR(ps.sampleScalarProperty(vy, 1 / 64., p), std::cos(p.y()),
                   1e-2);
       ASSERT_NEAR(ps.sampleScalarProperty(vx, 1 / 64., p), std::sin(p.x()),
                   1e-2);
     }
   }*/
}
