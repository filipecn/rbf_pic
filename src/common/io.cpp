#include <physics/particle_system.h>
#include <string>
#include "common/io.h"

namespace furoo {

bool IO::saveToKim(const char *filename, const ParticleSystem2 *ps) {
  char finalName[100];
  sprintf(finalName, "%s.xyz", filename);
  FILE *fp = fopen(finalName, "w+");
  if (!fp)
    return false;
  ps->iterateParticles([&](unsigned int i, Point2d p) {
    UNUSED_VARIABLE(i);
    fprintf(fp, "%lf %lf\n", p.x(), p.y());
  });
  fclose(fp);
  return true;
}

bool IO::loadFromCF(const char *filename,
                    SimulationDomain2 *dm,
                    unsigned fid) {
  char finalName[100];
  sprintf(finalName, "%s.cf", filename);
  FILE *fp = fopen(finalName, "r");
  if (!fp)
    return false;
  dm->iterateCellScalar(fid, [&](double &v, Point2d p) {
    UNUSED_VARIABLE(p);
    double _v = 0., x, y;
    fscanf(fp, "%lf %lf %lf\n", &x, &y, &_v);
    v = _v;
  });
  fclose(fp);
  return true;
}

bool IO::saveToCF(const char *filename,
                  const SimulationDomain2 *dm,
                  unsigned fid) {
  char finalName[100];
  sprintf(finalName, "%s.cf", filename);
  FILE *fp = fopen(finalName, "w+");
  if (!fp)
    return false;
  for (unsigned i = 0; i < dm->cellCount(); i++) {
    auto p = dm->cellCenterPosition(i);
    fprintf(fp, "%lf %lf %lf\n", p.x(), p.y(), dm->scalarAtCell(fid, i));
  }
  fclose(fp);
  return true;
}

bool IO::saveToQT(const char *filename, const SimQuadTree *qt) {
  char finalName[100];
  sprintf(finalName, "%s.qt", filename);
  FILE *fp = fopen(finalName, "w+");
  if (!fp)
    return false;
  auto tree = qt->tree();
  tree.traverse([&](const QuadTree::Node &n) -> bool {
    if (!n.isLeaf()) {
      fprintf(fp, "1 ");
      return true;
    }
    fprintf(fp, "0 ");
    return false;
  });
  fclose(fp);
  return true;
}

bool IO::loadFromQT(const char *filename, SimQuadTree *qt) {
  char finalName[100];
  sprintf(finalName, "%s.qt", filename);
  FILE *fp = fopen(finalName, "r");
  if (!fp)
    return false;
  QuadTree
      tree(BBox2D::make_unit_bbox(), [&](const QuadTree::Node &n) -> bool {
    UNUSED_VARIABLE(n);
    int r = 0;
    fscanf(fp, "%d", &r);
    return r == 1;
  });
  fclose(fp);
  qt->setTree(tree);
  return true;
}

} // furoo namespace
