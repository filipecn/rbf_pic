#define _USE_OPENMP // uncomment this line to use OpenMP
//#define _USE_TBB  // uncomment this line to use TBB (future work)
#include <furoo.h>
//#include <iomanip>

/*
// OpenMP simple test
#include <omp.h>
#include <stdio.h>

int
main()
{
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < 20; ++i)
#pragma omp critical
      printf("[%02d]: thread %d\n", i, omp_get_thread_num());
  }
  getchar();
  return 0;
}
*/

int main(int argc, char **argv) {
  if (argc < 3)
    std::cout << "custom_pic <configFile> <Dimensions> OPT[<frameDirectory "
                 "frameindex>]\n";
  int D = atoi(argv[2]);
  std::cerr << "simulation in " << D
            << " dimensions under configuration loaded from " << argv[1]
            << std::endl;
  if (D == 2) {
    furoo::CustomSimulation<2> sim;
    sim.loadFromFile(argv[1]);
    if (argc == 5)
      sim.loadFrame(argv[3], static_cast<size_t>(atoi(argv[4])));
    sim.run();
  } else {
    furoo::CustomSimulation<3> sim;
    sim.loadFromFile(argv[1]);
    if (argc == 5)
      sim.loadFrame(argv[3], static_cast<size_t>(atoi(argv[4])));
    try{
      std::cout << "will do sim.run now" << std::endl;
      sim.run();
    } catch (std::string &e){
      std::cerr << e << std::endl;
      return -1;
    }
  }
  return 0;
}
