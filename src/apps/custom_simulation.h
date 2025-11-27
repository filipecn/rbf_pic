#ifndef FUROO_CUSTOM_SIMULATION_H
#define FUROO_CUSTOM_SIMULATION_H

#include <apps/simulation_steps.h>
#include <blas/temporal_integrator.h>
#include <common/args_parser.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <physics/injector.h>
#include <physics/particle_system.h>
#include <physics/solid_box.h>
#include <string>
#include <structures/structure_interface.h>
#include <vector>
#include <omp.h>

namespace furoo {

/// Convenient class for setting up simulations with different algorithms built
/// from the set of simulation step class methods or user defined ones.
/// This class contains a cellgraph structure and several fields added to it.
/// It also holds a particle system with velocity properties.
template <int D> class CustomSimulation {
public:
  CustomSimulation();

  /// Loads simulation parameters and algorithm from a file
  /// The configuration file accepts the following variables:
  /// - size_t minLevel : min tree level (5)
  /// - size_t maxLevel : max tree level (7)
  /// - string "outputDir" : output directory ("output")
  /// - double "particleSpacing" : particles spacing (0.00390625)
  /// - int "frames" : total number of frames (60)
  /// - double "dt" : time step (0.001)
  /// - string "scene" : scene can be LDB, RDB, DDB and DROP ("DROP")
  /// - stringList "algorithm" : list of steps
  /// - doubleList "fluidBoxes" : list of boxes [BOX1min BOX1max ... BOXNmin
  /// BOXNmax]
  /// - doubleList "solidBoxes" : list of boxes [BOX1min BOX1max ... BOXNmin
  /// BOXNmax]
  /// - intList "iFluidBoxes" : [level^D BOX1min BOX1max ... BOXNmin BOXNmax]
  /// - intList "iSolidBoxes" : [level^D BOX1min BOX1max ... BOXNmin BOXNmax]
  /// - bool "graded" : make tree graded (true)
  /// - string "rbfKernel" : rbf kernel type (GAUSS, PHS)
  /// - string "base" : rbf kernel base type (LINEAR, QUADRATIC)
  /// - double "rbfRadius" : constant rbf radius (0.1)
  /// - bool "adaptiveRBFRadius" : adapt rbf radius based on stencil (false)
  /// \param filename path/to/file
  void loadFromFile(const std::string &filename);

  /// Stores simulation parameters into the file <filename>
  /// The configuration file is written as:
  /// - size_t minLevel : min tree level (5)
  /// - size_t maxLevel : max tree level (7)
  /// - string "outputDir" : output directory ("output")
  /// - double "particleSpacing" : particles spacing (0.00390625)
  /// - int "frames" : total number of frames (60)
  /// - double "dt" : time step (0.001)
  /// - string "scene" : scene can be LDB, RDB, DDB and DROP ("DROP")
  /// - stringList "algorithm" : list of steps
  /// - doubleList "fluidBoxes" : list of boxes [BOX1min BOX1max ... BOXNmin
  /// BOXNmax]
  /// - doubleList "solidBoxes" : list of boxes [BOX1min BOX1max ... BOXNmin
  /// BOXNmax]
  /// - intList "iFluidBoxes" : [level^D BOX1min BOX1max ... BOXNmin BOXNmax]
  /// - intList "iSolidBoxes" : [level^D BOX1min BOX1max ... BOXNmin BOXNmax]
  /// - bool "graded" : make tree graded (true)
  /// - string "rbfKernel" : rbf kernel type (GAUSS, PHS)
  /// - string "base" : rbf kernel base type (LINEAR, QUADRATIC)
  /// - double "rbfRadius" : constant rbf radius (0.1)
  /// - bool "adaptiveRBFRadius" : adapt rbf radius based on stencil (false)
  /// \param filename path/to/file
  void saveConfigFile(const std::string &filename);

  /// Injects particles based on a pre-defined scene. Available scenes are:
  /// "DROP" :  a drop falling into a tank of water
  /// "LDB" : left dam break
  /// "RDB" : right dam break
  /// "DDB" : double dam break
  /// \param scene scene identifier
  /// \param particleSpacing space between particles
  void loadScene(const std::string &scene, double particleSpacing);

  /// Injects particles based on a list of box regions described as
  /// BOX1min BOX1max ... BOXNmin BOXNmax
  /// where it point is a list of D doubles representing the coordinates
  /// \param boxes box list
  /// \param particleSpacing space between particles
  void loadScene(std::vector<double> boxes, double particleSpacing);

  /// Injects particles based on a list of box regions described as
  /// level^D BOX1min BOX1max ... BOXNmin BOXNmax
  /// where it point is a list of D doubles representing the coordinates
  /// \param boxes box list
  /// \param particleSpacing space between particles
  void loadScene(std::vector<int> boxes, std::vector<int> spheres,
                 double particleSpacing);

  /// Injects particles based on a list of box regions described as
  /// level^D BOX1min BOX1max ... BOXNmin BOXNmax
  /// where it point is a list of D doubles representing the coordinates
  /// \param boxes box list
  /// \param particleSpacing space between particles
  void loadImplicitScene(std::vector<int> boxes, double particleSpacing);

  /// Injects particles based on a list of box regions described as
  /// BOX1min BOX1max ... BOXNmin BOXNmax
  /// where it point is a list of D doubles representing the coordinates
  /// \param boxes box list
  void loadSceneSolids(std::vector<double> boxes);

  /// Injects particles based on a list of box regions described as
  /// level^D BOX1min BOX1max ... BOXNmin BOXNmax
  /// where it point is a list of D doubles representing the coordinates
  /// \param boxes box list
  void loadSceneSolids(std::vector<int> boxes);

  /// Loads a list of injectors described as
  /// Size^D Angle^(D-1)(degrees) Velocity^D InflowRate
  /// \param injectors list of injectors data
  void loadInjectors(std::vector<double> injectors);

  /// Load list of a given injector index, given as
  /// <Injector index> <Starting time> <End time>
  /// End time can be -1 for endless injection
  /// \param injectors
  void loadInjectorTime(std::vector<double> injectorTimes);

  /// Loads data from a input directory
  /// \param inputDirectory input directory
  /// \param frame frame index
  void loadFrame(std::string inputDirectory, size_t frame);

  /// This callback will be called every at each frame begin in order to compute
  /// the correct dt that will be passed to the steps
  /// \param cflCallback expected to return dx / (max vel)
  void setCFLCallback(const std::function<double()> &cflCallback);

  /// Add builtin step, can be one of these:
  /// "advectParticles"
  /// "updateTree"
  /// "markBoundariesOnFaces"
  /// "particlesToCells"
  /// "particlesToFaces"
  /// "addGravityToCells"
  /// "addGravityToFaces"
  /// "divergenceOnCellsFromFacesAndCellCenter"
  /// "pressureOnCells"
  /// "cellPressureToFacePressure"
  /// "pressureGradientOnFacesFromCells"
  /// "correctFaceVelocitiesFromFacePressureGradient"
  /// "reseedParticles"
  /// "facesToParticles"
  /// "saveAllData"
  /// \param name step string identifier
  void addStep(std::string name);

  /// Steps will be called in the order they were added
  /// \param name step name identifier
  /// \param callBack recieves dt as input parameter
  void addStep(std::string name, std::function<void(double)> callBack);

  /// \param dt frame time step
  void setTimeStep(double dt);

  /// \return frame time step
  double timeStep() const;

  /// Sets current frame, current subframe and curstep to 0
  /// NOTE: it does not change the state of the structure or particles, only the
  /// frame number
  void reset();

  /// Start/Stop simulating
  void toggle();

  /// \return  true if is currently simulating
  bool isRunning() const;

  /// run entire simulation
  void run();

  void run(const std::string &stepName);

  /// Executes the next step of the simulation and updates
  /// the frame number and cfl when necessary
  /// \return id of the next step
  size_t advanceFrame();

  /// Executes the next step of the simulation and updates
  /// the frame number and cfl when necessary
  /// \return id of the next step
  size_t advanceStep();

  /// Total number of steps added
  /// \return number of steps
  size_t stepsCount() const;

  /// Steps string identifier
  /// \param i step index
  /// \return step identifier
  std::string stepName(size_t i) const;

  /// The very last completed step
  /// \return current step (to be run)
  size_t curStep() const;

  /// The current frame being computed
  /// \return current frame (1 frame = all steps)
  size_t curFrame() const;

  /// Set currents frame number (resets the other fields)
  /// \param frame frame number
  void setCurFrameNumber(size_t frame);

  /// As cfl condition brakes dt, one frame may be simulated
  /// by multiple subframes of small pieces of dt
  /// \return current frame
  size_t curSubFrame() const;

  ParticleSystem<double, D> *particles();
  StructureInterface<D> *cellGraph();
  std::vector<std::shared_ptr<SolidInterface<D>>> &solids();
  int fieldId(const std::string &fieldName) const;
  const std::vector<std::string> &fieldNames() const;
  const std::vector<std::string> &cellFieldNames() const;
  const std::vector<std::string> &faceFieldNames() const;
  const std::vector<std::string> &vertexFieldNames() const;
  const std::vector<std::string> &particleFieldNames() const;

  /// Provides a list of strings representing the names of the steps in the
  /// current algorithm
  /// \param names [output] list of names
  void stepNames(std::vector<std::string> &names) const;

  size_t minLevel() const { return _minLevel; }
  size_t maxLevel() const { return _maxLevel; }
  std::string outputPath() const { return _outputPath; }
  void setOutputPath(std::string path) { _outputPath = path; }

private:
  void setupFields();
  std::map<std::string, size_t>
      _fieldIds; //!< map of field names to the respective ids in the structure
  std::vector<std::string> _fieldNames;   //!< list of structure field names
  std::vector<std::string> _pFieldsNames; //!< particle field names
  std::vector<std::string> _cFieldsNames; //!< cell centered field names
  std::vector<std::string> _vFieldsNames; //!< vertex centered field names
  std::vector<std::string> _fFieldsNames; //!< face centered field names
  std::vector<std::shared_ptr<SolidInterface<D>>> _solids; //!< list of solids
  std::vector<ParticleInjector<D>> _injectors; //!< list of injectors
  std::vector<Definitions::Boundary>
      _domainBoundaries; //!< boundary types of domain limits
  //   MLS<D> _interpolator = MLS<D>(Definitions::PolynomialType::QUADRATIC);
  RBFInterpolant<Point<double, D>> _interpolator =
    RBFInterpolant<Point<double, D>>{};
  EulerTemporalIntegrator<D> _integrator;
  std::shared_ptr<RadialKernel<Point<double, D>>> _rbfKernel;
  std::shared_ptr<DifferentialRBF<D>> _rbf;
  size_t _maxLevel = 5, _minLevel = 7;
  bool _adaptiveRBFRadius = false, _graded = true;
  double _rbfRadius = 0.1;
  std::string _outputPath = "";
  std::shared_ptr<StructureInterface<D>> _graph;
  std::shared_ptr<ParticleSystem<double, D>> _particles;
  std::map<std::string, std::function<void(double)>>
      _stepMap;               //!< map string function names to their definition
  bool _isSimulating = false; //!< simulation status
  size_t _frames = 1000;
  size_t _curSubFrame =
      0; //!< current sub frame (in case cfl breaks dt into pieces)
  size_t _curFrame = 0; //!< current frame
  size_t _curStep = 0;  //!< current step
  size_t _cflSteps = 0; //!< number of sub dts
  double _cflDT = 0.;   //!< dt defined by cfl
  double _dt = 0.001;   //!< time step
  std::function<double()>
      _cflCallback; //!< extepcted to return dx / (max vel)inl
  std::function<bool(size_t)>
      _useVerticalFluidFace; //!< filter for vertical fluid faces
  std::function<bool(size_t)>
      _useHorizontalFluidFace; //!< filter for horizontal fluid faces
  std::function<bool(size_t)>
      _useDepthFluidFace;                    //!< filter for depth fluid faces
  std::function<bool(size_t)> _useFluidFace; //!< filter for fluid faces
  std::function<bool(size_t)> _useFluidCell; //!< filter for fluid cells
  std::function<bool(size_t)> _useVelocityFace,
      _useVelocityCell; //!< filter for fluid/solid cells
  std::vector<std::pair<std::string, std::function<void(double)>>>
      _steps; //!< steps list
};

#include "custom_simulation.inl"

} // namespace furoo

#endif // FUROO_CUSTOM_SIMULATION_H
