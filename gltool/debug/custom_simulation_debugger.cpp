#include <furoo.h>
#include <cell_graph2_model.h>
#include <particle_system_model.h>

using namespace furoo;

#define D 2

std::map<std::string, size_t> modelFieldId;

void connectFields(CellGraphModel<D> &model, CustomSimulation<D> &sim) {
  auto graph = sim.cellGraph();
  {
    auto fieldNames = sim.cellFieldNames();
    for (auto &name : fieldNames) {
      int fieldId = sim.fieldId(name);
      std::string fieldType = graph->fieldDataTypeName(fieldId);
      if (fieldType == "double") {
        auto field = new FieldModel<double, D>(graph->template field<
            double>(fieldId));
        field->structureFieldId = fieldId;
        modelFieldId[name] = model.addCellField(field);
      } else if (fieldType == "int") {
        auto field = new FieldModel<int, D>(graph->template field<
            int>(fieldId));
        field->structureFieldId = fieldId;
        modelFieldId[name] = model.addCellField(field);
      } else if (fieldType == "Definitions::Material") {
        auto field = new FieldModel<Definitions::Material, D>(
            graph->field<Definitions::Material>(fieldId));
        field->structureFieldId = fieldId;
        field->isCategorical = true;
        field->colorMap[Definitions::Material::AIR] =
            aergia::COLOR_TRANSPARENT;
        field->colorMap[Definitions::Material::FLUID] = aergia::COLOR_BLUE;
        field->colorMap[Definitions::Material::SOLID] = aergia::COLOR_BLACK;
        modelFieldId[name] = model.addCellField(field);
      } else if (fieldType == "unsigned_char") {
        auto field = new FieldModel<unsigned char, D>(
            graph->field<unsigned char>(fieldId));
        field->structureFieldId = fieldId;
        field->isCategorical = true;
        field->colorMap[0] = aergia::COLOR_TRANSPARENT;
        field->colorMap[1] = aergia::COLOR_GREEN;
        modelFieldId[name] = model.addCellField(field);
      } else exit(-1);
    }
  }
  {
    auto fieldNames = sim.faceFieldNames();
    for (
      auto &name
        : fieldNames) {
      int fieldId = sim.fieldId(name);
      std::string fieldType = graph->fieldDataTypeName(fieldId);
      if (fieldType == "double") {
        auto field = new FieldModel<double, D>(graph->field<
            double>(fieldId));
        field->structureFieldId = fieldId;
        modelFieldId[name] = model.addFaceField(field);
      } else if (fieldType == "Definitions::Boundary") {
        auto field = new FieldModel<Definitions::Boundary, D>(
            graph->field<Definitions::Boundary>(fieldId));
        field->
            isCategorical = true;
        field->colorMap[Definitions::Boundary::NONE] =
            aergia::COLOR_TRANSPARENT;
        field->colorMap[Definitions::Boundary::NEUMANN] = aergia::COLOR_BLUE;
        field->colorMap[Definitions::Boundary::DIRICHLET] = aergia::COLOR_RED;
        field->structureFieldId = fieldId;
        modelFieldId[name] = model.addFaceField(field);
      } else if (fieldType == "Definitions::Material") {
        auto field = new FieldModel<Definitions::Material, D>(
            graph->field<Definitions::Material>(fieldId));
        field->
            isCategorical = true;
        field->colorMap[Definitions::Material::AIR] =
            aergia::COLOR_TRANSPARENT;
        field->colorMap[Definitions::Material::FLUID] = aergia::COLOR_BLUE;
        field->colorMap[Definitions::Material::SOLID] = aergia::COLOR_BLACK;
        field->structureFieldId = fieldId;
        modelFieldId[name] = model.addFaceField(field);
      } else exit(-1);
    }
  }
  {
    auto fieldNames = sim.vertexFieldNames();
    for (
      auto &name
        : fieldNames) {
      int fieldId = sim.fieldId(name);
      std::string fieldType = graph->fieldDataTypeName(fieldId);
      if (fieldType == "double")
        modelFieldId[name] =
            model.addVertexField(new FieldModel<double, D>(graph->field<
                double>(fieldId)
            ));
      else exit(-1);
    }
  }
}

int main() {
  aergia::SceneApp<> app(2000, 2000, "", D == 3);
  if(D == 2)
    app.addViewport2D(0,0,2000,2000);
  CustomSimulation<D> sim;
  sim.loadFromFile("config3d.txt");
  CellGraphModel<D> graphModel;
  ParticleSystemModel<D> particleSystemModel;
  particleSystemModel.set(sim.particles(), 0.001);
  graphModel.set(sim.cellGraph());
  connectFields(graphModel, sim);
  for (size_t i = 0; i < 7; i++)
    sim.advanceFrame();
  graphModel.updateCellsEdges();
  graphModel.updateCellRegions();
  graphModel.showCellText = true;
  app.keyCallback = [&](int key, int scancode, int action, int mods) {
    if (action == GLFW_RELEASE) {
      if (key == GLFW_KEY_W)
        graphModel.textSize += 0.0001;
      if (key == GLFW_KEY_Q)
        graphModel.textSize -= 0.0001;
      if (key == GLFW_KEY_V)
        graphModel.showVertexText = !graphModel.showVertexText;
      if (key == GLFW_KEY_F)
        graphModel.showFaceText = !graphModel.showFaceText;
      if (key == GLFW_KEY_C)
        graphModel.showCellText = !graphModel.showCellText;
      if (key == GLFW_KEY_S)
        graphModel.selectCellField(modelFieldId["cellSurfaceMaskField"]);
      if (key == GLFW_KEY_D)
        graphModel.selectCellField(modelFieldId["cellSurfaceTDistanceField"]);
      if (key == GLFW_KEY_M) {
        graphModel.selectCellField(modelFieldId["cellMaterialField"]);
        graphModel.selectFaceField(modelFieldId["faceMaterialField"]);
      }
      if (key == GLFW_KEY_P) {
        graphModel.selectCellField(modelFieldId["cellPressureField"]);
        graphModel.selectFaceField(modelFieldId["facePressureField"]);
      }
      if (key == GLFW_KEY_G) {
        graphModel.selectCellField(modelFieldId["cellDivergenceField"]);
        graphModel.selectFaceField(modelFieldId["facePressureGradientXField"]);
      }
      if (key == GLFW_KEY_V)
        graphModel.selectFaceField(modelFieldId["faceXVelocityField"]);
      if (key == GLFW_KEY_SPACE)
        sim.advanceStep();
      if (key == GLFW_KEY_ESCAPE)
        app.exit();
      graphModel.updateCellRegions();
      particleSystemModel.update();
      graphModel.updateCellsEdges();
      graphModel.updateFaceRegions();
    }
  };
  app.scene.add(&graphModel);
  app.scene.add(&particleSystemModel);
  app.run();
  return 0;
}
