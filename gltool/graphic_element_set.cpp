#include "graphic_element_set.h"

GraphicElementSet::GraphicElementSet(ponos::RawMesh *rawMesh) {
  _sceneMesh.reset(new aergia::SceneMesh(*rawMesh));
  setShader();
}

GraphicElementSet::GraphicElementSet(ponos::RawMeshSPtr &rawMesh) {
  _rawMesh = rawMesh;
  _sceneMesh.reset(new aergia::SceneMesh(*_rawMesh.get()));
  setShader();
}

void GraphicElementSet::setShader() {
  // initialize instances set
  const char *vs =
      "#version 440 core\n"
      // regular vertex attributes
      "layout (location = 0) in vec3 position;"
      "layout (location = 1) in vec2 texcoord;"
      // per instance attributes
      //                   "layout (location = 2) in vec3 pos;"  // instance
      //                   position "layout (location = 3) in float rad;" //
      //                   instance radius
      "layout (location = 2) in mat4 transform;" // instance color
      "layout (location = 6) in vec4 col;"       // instance color
      "layout (location = 7) in int isActive;"   // instance active
      // constants across draw
      "layout (location = 8) uniform mat4 view_matrix;"
      "layout (location = 9) uniform mat4 projection_matrix;"
      // output to fragment shader
      "out VERTEX {"
      "vec4 color;"
      "vec3 normal;"
      "vec2 uv;"
      "flat int isActive;"
      "} vertex;"
      "void main() {"
      // define model_matrix: instance transform
      "    mat4 model_matrix = transform;"
      //                   "    mat4 model_matrix;"
      //                   "    model_matrix[0] = vec4(rad, 0, 0, 0);"
      //                   "    model_matrix[1] = vec4(0, rad, 0, 0);"
      //                   "    model_matrix[2] = vec4(0, 0, rad, 0);"
      //                   "    model_matrix[3] = vec4(pos.x, pos.y, pos.z, 1);"
      "    mat4 model_view_matrix = view_matrix * model_matrix;\n"
      "    gl_Position = projection_matrix * model_view_matrix * "
      "vec4(position,1);"
      "   vertex.color = col;"
      "   vertex.isActive = isActive;"
      "}";
  const char *fs = "#version 440 core\n"
                   "in VERTEX { vec4 color; vec3 normal; vec2 uv; flat int "
                   "isActive; } vertex;"
                   "out vec4 outColor;"
                   "void main() {"
                   " if(vertex.isActive == 0) discard;"
                   "   outColor = vertex.color;"
                   "}";
  aergia::ShaderProgram shader(vs, nullptr, fs);
  shader.addVertexAttribute("position", 0);
  shader.addVertexAttribute("texcoord", 1);
  //  shader.addVertexAttribute("pos", 2);
  //  shader.addVertexAttribute("rad", 3);
  shader.addVertexAttribute("transform", 2);
  shader.addVertexAttribute("col", 6);
  shader.addVertexAttribute("isActive", 7);
  shader.addUniform("view_matrix", 8);
  shader.addUniform("projection_matrix", 9);
  // init instances
  _instances.reset(
      new aergia::InstanceSet(*_sceneMesh.get(), shader, _instancesCount));
  aergia::BufferDescriptor trans = aergia::create_array_stream_descriptor(16);
  trans.addAttribute("transform", 16, 0, trans.dataType);
  _posBuffer = _instances->add(trans);
  aergia::BufferDescriptor col = aergia::create_array_stream_descriptor(4);
  col.addAttribute("col", 4 /* r g b a */, 0, col.dataType);
  _colorBuffer = _instances->add(col);
  aergia::BufferDescriptor act =
      aergia::create_array_stream_descriptor(1, GL_UNSIGNED_INT);
  act.addAttribute("isActive", 1, 0, act.dataType);
  _activeBuffer = _activeBuffer = _instances->add(act);
}

void GraphicElementSet::draw(const aergia::CameraInterface *camera,
                             ponos::Transform t) {
  if (_instances && _instancesCount)
    _instances->draw(camera, t);
}

bool GraphicElementSet::intersect(const ponos::Ray3 &r, float *t) {
  return SceneObject::intersect(r, t);
}

size_t GraphicElementSet::count() const { return _instancesCount; }

void GraphicElementSet::resize(size_t n) {
  _instancesCount = n;
  _instances->resize(_instancesCount);
}

uint *GraphicElementSet::activeBuffer(size_t i) {
  return _instances->instanceU(_activeBuffer, i);
}
float *GraphicElementSet::colorBuffer(size_t i) {
  return _instances->instanceF(_colorBuffer, i);
}
float *GraphicElementSet::transformBuffer(size_t i) {
  return _instances->instanceF(_posBuffer, i);
}
