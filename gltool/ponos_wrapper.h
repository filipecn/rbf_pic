#ifndef FUROO_PONOS_WRAPPER_H
#define FUROO_PONOS_WRAPPER_H

#include <furoo.h>

#include <aergia/aergia.h>
#include <geometry/bbox.h>
#include <ponos/ponos.h>

namespace furoo {

template <int D> class SolidObject : public aergia::SceneMeshObject {
public:
  SolidObject(Mesh<D> *mesh) {
    this->rawMesh = new ponos::RawMesh();
    this->rawMesh->meshDescriptor.elementSize = D;
    this->rawMesh->meshDescriptor.count = mesh->elementCount();
    this->rawMesh->positionDescriptor.elementSize = D;
    this->rawMesh->positionDescriptor.count = mesh->vertexCount();
    for (size_t i = 0; i < mesh->vertexCount(); i++)
      for (int d = 0; d < D; d++)
        this->rawMesh->positions.emplace_back(mesh->vertex(i)[d]);
    for (size_t i = 0; i < mesh->elementCount(); i++)
      for (int d = 0; d < D; d++)
        this->rawMesh->indices.push_back(
            {mesh->element(i)[d], mesh->element(i)[d], mesh->element(i)[d]});
    this->rawMesh->splitIndexData();
    this->rawMesh->buildInterleavedData();
    const char *fs = "#version 440 core\n"
                     "out vec4 outColor;"
                     "void main(){"
                     " outColor = vec4(0,0,0,0.4);}";
    if (D == 2) {
      const char *vs = "#version 440 core\n"
                       "layout (location = 0) in vec2 position;"
                       "layout (location = 0) uniform mat4 mvp;"
                       "void main() {"
                       "gl_Position = mvp * vec4(position, 0, 1);}";
      shader.reset(new aergia::ShaderProgram(vs, nullptr, fs));
    } else {
      const char *vs = "#version 440 core\n"
                       "layout (location = 0) in vec3 position;"
                       "layout (location = 0) uniform mat4 mvp;"
                       "void main() {"
                       "gl_Position = mvp * vec4(position, 1);}";
      shader.reset(new aergia::ShaderProgram(vs, nullptr, fs));
    }
    shader->addVertexAttribute("position", 0);
    shader->addUniform("mvp", 0);
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    aergia::BufferDescriptor vd, id;
    aergia::create_buffer_description_from_mesh(*this->rawMesh, vd, id);
    vb.reset(new aergia::VertexBuffer(&this->rawMesh->interleavedData[0], vd));
    ib.reset(new aergia::IndexBuffer(&this->rawMesh->positionsIndices[0], id));
    vb->locateAttributes(*shader.get());
    //  shader->registerVertexAttributes(vb.get());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    //  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  }
  ~SolidObject() {}
  void draw(const aergia::CameraInterface *camera,
            ponos::Transform transform) override {
    glBindVertexArray(VAO);
    shader->begin();
    shader->setUniform("mvp",
                       ponos::transpose((camera->getProjectionTransform() *
                                         camera->getViewTransform() *
                                         camera->getModelTransform())
                                            .matrix()));
    glLineWidth(4.);
    glDrawElements((D == 2) ? GL_LINES : GL_TRIANGLES,
                   ib->bufferDescriptor.elementSize *
                       ib->bufferDescriptor.elementCount,
                   GL_UNSIGNED_INT, 0);
    shader->end();
    glBindVertexArray(0);
    glLineWidth(1.);
  }

  std::shared_ptr<aergia::ShaderProgram> shader; //!< shader applied on draw

private:
  GLuint VAO;
};

class Ponos {
public:
  static ponos::BBox2D bbox2D(const BBox2D &);
  static ponos::BBox2D bbox2D(const BBox2d &);
  static ponos::Point2 point2(const Point2d &);

  template <int D>
  static void drawRegion(BBox<double, D> region, aergia::Color color) {
    auto edgeColor = aergia::COLOR_BLACK;
    edgeColor.a = 0.9;
    aergia::draw_bbox(bbox2D(region), edgeColor, color);
  }

  template <int D>
  static void drawVector(Point<double, D> p, Vector<double, D> v) {
    aergia::draw_vector(ponos::Point2(p.x(), p.y()), ponos::vec2(v.x(), v.y()));
  }

  template <int D>
  static void drawBall(furoo::Point<double, D> c, double r,
                       aergia::Color color) {
    aergia::glColor(color);
    if (D == 2) {
      ponos::Circle circle(ponos::Point2(c.x(), c.y()), r);
      aergia::draw_circle(circle);
    }
  }
};

} // namespace furoo

#endif // FUROO_PONOS_WRAPPER_H
