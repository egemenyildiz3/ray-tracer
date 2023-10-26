#include "interpolate.h"
#include <glm/geometric.hpp>

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // TODO: implement this function.

    // Taken from chapter 4 from Real-Time Collison Detection by Christer Ericson

    glm::vec3 w0 = v1 - v0;
    glm::vec3 w1 = v2 - v0;
    glm::vec3 w2 = p - v0;

    float dot00 = glm::dot(w0, w0);
    float dot10 = glm::dot(w1, w0);
    float dot11 = glm::dot(w1, w1);
    float dot20 = glm::dot(w2, w0);
    float dot21 = glm::dot(w2, w1);

    // Using Cramers rule
    float alpha = (dot20 * dot11 - dot10 * dot21) / (dot00 * dot11 - dot10 * dot10);
    float beta = (dot00 * dot21 - dot20 * dot10) / (dot00 * dot11 - dot10 * dot10);
    float gamma = 1.0f - alpha - beta;

    return glm::vec3(alpha, beta, gamma);
}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    // TODO: implement this function.
    return glm::vec3(0.0);
}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
// TODO: implement this function.
    return glm::vec2(0.0);
}
