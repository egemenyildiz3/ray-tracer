#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include <texture.cpp>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // ...
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        Image map = state.scene.envMap;
        
        glm::vec3 dir = glm::normalize(ray.direction);
        glm::vec2 coords { std::fmod(1.0f, dir [0]), std::fmod(1.0f, dir[1]) };

        return sampleTextureNearest(map, coords);

        return glm::vec3 { 0, 1, 0 };
    } else {
        return glm::vec3(0.f);
    }
}

float areaAabb(AxisAlignedBox aabb) {
    float x = aabb.upper.x - aabb.lower.x;
    float y = aabb.upper.y - aabb.lower.y;
    float z = aabb.upper.z - aabb.lower.z;
    return x * y * z;
}

// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    // sort the primitives


    // initiate variables
    int sizePrim = primitives.size();
    int amountSplits = sizePrim / 4; // #bins = 1/4 of #primitives
    float sizePerBin = (aabb.upper[axis] - aabb.lower[axis]) / (amountSplits + 1);
    float minArea = std::numeric_limits<float>::max(); // a variable to check which area is the best
    float newArea;
    int MinIndex = sizePrim/2; // the index of the divide with the lowers area
    float distOnLine; // the actual dist of the bin divider
    int lastJ = 0; // variable to make shifting through the array easier and faster
    Primitive prim;
    AxisAlignedBox aabbLower;
    AxisAlignedBox aabbUpper;

    // loop through all bins
    for (int i = 0; i < amountSplits; i++) {
        distOnLine = (i + 1) * sizePerBin;
        
        // search for the end of the bin
        for (int j = lastJ; j < sizePrim; j++) {
            prim = primitives[j];
            // if the centroid is over the border, save the last primitive that was not over the border
            if (computePrimitiveCentroid(prim)[axis] > distOnLine) {
                // save the index of the next primitive to test
                lastJ = j;
                if (lastJ == 0) {
                    lastJ = 1;
                    break;
                }
                
                // get the areas of the two bins that you made
                aabbLower = computeSpanAABB(primitives.subspan(0, lastJ));
                aabbUpper = computeSpanAABB(primitives.subspan(lastJ, primitives.size() - lastJ));
                newArea = areaAabb(aabbLower) + areaAabb(aabbUpper);

                // see if it is better than the last bins
                if (newArea < minArea) {
                    minArea = newArea;
                    MinIndex = lastJ;
                }
                break;
            }
        }
    }
    return MinIndex;
}