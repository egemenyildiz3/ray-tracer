#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include <texture.cpp>
#include <iostream>

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

int fac(int n)
{
    if (n == 1 || n == 0) {
        return 1;
    }
    return n * fac(n - 1);
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

    Screen bloomy(image.resolution(), true);

    for (int x = 0; x < bloomy.resolution().x; x++) {
        for (int y = 0; y < bloomy.resolution().y; y++) {
            glm::vec3 currentPixel = image.pixels()[image.indexAt(x, y)];
            if (currentPixel.r >= 0.8f || currentPixel.g >= 0.8f || currentPixel.b >= 0.8f) {
                bloomy.setPixel(x, y, currentPixel);
            }
        }
    }

    Screen blurry(bloomy.resolution(), true);
    int wh = 10;
    std::vector<float> detector;

    for (int a = 0; a < wh; a++) {
        float note = fac(wh - 1) / ((float)fac(wh - a - 1) * (float)fac(a));
        detector.push_back(note);
    }

    float add = 0.0f;
    for (const auto& x : detector) {
        add = add + x;
    }

    for (auto& x : detector) {
        x = x / add;
    }

    Screen newPhoto(bloomy.resolution(), true);

    for (int m = 0; m < bloomy.resolution().y; m++) {
        for (int n = 0; n < bloomy.resolution().x; n++) {
            float re = 0.0f;
            float gr = 0.0f;
            float b = 0.0f;

            for (int c = 0; c < wh; c++) {
                int q = m - (wh-1) / 2 + c;
                if (q < image.resolution().x && q >= 0) {
                    glm::vec3 recent = image.pixels()[image.indexAt(q, n)];
                    re = re + recent.r * detector[c];
                    gr = gr + recent.g * detector[c];
                    b = b + recent.b * detector[c];
                }
            }
            newPhoto.setPixel(m, n, glm::vec3(re, gr, b));
        }
    }
    for (int a = 0; a < bloomy.resolution().y; a++) {
        for (int n = 0; n < bloomy.resolution().x; n++) {
            float re = 0.0f;
            float gr = 0.0f;
            float b = 0.0f;

            for (int c = 0; c < wh; c++) {
                int q = b - (wh-1) / 2 + c;
                if (q >= 0 && q < image.resolution().y) {
                    glm::vec3 recent = newPhoto.pixels()[newPhoto.indexAt(a, q)];
                    re = re + recent.r * detector[c];
                    gr = gr + recent.g * detector[c];
                    b = b + recent.b * detector[c];
                }
            }
            blurry.setPixel(a, b, glm::vec3(re, gr, b));
        }
    }

    for (int a = 0; a < image.resolution().y; a++) {
        for (int b = 0; b < image.resolution().x; b++) {
            glm::vec3 newPix = blurry.pixels()[blurry.indexAt(a, b)] + image.pixels()[image.indexAt(a, b)];
            for (int m = 0; m < 3; m++) {
                if (newPix[m] > 1) {
                    newPix[m] = 1;
                }
            }
            image.setPixel(a, b, newPix);
        }
    }
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
    } else {
        return glm::vec3(0.f);
    }
}

float areaAabb(const AxisAlignedBox aabb) {
    float x = aabb.upper.x - aabb.lower.x;
    float y = aabb.upper.y - aabb.lower.y;
    float z = aabb.upper.z - aabb.lower.z;
    return x * y * z;
}

struct Sorter {
    Sorter(int axis) { this->axis = axis; }
    bool operator() (BVH::Primitive a, BVH::Primitive b) {
        return computePrimitiveCentroid(a)[axis] < computePrimitiveCentroid(b)[axis];
    }

    int axis;
};

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

    // sort the primitives on basis of the right axis, low to high, internal sort, like you have to sort the array 

    // initiate variables
    int amountPrim = primitives.size();
    int amountSplits = 9; // amountPrim / 4; // #bins = 1/4 of #primitives
    float sizePerBin = (aabb.upper[axis] - aabb.lower[axis]) / (amountSplits + 1);
    float minArea = std::numeric_limits<float>::max(); // a variable to check which area is the best
    float newArea;
    int MinIndex = amountPrim/2; // the index of the divide with the lowers area
    float distOnLine; // the actual dist of the bin divider
    int lastJ = 0; // variable to make shifting through the array easier and faster
    Primitive prim;
    AxisAlignedBox aabbLower;
    AxisAlignedBox aabbUpper;

    std::sort(primitives.begin(), primitives.end(), Sorter(axis));
    //std::binary(primitives.begin(), primitives.end(), Sorter(axis));

    // loop through all bins
    for (int i = 0; i < amountSplits; i++) {
        distOnLine = aabb.lower[axis] + (i + 1) * sizePerBin;
        
        // search for the end of the bin
        for (int j = lastJ; j < amountPrim; j++) {
            prim = primitives[j];
            // if the centroid is over the border, save the last primitive that was not over the border
            if (computePrimitiveCentroid(prim)[axis] > distOnLine) {
                // save the index of the next primitive to test
                lastJ = j-1;
                if (lastJ == -1) {
                    lastJ = 0;
                    break;
                }
                
                // get the areas of the two bins that you made
                aabbLower = computeSpanAABB(primitives.subspan(0, j));
                aabbUpper = computeSpanAABB(primitives.subspan(j, amountPrim - j));
                newArea = areaAabb(aabbLower) + areaAabb(aabbUpper);

                // see if it is better than the last bins
                if (newArea < minArea) {
                    minArea = newArea;
                    MinIndex = lastJ+1;
                }
                break;
            }
        }
    }
    //std::cout << MinIndex << "\n";
    return MinIndex;
}