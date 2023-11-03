#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include <draw.h>
#include <iostream>
#include <texture.cpp>
#include <iostream>


glm::vec3 avg(std::vector<glm::vec3> input)
{
    glm::vec3 output {};
    for (glm::vec3 i : input) {
        output += (1.0f / input.size()) * i;
    }
    return output;
}

void printVector(glm::vec3 input, std::string pre = "")
{
    std::cout << pre << "[ " << input[0] << ", " << input[1] << ", " << input[2] << " ]\n";
}

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

    float offset = features.extra.offsetFocal;
    int amountSamples = 20;
    float lensRadius = 1.0f;

#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            // Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
            // Note; we seed the sampler for consistenct behavior across frames
            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };

            std::vector<glm::vec3> color;
            // do amountSamples random offsets from the base of the ray, move them on the edge of a circle
            for (int i = 0; i < amountSamples; i++) {
                glm::vec3 rand = { state.sampler.next_2d() * lensRadius * 2.0f - 1.0f, 0 }; // make a random variable that goes from -lenssize to + lenssize
                auto rays = generatePixelRays(state, camera, { x , y }, screen.resolution()); // this generates rays from the camera to the xy?
                for (int j = 0; j < rays.size(); j++) { // it could be that it made more rays, jittering samples
                    auto i = rays[j];
                    i.origin.x -= rand[0]; // ofset the origin on a circle and scale
                    i.origin.y -= rand[1];
                    i.direction.operator+=(rand/offset);
                    rays[j] = i;
                }
                color.push_back(renderRays(state, rays)); // add the new colors to a array
            }
            
            auto L = avg(color); // avg that array
            screen.setPixel(x, y, L);
        }
    }
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

//a simple factorial function
int fac(int n)
{
    if (n == 1 || n == 0) {
        return 1;
    }
    return n * fac(n - 1);
}

void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    Screen bloomy(image.resolution(), true);

    //iterate over pixels o image and check if RGB has a value greater than or equal to 0.8.
    //If that is the case, bloomy is set to the same color. 
    //This part is to select the "bright" parts.
    for (int a = 0; a < bloomy.resolution().x; a++) {
        for (int b = 0; b < bloomy.resolution().y; b++) {
            glm::vec3 pix = image.pixels()[image.indexAt(a, b)];
            if (pix.r >= 0.8f || pix.g >= 0.8f || pix.b >= 0.8f) {
                bloomy.setPixel(a, b, pix);
            }
        }
    }

    // calculate the "checks" for a Gaussian blur filter for blurring the bright areas.
    int limit = 10;
    std::vector<float> checks;
    for (int x = 0; x < limit; x++) {
        float combination = (float)fac(limit - 1) / ((float)fac(limit - x - 1) * (float)fac(x));
        checks.push_back(combination);
    }

    //a quick normalization
    float add = 0.0f;
    for (const auto& x : checks) {
        add += x;
    }
    for (auto& x : checks) {
        x /= add;
    }
    
    //horizontal blurring
    Screen newPhoto(bloomy.resolution(), true);

    for (int m = 0; m < bloomy.resolution().y; m++) {
        for (int n = 0; n < bloomy.resolution().x; n++) {
            float r = 0.0f, g = 0.0f, b = 0.0f;

            //to apply a weighted average
            for (int w = 0; w < limit; w++) {
                //calculating the index
                int q = m - (limit - 1) / 2 + w;
                // to ensure it is not out of bounds
                if (q < bloomy.resolution().x && q >= 0) {
                    // updating the RGB values from our calculations and the bloomy
                    glm::vec3 pixi = bloomy.pixels()[bloomy.indexAt(q, n)];
                    r = r + pixi.r * checks[w];
                    g = g + pixi.g * checks[w];
                    b = b + pixi.b * checks[w];
                }
            }
            newPhoto.setPixel(m, n, glm::vec3(r, g, b));
        }
    }

    Screen blurry(bloomy.resolution(), true);

    //vertical blurring 
    for (int m = 0; m < bloomy.resolution().y; m++) {
        for (int n = 0; n < bloomy.resolution().x; n++) {
            float r = 0.0f, g = 0.0f, b = 0.0f;

            //to apply a weighted average
            for (int w = 0; w < limit; w++) {
                // calculating the index
                int q = n - (limit - 1) / 2 + w;
                //to ensure it is not out of bounds
                if (q >= 0 && q < bloomy.resolution().y) {
                    //updating the RGB values from our calculations and the newPhoto
                    glm::vec3 pixi = newPhoto.pixels()[newPhoto.indexAt(m, q)];
                    r = r + pixi.r * checks[w];
                    g = g + pixi.g * checks[w];
                    b = b + pixi.b * checks[w];
                }
            }
            blurry.setPixel(m, n, glm::vec3(r, g, b));
        }
    }

    //mix the original with the final result in order to obtain the bloom effect
    for (int a = 0; a < image.resolution().y; a++) {
        for (int b = 0; b < image.resolution().x; b++) {
            glm::vec3 pixi = blurry.pixels()[blurry.indexAt(a, b)] + image.pixels()[image.indexAt(a, b)];
            //ensure no RBG values over 1
            for (int m = 0; m < 3; m++) {
                if (pixi[m] > 1) {
                    pixi[m] = 1;
                }
            }
            image.setPixel(a, b, pixi);
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
    // ...

    float numSamples = state.features.extra.numGlossySamples;
    float n = hitInfo.material.shininess / 64;

    Ray r = generateReflectionRay(ray, hitInfo);
    drawRay(r, { 0, 1, 0 });

    glm::vec3 w = glm::normalize(r.direction);
    glm::vec3 t = { 1, 0, 0 };
    if (r.direction.y == r.direction.z == 0)
        t = { 0, 1, 0 };
    glm::vec3 v = glm::normalize(glm::cross(r.direction, t));
    glm::vec3 u = glm::cross(w, v);

    // We use the method from the book in 14.4.1 to get a random vector direction uniformly distributed over a hemisphere.
    for (int i = 0; i < numSamples; i++) {
        glm::vec2 xi = state.sampler.next_2d();
        float radius = glm::sqrt(xi[0]) * n;
        float phi = glm::two_pi<float>() * xi[1];

        float a = glm::sin(phi) * radius; 
        float b = glm::cos(phi) * radius; 

        r.direction = w + a * u + b * v;

        if (glm::dot(hitInfo.normal, r.direction) < 0)
            continue;

        hitColor += hitInfo.material.ks * renderRay(state, r, rayDepth + 1) / numSamples;
    }
}

int sign(float v) {
    if (v < 0)
        return -1;
    if (v == 0)
        return 0;
    return 1;
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, const Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        float x = ray.direction.x;
        float y = ray.direction.y;
        float z = ray.direction.z;
        float absX = glm::abs(x);
        float absY = glm::abs(y);
        float absZ = glm::abs(z);

        glm::vec2 coords;

        int i = 0;
        if (absX >= glm::max(absY, absZ)) {
            coords = glm::vec2 { -sign(x) * z, -y } / absX;
            if (x < 0)
                i = 0;
            else
                i = 3;
        }
        if (absY >= glm::max(absX, absZ)) {
            coords = glm::vec2 { x, sign(y) * z } / absY;
            if (y < 0)
                i = 1;
            else
                i = 4;
        }
        if (absZ >= glm::max(absX, absY)) {
            coords = glm::vec2 { sign(z) * x, -y } / absZ;
            if (z < 0)
                i = 2;
            else
                i = 5;
        }

        const std::vector<Image>& sides = state.scene.envMap;
        coords = (1.0f + coords) / 2.0f;
        coords[1] = 1 - coords[1];
        if (sides.size() == 6 && coords.x > 0 && coords.y > 0 && coords.x < 1 && coords.y < 1)
            return sampleTextureNearest(sides[i], coords);
        return { 0, 0, 0 };
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
