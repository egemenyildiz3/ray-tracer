// Put your includes here
#include "bvh.h"
#include "render.h"
#include "sampler.h"
#include "scene.h"
#include "shading.h"
#include <limits>
#include "interpolate.h"
#include "extra.h"

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <catch2/catch_all.hpp>
#include <glm/glm.hpp>
#include <iostream>
DISABLE_WARNINGS_POP()

// In this file you can add your own unit tests using the Catch2 library.
// You can find the documentation of Catch2 at the following link:
// https://github.com/catchorg/Catch2/blob/devel/docs/assertions.md
//
// These tests are only to help you verify that your code is correct.
// You don't have to hand them in; we will not consider them when grading.
//

void printVector(glm::vec3 input, std::string pre = ""){
    std::cout << pre << "[ " << input[0] << ", " << input[1] << ", " << input[2] << " ]\n";
}

bool floatEqual(float a, float b, float epsilon = 0.001)
{
    return (a - b < epsilon);
};

// Add your tests here, if you want :D
TEST_CASE("StudentTest")
{
    // Add your own tests here...
    SECTION("Bayocentric coordinates")
    {
        glm::vec3 a { .5, .5, 0 };
        glm::vec3 b { 1, 1, 0 };
        glm::vec3 c { -2, 1, 0 };

        glm::vec3 p { 0.2123, 0.75, 0 };

        glm::vec3 sol = computeBarycentricCoord(a, b, c, p);

        printVector(sol);

        glm::vec3 check = sol.x * a + sol.y * b + sol.z * c;

        printVector(check, "Actual");
        printVector(p, "Expected");

        CHECK((floatEqual(check.x, p.x) && floatEqual(check.y, p.y) && floatEqual(check.z, p.z)));
    }
    
    SECTION("computePrimitiveAABB")
    {
        BVHInterface::Primitive triangle;
        triangle.v0 = { { 0.0f, 0.0f, 0.0f } };
        triangle.v1 = { { 1.0f, 0.0f, 0.0f } };
        triangle.v2 = { { 0.0f, 1.0f, 0.0f } };
        AxisAlignedBox v = { .lower = { 0, 0, 0 }, .upper = { 1, 1, 0 } };
        AxisAlignedBox box = computePrimitiveAABB(triangle);
        CHECK(box.lower ==  v.lower);
        CHECK(box.upper == v.upper);
    }

    SECTION("computeSpanAABB") 
    {
        BVHInterface::Primitive triangle0 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 1.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 1.0f } }
        };
        BVHInterface::Primitive triangle1 = {
            .v0 = { { 7.0f, 1.0f, 0.0f } },
            .v1 = { { 1.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle2 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 1.0f, 4.0f, 0.0f } },
            .v2 = { { 4.0f, 1.0f, 0.0f } }
        };
        std::vector triangles = { triangle0, triangle1, triangle2 };
        AxisAlignedBox v = { .lower = { 0, 0, 0 }, .upper = { 7, 4, 1 } };
        AxisAlignedBox box = computeSpanAABB(triangles);

        printVector(box.upper);
        printVector(v.upper);
        CHECK(box.lower == v.lower);
        CHECK(box.upper == v.upper);
    }

    SECTION("computePrimitiveCentroid")
    {
        BVHInterface::Primitive triangle = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 1.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 0.0f } }
        };
        glm::vec3 v = { 1.0f/3.0f,1.0f/3.0f,0 };
        glm::vec3 center = computePrimitiveCentroid(triangle);
        CHECK(v == center);
    }

    SECTION("computeAABBLongestAxis") 
    {
        AxisAlignedBox aabb { .lower = { 0,1,-2 }, .upper = {10,10,10} };
        uint32_t v = 2;
        uint32_t longest = computeAABBLongestAxis(aabb);
        CHECK(v == longest);
    }

    SECTION("splitPrimitivesByMedian")
    {
        AxisAlignedBox box {
            .lower = { 0, 0, 0 },
            .upper = { 10, 10, 10 }
        };
        int axis = 0; 
        BVHInterface::Primitive triangle0 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 0.5f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle1 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 5.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle2 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 1.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 0.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle3 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 2.0f, 0.0f, 5.0f } },
            .v2 = { { 0.0f, 1.0f, 10.0f } }
        };
        BVHInterface::Primitive triangle4 = {
            .v0 = { { 0.0f, 7.0f, 0.0f } },
            .v1 = { { 1.0f, 0.0f, 0.0f } },
            .v2 = { { 8.0f, 1.0f, 0.0f } }
        };
        std::vector<BVHInterface::Primitive> triangles { triangle0, triangle1, triangle2, triangle3, triangle4 };
        
        size_t split = splitPrimitivesByMedian(box, axis, triangles);
        CHECK(split == 3);
        bool firstHalve = (triangles[0].operator==(triangle0) || triangles[0].operator==(triangle2) || triangles[0].operator==(triangle3))
                && triangles[1].operator==(triangle0) || triangles[1].operator==(triangle2) || triangles[1].operator==(triangle3) 
                && triangles[2].operator==(triangle0) || triangles[2].operator==(triangle2) || triangles[2].operator==(triangle3);
        CHECK(firstHalve == true);
        bool secondHalve = triangles[3].operator==(triangle1) || triangles[3].operator==(triangle4) 
                && triangles[4].operator==(triangle1) || triangles[4].operator==(triangle4);
        CHECK(secondHalve == true);
    }

    SECTION("splitPrimitivesBySAHBin_sorted")
    {
        AxisAlignedBox box {
            .lower = { 0, 0, 0 },
            .upper = { 10, 10, 10 }
        };
        int axis = 0;
        BVHInterface::Primitive triangle0 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 0.5f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle1 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 1.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 0.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle2 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 2.0f, 0.0f, 5.0f } },
            .v2 = { { 0.0f, 1.0f, 10.0f } }
        };
        BVHInterface::Primitive triangle3 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 5.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle4 = {
            .v0 = { { 10.0f, 7.0f, 0.0f } },
            .v1 = { { 7.0f, 0.0f, 0.0f } },
            .v2 = { { 8.0f, 1.0f, 0.0f } }
        };
        std::vector<BVHInterface::Primitive> triangles { triangle0, triangle1, triangle2, triangle3, triangle4 };

        size_t split = splitPrimitivesBySAHBin(box, axis, triangles);
        CHECK(split == 4);
    }

        SECTION("splitPrimitivesBySAHBin_unsorted")
    {
        AxisAlignedBox box {
            .lower = { 0, 0, 0 },
            .upper = { 10, 10, 10 }
        };
        int axis = 0;
        BVHInterface::Primitive triangle0 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 0.5f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle1 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 1.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 0.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle2 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 2.0f, 0.0f, 5.0f } },
            .v2 = { { 0.0f, 1.0f, 10.0f } }
        };
        BVHInterface::Primitive triangle3 = {
            .v0 = { { 0.0f, 0.0f, 0.0f } },
            .v1 = { { 5.0f, 0.0f, 0.0f } },
            .v2 = { { 0.0f, 1.0f, 0.0f } }
        };
        BVHInterface::Primitive triangle4 = {
            .v0 = { { 10.0f, 7.0f, 0.0f } },
            .v1 = { { 7.0f, 0.0f, 0.0f } },
            .v2 = { { 8.0f, 1.0f, 0.0f } }
        };
        //unsorted
        std::vector<BVHInterface::Primitive> triangles { triangle3, triangle1, triangle4, triangle0, triangle2 };

        size_t split = splitPrimitivesBySAHBin(box, axis, triangles);
        CHECK(split == 4);
    }
}

// The below tests are not "good" unit tests. They don't actually test correctness.
// They simply exist for demonstrative purposes. As they interact with the interfaces
// (scene, bvh_interface, etc), they allow you to verify that you haven't broken
// our grading interface. They should compile without changes. If they do
// not compile, neither will our grading tests!
TEST_CASE("InterfaceTest")
{
    // Setup a RenderState object with some defaults
    Features features = {
        .enableShading = true,
        .enableAccelStructure = false, // BVH is not actually active r.n.
        .shadingModel = ShadingModel::Lambertian
    };
    Scene scene = loadScenePrebuilt(SceneType::CornellBox, DATA_DIR);
    BVH bvh(scene, features);
    RenderState state = { .scene = scene, .features = features, .bvh = bvh, .sampler = {} };

    SECTION("BVH generation")
    {
        // There's something in here?
        CHECK(!state.bvh.primitives().empty());
    }

    SECTION("BVH traversal")
    {
        Ray ray = { .origin = glm::vec3(0), .direction = glm::vec3(1) };
        HitInfo hitInfo;

        // Hit something?
        CHECK(state.bvh.intersect(state, ray, hitInfo));
        CHECK(ray.t != std::numeric_limits<float>::max());
    }

    SECTION("Hit shading")
    {
        Ray ray = { .origin = glm::vec3(0), .direction = glm::vec3(1) };
        HitInfo hitInfo;
        state.bvh.intersect(state, ray, hitInfo);

        // Shaded something?
        glm::vec3 Lo = computeShading(state, ray.direction, -ray.direction, glm::vec3(1), hitInfo);
        CHECK(glm::any(glm::notEqual(Lo, glm::vec3(0))));
    }
}
