#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()


// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.
    position = glm::vec3(0.0);
    color = glm::vec3(0.0);
    position = light.endpoint0 + sample * (light.endpoint1 - light.endpoint0);
    color = light.color0 + sample * (light.color1 - light.color0);
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    position = light.v0 + sample.x * light.edge01 + sample.y * light.edge02;
    const float a = 1.0f - sample.x;
    const float b = 1.0f - sample.y;
    color = a * b * light.color0 + a * sample.y * light.color2 + sample.x * sample.y * light.color3 + sample.x * b * light.color1;
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        return true;
    } else {
        Ray shadowRay;
        shadowRay.origin = ray.origin + ray.t * ray.direction + glm::normalize(lightPosition - (ray.origin + ray.t * ray.direction)) * 0.0001f;
        shadowRay.direction = glm::normalize(lightPosition - (ray.origin + ray.t * ray.direction));
        HitInfo sh;
        state.bvh.intersect(state, shadowRay, sh);
        float a = glm::distance(shadowRay.origin, lightPosition);
        if (shadowRay.t < a) {
            return false;
        } else {    
            return true;
        }
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo)) {
            return lightColor;
    }
    bool check = false;
    glm::vec3 newCol = lightColor;
    glm::vec3 pos = (ray.origin + ray.t * ray.direction) - lightPosition;
    glm::vec3 vec1 = { 1.0f, 1.0f, 1.0f };

    Ray ray1;
    HitInfo newHit = hitInfo;
    Ray newRay = ray;
    glm::vec3 newDir = ray.direction;
       do {
            glm::vec3 shDir = glm::normalize(lightPosition - (newRay.origin + newRay.t * newRay.direction));

            ray1.origin = newRay.origin + newRay.t * newRay.direction + shDir * 0.001f; 
            ray1.direction = shDir;
            pos = lightPosition - ray1.origin;
            
            ray1.t = FLT_MAX;
            check = state.bvh.intersect(state, ray1, newHit);
            glm::vec3 updateK = sampleMaterialKd(state, newHit);
            
            if (glm::length(pos) <= ray1.t + 0.001f)
            break;

            if (!check)
                {}
            else {
            vec1 = vec1 * ( updateK * (1.0f - newHit.material.transparency));
            }
            newHit = hitInfo;
            newRay = ray1;

       }while (newRay.t + 0.001f <= glm::length(pos));

    glm::vec3 ret = vec1 * lightColor;
    return ret;
}

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// +6
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    if (visibilityOfLightSample(state, light.position, light.color, ray, hitInfo) != glm::vec3 { 0, 0, 0 }) {
        glm::vec3 p = ray.origin + ray.t * ray.direction;
        glm::vec3 l = glm::normalize(light.position - p);
        glm::vec3 v = -ray.direction;
        return visibilityOfLightSample(state, light.position, light.color, ray, hitInfo) * computeShading(state, v, l, light.color, hitInfo);
    }
    return glm::vec3 { 0, 0, 0 };
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    glm::vec3 total(0.0f);

    for (uint32_t i = 0; i < numSamples; ++i) {
        float sample = state.sampler.next_1d();
        glm::vec3 posLight;
        glm::vec3 colLight;
        sampleSegmentLight(sample, light, posLight, colLight);
        PointLight newLight;
        newLight.position = posLight;
        newLight.color = colLight;
        total += computeContributionPointLight(state, newLight, ray, hitInfo);
    }

    return total;
}

// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    if (numSamples == 0) {
        return glm::vec3 { 0.0f };
    }

    glm::vec3 point = ray.origin + ray.t * ray.direction;
    glm::vec3 dvec = -ray.direction;

    glm::vec3 add { 0.0f };
    glm::vec3 m { 0.0f };
    glm::vec3 n { 0.0f };

    for (int i = 0; i < numSamples; i++) {
        glm::vec2 curpos = state.sampler.next_2d();
        sampleParallelogramLight(curpos, light, m, n);
        glm::vec3 lightColourAtTarget = visibilityOfLightSample(state, m, n, ray, hitInfo);
        if (lightColourAtTarget != glm::vec3 { 0.0f }) {
            add += computeShading(state, dvec, glm::normalize(m - point), lightColourAtTarget, hitInfo);
        }
    }
    float numSaFloat = static_cast<float>(numSamples);
    glm::vec3 retu = add / numSaFloat;
    return retu;
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}
