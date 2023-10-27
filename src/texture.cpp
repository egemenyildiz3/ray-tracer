#include "texture.h"
#include "render.h"
#include <framework/image.h>


glm::vec3 getPixel(Image image, int& i, int& j) {
    return image.pixels[i + j * image.width];
}

glm::vec3 getPixel(Image image, glm::vec2& p)
{
    return image.pixels[p.x + p.y * image.width];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    
    assert(texCoord.x >= 0 && texCoord.x <= 1);
    assert(texCoord.y >= 0 && texCoord.y <= 1);
    assert(image.pixels.size() == image.width * image.height);

    int i = glm::round(texCoord.x * image.width  - 0.5);
    int j = glm::round(texCoord.y * image.height - 0.5);

    return getPixel(image, i, j);
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    // Brightspace video: https://brightspace.tudelft.nl/d2l/le/content/595314/viewContent/3512131/View

    assert(texCoord.x >= 0 && texCoord.x <= 1);
    assert(texCoord.y >= 0 && texCoord.y <= 1);
    assert(image.pixels.size() == image.width * image.height);

    int x1 = glm::floor(texCoord.x * image.width );
    int x2 = glm::ceil( texCoord.x * image.width );
    int y1 = glm::floor(texCoord.y * image.height);
    int y2 = glm::ceil( texCoord.y * image.height);

    glm::vec2 p1 = { x1, y1 };
    glm::vec2 p2 = { x2, y1 };
    glm::vec2 p3 = { x1, y2 };
    glm::vec2 p4 = { x2, y2 };

    float a = texCoord.x * image.width  - x1; 
    float b = texCoord.y * image.height - y1;

    glm::vec3 c1 = getPixel(image, p1) * (1-a) * (1-b);
    glm::vec3 c2 = getPixel(image, p2) * ( a ) * (1-b);
    glm::vec3 c3 = getPixel(image, p3) * (1-a) * ( b );
    glm::vec3 c4 = getPixel(image, p4) * ( a ) * ( b );

    return c1 + c2 + c3 + c4;
}