#include "texture.h"
#include "render.h"
#include <framework/image.h>


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

    int i = glm::round(texCoord.x * image.width - 0.5);
    int j = glm::round((1 - texCoord.y) * image.height - 0.5);

    return image.pixels[i + j * image.width];
}

glm::vec3 getPixel(const Image& image, const glm::vec2& p)
{
    return image.pixels[p.x + p.y * image.width];
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

    float texX = texCoord.x; 
    float texY = texCoord.y;
    if (texX == 0)
        texX++;
    if (texY == 0)
        texY++;
    if (texX == image.width)
        texX--;
    if (texY == image.height)
        texY--;
        
    const glm::vec2 imgCoords = { texX * image.width - 0.5, (1 - texY) * image.height - 0.5 };

    int x1 = glm::floor(imgCoords.x);
    int y1 = glm::floor(imgCoords.y);

    const glm::vec2 p1 = { x1,     y1 };
    const glm::vec2 p2 = { x1 + 1, y1 };
    const glm::vec2 p3 = { x1,     y1 + 1 };
    const glm::vec2 p4 = { x1 + 1, y1 + 1 };

    float a = imgCoords.x - x1; 
    float b = imgCoords.y - y1;

    glm::vec3 res;
    res += getPixel(image, p1) * (1-a) * (1-b);
    res += getPixel(image, p2) * ( a ) * (1-b);
    res += getPixel(image, p3) * (1-a) * ( b );
    res += getPixel(image, p4) * ( a ) * ( b );

    return res;
}