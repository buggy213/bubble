#include "util/image.h"
#include "util/lodepng.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

float float_from_bytes(char bytes[4]) {
    float output;

    *((char*)(&output) + 3) = bytes[3];
    *((char*)(&output) + 2) = bytes[2];
    *((char*)(&output) + 1) = bytes[1];
    *((char*)(&output) + 0) = bytes[0];

    return output;
}

CGL::HDRImageBuffer load_hdr_pfm(std::ifstream &hdr_file) {
    CGL::HDRImageBuffer buffer;
    std::string s;
    std::getline(hdr_file, s);
    if (s != "PF") {
        std::cerr << "unexpected header: " << s << std::endl;
        std::abort();
    }

    std::getline(hdr_file, s);
    int split = s.find(' ');
    std::string width_str = s.substr(0, split);
    std::string height_str = s.substr(split+1);
    int width = std::stoi(width_str);
    int height = std::stoi(height_str);
    buffer.resize(width, height);

    std::getline(hdr_file, s);
    if (s != "-1.0") {
        std::cerr << "expect little-endian" << std::endl;
        std::abort();
    }

    for (int i = 0; i < height; i += 1) {
        for (int j = 0; j < width; j += 1) {
            int index = i * width + j;
            
            char bytes[4];
            hdr_file.read(bytes, 4);
            buffer.data[index].r = float_from_bytes(bytes);
            hdr_file.read(bytes, 4);
            buffer.data[index].g = float_from_bytes(bytes);
            hdr_file.read(bytes, 4);
            buffer.data[index].b = float_from_bytes(bytes);
        }
    }

    return buffer;    
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "no filename / output filename provided" << std::endl;
        return 1;
    }

    std::string filename {argv[1]};
    std::ifstream pfm_file {filename};

    if (!pfm_file.is_open()) {
        std::cerr << "unable to open file " << filename << std::endl;
        return 1;
    }

    CGL::HDRImageBuffer image = load_hdr_pfm(pfm_file);
    CGL::ImageBuffer sdr_image;
    sdr_image.resize(image.w, image.h);
    image.toColor(sdr_image, 0, 0, image.w, image.h);
    
    uint32_t *frame = &sdr_image.data[0];
    size_t w = sdr_image.w;
    size_t h = sdr_image.h;
    uint32_t *frame_out = new uint32_t[w * h];
    for (size_t i = 0; i < h; ++i) {
        memcpy(frame_out + i * w, frame + (h - i - 1) * w, 4 * w);
    }

    for (size_t i = 0; i < w * h; ++i) {
        frame_out[i] |= 0xFF000000;
    }

    std::string output_filename {argv[2]};
    
    lodepng::encode(output_filename, (unsigned char *)frame_out, w, h);

    return 0;
}