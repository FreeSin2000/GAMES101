//
// Created by LEI XU on 4/27/19.
//

#ifndef RASTERIZER_TEXTURE_H
#define RASTERIZER_TEXTURE_H
#include "global.hpp"
#include <eigen3/Eigen/Eigen>
#include <opencv2/opencv.hpp>
class Texture{
private:
    cv::Mat image_data;

public:
    Texture(const std::string& name)
    {
        image_data = cv::imread(name);
        cv::cvtColor(image_data, image_data, cv::COLOR_RGB2BGR);
        width = image_data.cols;
        height = image_data.rows;
    }

    int width, height;

    Eigen::Vector3f getColor(float u, float v)
    {
        auto u_img = u * width;
        auto v_img = (1 - v) * height;
        auto color = image_data.at<cv::Vec3b>(v_img, u_img);
        return Eigen::Vector3f(color[0], color[1], color[2]);
    }

    Eigen::Vector3f getColorBilinear(float u, float v)
    {
        auto u_img = u * width;
        auto v_img = (1 - v) * height;
        
        int u0 = std::max(0, (int)std::floor(u_img - 0.5));
        int u1 = std::min(width - 1, (int)std::ceil(u_img));

        int v0 = std::max(0, (int)std::floor(v_img - 0.5));
        int v1 = std::min(height - 1, (int)std::ceil(v_img));

        auto color = image_data.at<cv::Vec3b>(v_img, u_img);
        auto to_eigen = [](cv::Vec3b c) {
            return Eigen::Vector3f(c[0], c[1], c[2]);
        };
        auto c00 = to_eigen(image_data.at<cv::Vec3b>(v0, u0));
        auto c01 = to_eigen(image_data.at<cv::Vec3b>(v0, u1));
        auto c10 = to_eigen(image_data.at<cv::Vec3b>(v1, u0));
        auto c11 = to_eigen(image_data.at<cv::Vec3b>(v1, u1));


        float alpha = v_img - v0;
        float beta = u_img - u0;
        auto c0 = (c00 * alpha + (1.0 - alpha) * c10);
        auto c1 = (c01 * alpha + (1.0 - alpha) * c11);

        return  c0 * beta + (1.0 - beta) * c1;
    }

};
#endif //RASTERIZER_TEXTURE_H
