#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <algorithm> 
#include "Cell.h"
#include "common_macro.hpp"

#define TOT_GEOM_FEATURES 29

class Feature {
  public:
    float linearity;
    float planarity;
    float sphericity;
    float omnivariance;
    float anisotropy;
    float eigenentropy;
    float sum_of_eigenvalues;
    float curvature;
    float angle;
    float goodness_of_fit;
    float roughness;
    float nvx;
    float nvy;
    float nvz;
    float unevenness;
    float surface_density;
    float z_diff;
    float intensity_ave;
    float intensity_variance;
    float intensity_diff;
    float max_intensity;
    float min_intensity;
    float z_ave;
    float z_variance;
    float z_std;
    float z_min;
    float z_max;
    float radius1_intensity_distribution;
    float radius2_intensity_distribution;


    double cx, cy, cz, d1, d2, d3, d, normal_magnitude, numpoints;
    double numpoints_inverse;
    double a11, a12, a13, a22, a23, a33, ax, ay, az;

    Eigen::MatrixXd matA1;
    Eigen::Vector3d *p;

    std::vector<float> derived_features;

    void computeCorrelationMatrix(std::vector<int> &points_idx, std::vector<Eigen::Vector3d> &points);

    // Feature();

    std::string toString();

    std::vector<float> toVector();

    std::vector<float> toVectorTransformed();
    void toVectorTransformed(std::vector<float> &feature);

    void toFile(std::ofstream &out);

    int fromFileLine(std::ifstream &in, int derived_features_num);
    int ignoreFeatureFromFile(std::ifstream &in, int derived_features_num);
    int computeFeatures(Cell *cell, Eigen::MatrixXd &scene_normal, std::vector<Eigen::Vector3d> &points, std::vector<float> &intensities);

};
