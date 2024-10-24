#include "Feature.h"

void Feature::computeCorrelationMatrix(std::vector<int> &points_idx, std::vector<Eigen::Vector3d> &points) {
    
    cx=0; cy=0; cz=0;  
    a11 = 0, a12 = 0, a13 = 0, a22 = 0, a23 = 0, a33 = 0;
    
    // compute eigenvalues and eigenvectors
    matA1 = Eigen::MatrixXd::Zero(3, 3);

    for (auto p_idx: points_idx) {
        p = &(points[p_idx]);
        cx += (*p)(0);
        cy += (*p)(1);
        cz += (*p)(2);
    }

    cx *= numpoints_inverse; cy *= numpoints_inverse;  cz *= numpoints_inverse;

    for (auto p_idx: points_idx) {
        p = &(points[p_idx]);
        ax = (*p)(0) - cx;
        ay = (*p)(1) - cy;
        az = (*p)(2) - cz;

        a11 += ax * ax; a12 += ax * ay; a13 += ax * az;
        a22 += ay * ay; a23 += ay * az;
        a33 += az * az;
    }
    a11 *= numpoints_inverse; a12 *= numpoints_inverse; a13 *= numpoints_inverse; 
    a22 *= numpoints_inverse; a23 *= numpoints_inverse; a33 *= numpoints_inverse;

    matA1 << a11, a12, a13, a12, a22, a23, a13, a23, a33;
}

std::string Feature::toString() {
  std::stringstream ss("");

  ss << std::to_string(linearity) << " ";
  ss << std::to_string(planarity) << " ";
  ss << std::to_string(sphericity) << " ";
  ss << std::to_string(omnivariance) << " ";
  ss << std::to_string(anisotropy) << " ";
  ss << std::to_string(eigenentropy) << " ";
  ss << std::to_string(sum_of_eigenvalues) << " ";
  ss << std::to_string(curvature) << " ";
  ss << std::to_string(angle) << " ";
  ss << std::to_string(goodness_of_fit) << " ";
  ss << std::to_string(roughness) << " ";
  ss << std::to_string(nvx) << " ";
  ss << std::to_string(nvy) << " ";
  ss << std::to_string(nvz) << " ";
  ss << std::to_string(unevenness) << " ";
  ss << std::to_string(surface_density) << " ";
  ss << std::to_string(z_diff) << " ";
  ss << std::to_string(intensity_ave) << " ";
  ss << std::to_string(intensity_variance) << " ";
  ss << std::to_string(intensity_diff) << " ";
  ss << std::to_string(max_intensity) << " ";
  ss << std::to_string(min_intensity) << " ";
  ss << std::to_string(z_ave) << " ";
  ss << std::to_string(z_variance) << " ";
  ss << std::to_string(z_max) << " ";
  ss << std::to_string(z_min) << " ";
  ss << std::to_string(radius1_intensity_distribution) << " ";
  ss << std::to_string(radius2_intensity_distribution) << " ";


  ss << " derived: ";
  for (int i=0; i<(int)derived_features.size(); i++)
    ss << " " << std::to_string(derived_features[i]);

  return ss.str();
}

std::vector<float> Feature::toVector() {
  std::vector<float> feature = std::vector<float>(TOT_GEOM_FEATURES + (int)derived_features.size());
  feature[0] = linearity;
  feature[1] = planarity;
  feature[2] = sphericity;
  feature[3] = omnivariance;
  feature[4] = anisotropy;
  feature[5] = eigenentropy;
  feature[6] = sum_of_eigenvalues;
  feature[7] = curvature;
  feature[8] = angle;
  feature[9] = goodness_of_fit;
  feature[10] = roughness;
  feature[11] = nvx;
  feature[12] = nvy;
  feature[13] = nvz;
  feature[14] = unevenness;
  feature[15] = surface_density;
  feature[16] = z_diff;
  feature[17] = intensity_ave;
  feature[18] = intensity_variance;
  feature[19] = intensity_diff;
  feature[20] = max_intensity;
  feature[21] = min_intensity;
  feature[22] = z_ave;
  feature[23] = z_variance;
  feature[24] = z_std;
  feature[25] = z_max;
  feature[26] = z_min;
  feature[27] = radius1_intensity_distribution;
  feature[28] = radius2_intensity_distribution;



  for (int i=0; i<(int)derived_features.size(); i++)
    feature[TOT_GEOM_FEATURES+i] = derived_features[i];
  return feature;
}

std::vector<float> Feature::toVectorTransformed() {
  std::vector<float> feature = toVector();
  for (int i=0; i<TOT_GEOM_FEATURES; i++) // TOT_GEOM_FEATURES because derived feature are already transformed
    feature[i] = (double)log( fabs((double)feature[i]) + 1e-4 );
  return feature;
}

void Feature::toVectorTransformed(std::vector<float> &feature) {
  // feature[0] = linearity;
  // feature[1] = planarity;
  // feature[2] = sphericity;
  // feature[3] = omnivariance;
  // feature[4] = anisotropy;
  // feature[5] = eigenentropy;
  // feature[6] = sum_of_eigenvalues;
  // feature[7] = curvature;
  // feature[8] = angle;
  // feature[9] = goodness_of_fit;
  // feature[10] = roughness;
  // feature[11] = nvx;
  // feature[12] = nvy;
  // feature[13] = nvz;
  // feature[14] = unevenness;
  // feature[15] = surface_density;
  // feature[16] = z_diff;
  // for (int i=0; i<(int)derived_features.size(); i++)
  //   feature[TOT_GEOM_FEATURES+i] = derived_features[i];
  // for (int i=0; i<TOT_GEOM_FEATURES; i++) // TOT_GEOM_FEATURES because derived feature are already transformed
  //   feature[i] = log( abs(feature[i]) + 1e-4 );
}

void Feature::toFile(std::ofstream &out) {
  out.write( reinterpret_cast<const char*>( &(linearity) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(planarity) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(sphericity) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(omnivariance) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(anisotropy) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(eigenentropy) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(sum_of_eigenvalues) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(curvature) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(angle) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(goodness_of_fit) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(roughness) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(nvx) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(nvy) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(nvz) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(unevenness) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(surface_density) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(z_diff) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(intensity_ave) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(intensity_variance) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(intensity_diff) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(max_intensity) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(min_intensity) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(z_ave) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(z_variance) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(z_std) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(z_max) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(z_min) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(radius1_intensity_distribution) ), sizeof( float ));
  out.write( reinterpret_cast<const char*>( &(radius2_intensity_distribution) ), sizeof( float ));



  for (int i=0; i<(int)derived_features.size(); i++)
    out.write( reinterpret_cast<const char*>( &(derived_features[i]) ), sizeof( float ));

}

int Feature::fromFileLine(std::ifstream &in, int derived_features_num) {
  if (!in.read( reinterpret_cast< char*>( &(linearity) ), sizeof( float ))) return 0;
    in.read( reinterpret_cast< char*>( &(planarity) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(sphericity) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(omnivariance) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(anisotropy) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(eigenentropy) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(sum_of_eigenvalues) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(curvature) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(angle) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(goodness_of_fit) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(roughness) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(nvx) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(nvy) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(nvz) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(unevenness) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(surface_density) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(z_diff) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(intensity_ave) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(intensity_variance) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(intensity_diff) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(max_intensity) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(min_intensity) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(z_ave) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(z_variance) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(z_std) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(z_max) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(z_min) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(radius1_intensity_distribution) ), sizeof( float ));
    in.read( reinterpret_cast< char*>( &(radius2_intensity_distribution) ), sizeof( float ));



    if (derived_features_num<=0) return 1;
    derived_features.resize(derived_features_num); 
    for (int i=0; i<derived_features_num; i++)
      in.read( reinterpret_cast< char*>( &(derived_features[i]) ), sizeof( float ));


  return 1;
}

int Feature::ignoreFeatureFromFile(std::ifstream &in, int derived_features_num) {
  if (!in.read( reinterpret_cast< char*>( &(linearity) ), sizeof( float ))) return 0;
  
  in.ignore(sizeof(float)*(TOT_GEOM_FEATURES - 1));
  if (derived_features_num<=0) return 1;
  
  in.ignore(sizeof(float)*derived_features_num);

  return 1;
}



int Feature::computeFeatures(Cell *cell, Eigen::MatrixXd &scene_normal, std::vector<Eigen::Vector3d> &points, std::vector<float> &intensities) {

  cx = 0; cy = 0; cz = 0;
  
  linearity = 0;
  planarity = 0;
  sphericity = 0;
  omnivariance = 0;
  anisotropy = 0;
  eigenentropy = 0;
  sum_of_eigenvalues = 0;
  curvature = 0;
  angle = 0;
  goodness_of_fit = 0;
  roughness = 0;
  nvx = 0;
  nvy = 0;
  nvz = 0;
  unevenness = 0;
  surface_density = 0;
  z_diff = 0;
  intensity_ave = 0;
  intensity_variance = 0;
  intensity_diff = 0;
  max_intensity = intensities[cell->points_idx[0]];
  min_intensity = intensities[cell->points_idx[0]];
  z_ave = 0;
  z_variance = 0;
  z_std = 0;
  z_max = 0;
  z_min = 0;
  radius1_intensity_distribution = 0;
  radius2_intensity_distribution = 0;



  float sum_intensity = 0.0;


  
  numpoints = cell->points_idx.size();
  numpoints_inverse = 1.0/numpoints;
  
  computeCorrelationMatrix(cell->points_idx, points);


  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(matA1);
  if (eigensolver.info() != Eigen::Success) {
    return 0;
  }
  // store eigenvalues in an easy way
  d1 = eigensolver.eigenvalues()(0);
  d2 = eigensolver.eigenvalues()(1);
  d3 = eigensolver.eigenvalues()(2);

  if (d1<1e-16)  {
    d1=1e-16;
  }

  if (d2<1e-16) d2=1e-16;

  if (d3<1e-16) d3=1e-16;

  nvx = eigensolver.eigenvectors()(0, 0);
  nvy = eigensolver.eigenvectors()(1, 0);
  nvz = eigensolver.eigenvectors()(2, 0);

  double d1_inverse = 1.0 / d1;          

  /// COVARIANCE-BASED
  linearity    = (d1 - d2)  * d1_inverse; 
  planarity    = (d2 - d3)  * d1_inverse;
  sphericity   = d3  * d1_inverse;       
  omnivariance = std::cbrt(d1*d2*d3);
  anisotropy   = (d1 - d3) * d1_inverse;  
  eigenentropy = - ( d1*std::log(d1) + d2*std::log(d2) + d3*std::log(d3) );
  sum_of_eigenvalues = d1 + d2 + d3;
  curvature    = d3 / sum_of_eigenvalues;


  d = - ( nvx * cx + nvy * cy + nvz * cz );
  normal_magnitude = eigensolver.eigenvectors().col(0).norm();
  if (!normal_magnitude) {
    // std::cout << "normal_magnitude\n";
    return 0;
  }

  angle = std::acos(nvz); 
  
  double min = 10000000, max = -10000000;

  double normal_magnitude_inverse = 1.0 / normal_magnitude; 

  Eigen::Vector3d max_intensity_point, median_intensity_point, min_intensity_point;
  
  for (auto p_idx : cell->points_idx) {
      p = &(points[p_idx]);
      goodness_of_fit +=
              std::abs((nvx*(*p)(0) + nvy*(*p)(1) + nvz*(*p)(2)) + d) * normal_magnitude_inverse;
      roughness += NUM2SQ((*p)(2) - cz);

      double z = SCALAR_PRODUCT_2p((*p), scene_normal) * scene_normal(2);

      if (z>max) max = z;
      else if (z<min) min = z;

      z_ave += z;


      float current_intensity = intensities[p_idx];
      intensity_ave += intensities[p_idx];
      if (current_intensity < min_intensity) {
          min_intensity = current_intensity;
          min_intensity_point = points[p_idx];
      }
      if (current_intensity > max_intensity) {
          max_intensity = current_intensity;
          max_intensity_point = points[p_idx];
      }
  }
  intensity_ave *= numpoints_inverse;
  intensity_diff = max_intensity - min_intensity; 
  
  roughness *= numpoints_inverse;
  unevenness = normal_magnitude * numpoints_inverse;
  surface_density = numpoints / cell->area_inverse;
  z_diff = max-min;
  z_max = max;
  z_min = min;
  z_ave /= numpoints;

  double radius1 = (max_intensity_point - min_intensity_point).norm();
  float sum_intensity1 = 0;
  float sum_intensity2 = 0;
  int count_1 = 0;
  int count_2 = 0;
  for (auto p_idx : cell->points_idx) {

      float diff = intensities[p_idx] - intensity_ave;
      intensity_variance += diff * diff;
      p = &(points[p_idx]);
      double z_current = SCALAR_PRODUCT_2p((*p), scene_normal) * scene_normal(2);
      double z_diff = z_current - z_ave;
      z_variance += z_diff * z_diff;

      Eigen::Vector3d diff1 = points[p_idx] - max_intensity_point;
      Eigen::Vector3d diff2 = points[p_idx] - min_intensity_point;
      if (diff1.norm() <= radius1) {
      count_1++;
      sum_intensity1 += intensities[p_idx];
      }
     if (diff2.norm() <= radius1) {
      count_2++;
      sum_intensity2 += intensities[p_idx];
      }

  }
  intensity_variance /= cell->points_idx.size();
  z_variance /= numpoints; 
  z_std = std::sqrt(z_variance); 
  radius1_intensity_distribution = (sum_intensity1 / static_cast<float>(std::max(count_1, 1)));
  radius2_intensity_distribution = (sum_intensity2 / static_cast<float>(std::max(count_2, 1)));
  return 1;
}






