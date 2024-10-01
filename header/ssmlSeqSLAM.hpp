/**
 * @class SSMLSeqSLAM
 * @brief A class for performing SeqSLAM on spherical images.
 *
 * This class provides methods to load datasets, preprocess images, and compute
 * SeqSLAM matches using different spherical representations. It supports three
 * types of spherical mappings: DIRECT_MAPPING, RESIZED_MAPPING, and VDSIL_MAPPING.
 *
 * The class uses LibPeR to perform the spherical image mapping.
 * The class uses OpenCV for image processing and Boost for filesystem operations.
 * It also integrates with the OpenSeqSLAM library for sequence matching.
 *
 * Key functionalities include:
 * - Loading datasets from CSV files.
 * - Extracting spherical intensities using different mapping methods.
 * - Computing SeqSLAM matches and saving results to CSV files.
 *
 */

#include <boost/filesystem/operations.hpp>
#include <iomanip>
#include <iostream>
#include <numeric>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <per/prRegularlySampledCSImage.h>
#include <per/prStereoModel.h>
#include <per/prStereoModelXML.h>
#include <per/prVoronoiIcosahedronImageMapping.h>

#include <boost/filesystem.hpp>

#include <visp/vpColVector.h>
#include <visp/vpDisplayX.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageTools.h>
#include <visp/vpTime.h>

#include "../header/OpenSeqSLAM.h"

#include <fstream>
#include <sstream>
#include <string>

#ifndef SSMLSEQSLAM_HPP
#define SSMLSEQSLAM_HPP

// Types of available spherical representations
typedef enum
{
  DIRECT_MAPPING,
  RESIZED_MAPPING,
  UNIPHORM_MAPPING
} prMethodType;

// Class for performing SeqSLAM on spherical images
class SSMLSeqSLAM
{
public:
  SSMLSeqSLAM() { isSphereInit = false; }
  void setSavePath(std::string path) { savePath = path; }

  void loadDatasetCV(std::string path)
  {

    std::string str_buf;
    std::string str_conma_buf;

    std::vector<std::vector<std::string>> img_paths;

    std::ifstream ifs_csv_file(path);

    while (getline(ifs_csv_file, str_buf))
    {
      std::istringstream i_stream(str_buf);

      std::vector<std::string> temp_line;

      while (getline(i_stream, str_conma_buf, ','))
      {
        temp_line.push_back(str_conma_buf);
      }
      img_paths.push_back(temp_line);
    }

    cv::Mat imgCurrCV;

    if (dataset_base_cv.size() > 0 || dataset_query_cv.size() > 0)
    {
      dataset_base_cv.clear();
      dataset_query_cv.clear();
    }

    int newWidth = 1000;
    int newHeight = 500;

    for (int i = 0; i < img_paths.size(); i++)
    {

      imgCurrCV = cv::imread(datasetPath + img_paths[i][0] + "/" +
                                 img_paths[i][2] + ".jpg",
                             cv::IMREAD_GRAYSCALE);
      cv::resize(imgCurrCV, imgCurrCV, cv::Size(newWidth, newHeight));
      dataset_base_cv.push_back(imgCurrCV);

      imgCurrCV = cv::imread(datasetPath + img_paths[i][1] + "/" +
                                 img_paths[i][3] + ".jpg",
                             cv::IMREAD_GRAYSCALE);
      cv::resize(imgCurrCV, imgCurrCV, cv::Size(newWidth, newHeight));
      dataset_query_cv.push_back(imgCurrCV);

      cv::imshow("img1", dataset_base_cv[i]);
      cv::imshow("img2", dataset_query_cv[i]);
      cv::waitKey(1);
    }
  }

  int computeOriginalSeqSLAM()
  {
    /* Find the matches */
    OpenSeqSLAM seq_slam;
    vector<cv::Mat> preprocessed_1 = seq_slam.preprocess(dataset_base_cv);
    vector<cv::Mat> preprocessed_2 = seq_slam.preprocess(dataset_query_cv);

    Mat matches = seq_slam.apply(preprocessed_1, preprocessed_2);

    std::string outResultFile =
        "../analysis/" + savePath + "/original/original_seqSLAM_equirect.csv";
    std::ofstream outResult;
    outResult.open(outResultFile);
    outResult << "score,index" << std::endl;

    char temp[100];
    float threshold = 0.9;

    float *index_ptr = matches.ptr<float>(0);
    float *score_ptr = matches.ptr<float>(1);

    double mean_score = 0.;
    double n = 0.;

    for (int x = 0; x < dataset_base_cv.size(); x++)
    {
      int index = static_cast<int>(index_ptr[x]);

      std::cout << score_ptr[x] << ", ";

      outResult << score_ptr[x] << "," << index << std::endl;
      if (x > 5 && x < dataset_base_cv.size() - 5)
      {
        mean_score += score_ptr[x];
        n++;
      }
    }

    std::cout << "mean score: " << mean_score / n << std::endl;

    outResult.close();
    return 0;
  }

  std::vector<Eigen::VectorXd>
  extractSphericalIntensities(std::vector<cv::Mat> set)
  {

    unsigned long nbSamples = 0;

    std::vector<Eigen::VectorXd> intensities;
    switch (this->methodType)
    {
    case DIRECT_MAPPING:
      nbSamples = delauMapper->nbSamples;
      break;

    case RESIZED_MAPPING:
      nbSamples = delauMapper->nbSamples;
      break;

    case VDSIL_MAPPING:
      nbSamples = voroMapper->nbSamples;
      break;
    }

    ecam.setPrincipalPoint(set[0].cols * 0.5, set[0].rows * 0.5);
    ecam.setPixelRatio(set[0].cols * 0.5 / M_PI,
                       set[0].rows * 0.5 / (M_PI * 0.5));

    Eigen::VectorXd localIntensity(nbSamples);

    vpImage<unsigned char> I_req_temp, Mask_temp;

    if (this->methodType == RESIZED_MAPPING && !isInit)
    {
      scaleFactor =
          sqrt((double)set[0].rows * (double)set[0].cols / (double)nbSamples);
      resizedWidth = round((double)set[0].cols / scaleFactor);
      resizedHeight = round((double)set[0].rows / scaleFactor);

      if (resizedHeight % 2 != 0)
      {
        resizedHeight++;
      }
      if (resizedWidth % 2 != 0)
      {
        resizedWidth++;
      }

      std::cout << "resize width x height " << resizedWidth << " x "
                << resizedHeight << std::endl;

      ecam.setPrincipalPoint(resizedWidth * 0.5, resizedHeight * 0.5);
      ecam.setPixelRatio(resizedWidth * 0.5 / M_PI,
                         resizedHeight * 0.5 / (M_PI * 0.5));

      Mask_temp.resize(resizedHeight, resizedWidth);

      vpImageTools::resize(Mask, Mask_temp, vpImageTools::INTERPOLATION_AREA);
      Mask.resize(resizedHeight, resizedWidth);
      Mask = Mask_temp;

      isInit = true;
    }

    for (unsigned long i = 0; i < set.size(); i++)
    {

      vpImageConvert::convert(set[i], I_req_temp);

      switch (this->methodType)
      {
      case DIRECT_MAPPING:
        delauMapper->buildFromEquiRect(I_req_temp, ecam, &Mask);
        break;

      case RESIZED_MAPPING:
        cv::resize(set[i], set[i], cv::Size(resizedWidth, resizedHeight));
        delauMapper->buildFromEquiRect(I_req_temp, ecam, &Mask);
        break;

      case UNIPHORM_MAPPING:
        voroMapper->buildFromEquiRect(I_req_temp, ecam, &Mask);
        break;
      }

      for (unsigned long s = 0; s < nbSamples; s++)
      {

        switch (this->methodType)
        {
        case DIRECT_MAPPING:
          localIntensity(s) = (double)delauMapper->bitmap[s];
          break;

        case RESIZED_MAPPING:
          localIntensity(s) = (double)delauMapper->bitmap[s];
          break;

        case UNIPHORM_MAPPING:
          localIntensity(s) = (double)voroMapper->bitmap[s];
          break;
        }
      }
      intensities.push_back(localIntensity);
    }

    return intensities;
  }

  int compute(int subdivLevel, prMethodType methodType)
  {

    std::cout << dataset_base_cv.size() << ", " << dataset_query_cv.size()
              << std::endl;
    std::cout << dataset_base_cv[0].size() << std::endl;

    std::cout << "subdivLevel: " << subdivLevel << std::endl;

    this->methodType = methodType;
    isInit = false;

    voroMapper =
        new prVoronoiIcosahedronImageMapping<unsigned char>(subdivLevel);
    delauMapper = new prRegularlySampledCSImage<unsigned char>(subdivLevel);

    Mask.resize(dataset_base_cv[0].rows, dataset_base_cv[0].cols);
    Mask = 255;

    std::vector<Eigen::VectorXd> intensities1 =
        extractSphericalIntensities(dataset_base_cv);
    std::vector<Eigen::VectorXd> intensities2 =
        extractSphericalIntensities(dataset_query_cv);

    std::cout << "spherical intensities extracted!" << std::endl;

    ecam.setPrincipalPoint(dataset_base_cv[0].cols * 0.5,
                           dataset_base_cv[0].rows * 0.5);
    ecam.setPixelRatio(dataset_base_cv[0].cols * 0.5 / M_PI,
                       dataset_base_cv[0].rows * 0.5 / (M_PI * 0.5));

    /* Find the matches */
    OpenSeqSLAM seq_slam;
    Mat matches = seq_slam.apply(intensities1, intensities2);

    std::cout << "matches extracted" << std::endl;

    std::string suffixPath;

    switch (methodType)
    {
    case DIRECT_MAPPING:
      boost::filesystem::create_directories("../analysis/" + savePath +
                                            "/direct");
      suffixPath =
          "direct/direct_mapping_" + std::to_string(subdivLevel) + "_sub";
      break;

    case RESIZED_MAPPING:
      boost::filesystem::create_directories("../analysis/" + savePath +
                                            "/resize");
      suffixPath =
          "resize/resize_mapping_" + std::to_string(subdivLevel) + "_sub";
      break;

    case VDSIL_MAPPING:
      boost::filesystem::create_directories("../analysis/" + savePath +
                                            "/voronoi");
      suffixPath =
          "voronoi/voronoi_mapping_" + std::to_string(subdivLevel) + "_sub";
      break;
    }

    std::string outResultFile =
        "../analysis/" + savePath + "/" + suffixPath + ".csv";
    std::ofstream outResult;
    outResult.open(outResultFile);
    outResult << "score,index" << std::endl;

    char temp[100];
    float threshold = 0.99;

    float *index_ptr = matches.ptr<float>(0);
    float *score_ptr = matches.ptr<float>(1);

    double mean_score = 0.;
    double n = 0;

    for (int x = 0; x < dataset_base_cv.size(); x++)
    {
      int index = static_cast<int>(index_ptr[x]);

      outResult << score_ptr[x] << "," << index << std::endl;
      if (x > 5 && x < dataset_base_cv.size() - 5)
      {
        mean_score += score_ptr[x];
        n++;
      }
    }

    std::cout << "mean score: " << mean_score / n << std::endl;

    outResult.close();
    return 0;
  }

  void setDataPath(std::string path) { datasetPath = path; }

private:
  void resizeStereoCam(prStereoModel &stereoCam, double scaleFactor)
  {
    for (int camNum = 0; camNum < 2; camNum++)
    {
      double u0 = ((prOmni *)(stereoCam.sen[camNum]))->getu0() / scaleFactor;
      double v0 = ((prOmni *)(stereoCam.sen[camNum]))->getv0() / scaleFactor;
      double au = ((prOmni *)(stereoCam.sen[camNum]))->getau() / scaleFactor;
      double av = ((prOmni *)(stereoCam.sen[camNum]))->getav() / scaleFactor;

      ((prOmni *)(stereoCam.sen[camNum]))->setPrincipalPoint(u0, v0);
      ((prOmni *)(stereoCam.sen[camNum]))->setPixelRatio(au, av);
    }
  }

  double scaleFactor;
  prMethodType methodType;
  std::string calibPath;
  prEquirectangular ecam;

  bool isInit, isSphereInit;
  int resizedWidth, resizedHeight;

  vpImage<unsigned char> Mask;
  prVoronoiIcosahedronImageMapping<unsigned char> *voroMapper;
  prRegularlySampledCSImage<unsigned char> *delauMapper;

  std::vector<cv::Mat> dataset_base_cv, dataset_query_cv;

  std::string savePath, datasetPath;
};

#endif // SSMLSEQSLAM_HPP
