#ifndef IMAGE_OPERATIONS_H
#define IMAGE_OPERATIONS_H

#include "../structures/outputchannelinfo.h"
#include "../io/imagefilename.h"

#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include <vector>
#include <string>

class Settings;

class ImageOperations {
 public:
  static void FitBeamSize(const Settings& settings, double& bMaj, double& bMin,
                          double& bPA, const aocommon::Image& image,
                          double beamEstimate);

  static void DetermineBeamSize(const Settings& settings, double& bMaj,
                                double& bMin, double& bPA, double& bTheoretical,
                                const aocommon::Image& image,
                                double initialEstimate);

  static void MakeMFSImage(const Settings& settings,
                           const std::vector<OutputChannelInfo>& infoPerChannel,
                           OutputChannelInfo& mfsInfo,
                           const std::string& suffix, size_t intervalIndex,
                           aocommon::PolarizationEnum pol,
                           ImageFilenameType image_type,
                           std::optional<size_t> directionIndex = std::nullopt);

  static void RenderMFSImage(const Settings& settings,
                             const OutputChannelInfo& mfsInfo,
                             size_t intervalIndex,
                             aocommon::PolarizationEnum pol, bool isImaginary,
                             bool isPBCorrected);
};

#endif
