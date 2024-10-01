#include "header/ssmlSeqSLAM.hpp"
#include <vector>

int main(int argc, char **argv)
{

    SSMLSeqSLAM mapillary;

    std::string basePath = "path/to/sequences/ids/";

    // Paris streets
    std::string baseId = "i2PZwEzhAWag8tl9xJcQXo";
    std::string queryId = "gsOUFN18HjZlzKIhm4C7EX";

    mapillary.setDataPath("path/to/images/");
    mapillary.loadDatasetCV(basePath + baseId + "_" + queryId +
                            "_cured_data_pair_sequence.txt");

    // To perform SeqSLAM original implementation
    // mapillary.computeOriginalSeqSLAM();

    // To perform SeqSLAM on spherical images
    mapillary.compute(3, prMethodType::VDSIL_MAPPING);

    return 0;
}
