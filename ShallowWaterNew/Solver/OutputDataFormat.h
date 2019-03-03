#ifndef OUTPUTDATAFORMAT_H
#define OUTPUTDATAFORMAT_H
//-----------------------------//
#include <iostream>
#include <fstream>
#include <vector>
//-----------------------------//
template <typename T>
class OutputDataFormat {
public:
    explicit OutputDataFormat(std::string&& PathP) : Path(std::move(PathP)) {}
    ~OutputDataFormat() = default;

    //----------//

    void SetGridSize(unsigned long xSizeP, unsigned long ySizeP) {
        xSize = xSizeP;
        ySize = ySizeP;
    }
    void SetXRange(unsigned long xMinP, unsigned long xMaxP) {
        xMin = xMinP;
        xMax = xMaxP;
    }
    void SetYRange(unsigned long yMinP, unsigned long yMaxP) {
        yMin = yMinP;
        yMax = yMaxP;
    }
    void SetFrameNumber(unsigned long NumberP) {
        FrameNumber = NumberP;
    }

    void AppendData(const std::vector <std::vector <T>>& DataP) {
        Data.push_back(DataP);
    }

    void Open() {
        File.open(Path);

        File    << xSize    << "\t" << ySize    << "\t"
                << xMin     << "\t" << xMax     << "\t"
                << yMin     << "\t" << yMax     << "\t"
                << FrameNumber
                << std::endl;
    }
    void Close() {
        File.close();
    }
private:
    std::string Path;

    unsigned long xSize = 1;
    unsigned long ySize = 1;

    unsigned long FrameNumber = 1;

    unsigned long xMin = 0;
    unsigned long xMax = 1;

    unsigned long yMin = 0;
    unsigned long yMax = 1;

    std::ofstream File;
    std::vector <std::vector <T>> Data;
};
//-----------------------------//
#endif
