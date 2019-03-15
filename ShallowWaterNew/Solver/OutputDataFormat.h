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
    OutputDataFormat() = default;
    OutputDataFormat(   unsigned long xSizeP,   unsigned long ySizeP,
                        long xMinP,             long xMaxP,
                        long yMinP,             long yMaxP) :
                        xSize(xSizeP),          ySize(ySizeP),
                        xMin(xMinP),            xMax(xMaxP),
                        yMin(yMinP),            yMax(yMaxP) {}
    ~OutputDataFormat() = default;

    //----------//

    void SetGridSize(unsigned long xSizeP, unsigned long ySizeP) {
        xSize = xSizeP;
        ySize = ySizeP;
    }
    void SetXRange(long xMinP, long xMaxP) {
        xMin = xMinP;
        xMax = xMaxP;
    }
    void SetYRange(long yMinP, long yMaxP) {
        yMin = yMinP;
        yMax = yMaxP;
    }
    void SetFrameNumber(unsigned long NumberP) {
        FrameNumber = NumberP;
        FrameNumberManual = true;
    }

    void AppendData(const std::vector <std::vector <T>>& DataP) {
        Data.push_back(DataP);

        if (!FrameNumberManual) {
            FrameNumber++;
        }
    }

    void Open(const std::string& PathP) {
        File.open(PathP);

        File    << xSize        << "\t"     << ySize    << "\t"
                << xMin         << "\t"     << xMax     << "\t"
                << yMin         << "\t"     << yMax     << "\t"
                << FrameNumber  << "\n";

    }
    void Close() {
        for (auto LineI : Data) {
            for (auto NumberI : LineI) {
                File << NumberI << "\t";
            }

            File << "\n";
        }

        File.close();
    }
private:
    unsigned long xSize = 1;
    unsigned long ySize = 1;

    unsigned long FrameNumber = 0;
    bool FrameNumberManual = false;

    long xMin = 0;
    long xMax = 1;

    long yMin = 0;
    long yMax = 1;

    std::ofstream File;
    std::vector <std::vector <T>> Data;
};
//-----------------------------//
#endif
