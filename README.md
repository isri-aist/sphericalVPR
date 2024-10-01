# sphericalVPR

Visual Place Recognition SeqSLAM implementation for omnidirectional images mapped on the sphere

Based on the original works of *Milford et al.* [1]

This short example code shows how to use [LibPeR](https://github.com/PerceptionRobotique/libPeR_base) to map omnidirectional images on a spherical mesh and perform **visual place recognition**.

To download the example sequences, please follow the instructions available [here](https://github.com/isri-aist/MOVPR_downloader).

Available spherical representations are the followings:

- UniphorM
- Resized squares
- Direct projection

## How to use

Fill in required path in the example `main.cpp` file, then build the project using:

```bash
mkdir build
cd build
cmake ..
make
./main
```

## Dependencies

- Eigen
- Boost (for file and directory handling)
- OpenCV (tested with 4.2)
- LibPeR (tested with 0.5.0)
- ViSP (tested with 3.5)

## Credits

```
This software was developed at:
CNRS - AIST JRL (Joint Robotics Laboratory)
1-1-1 Umezono, Tsukuba, Ibaraki
Japan

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Description:
Implementation of the SeqSLAM algorithm for omnidirectional images mapped on the sphere

Authors:
Antoine ANDRE, Guillaume CARON

```

[1] Milford, M. J., & Wyeth, G. F. (2012, May). SeqSLAM: Visual route-based navigation for sunny summer days and stormy winter nights. In 2012 IEEE international conference on robotics and automation (pp. 1643-1649). IEEE.
