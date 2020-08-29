# Source Camera Verification from Strongly Stabilized Videos

This repository contains two newly generated datasets, pre-calculated transforms of VISION dataset, and the implementation for the proposed approach in Ref. [1].

## Implementation

You can open ‘VISIONDataSetConstrained3LevelHGS.m’ in Matlab. This function is designed for recreating results of VISION dataset in Ref. [1]. It also contains the implementation of the constrained 3 levels hierarchical grid search (HGS).

### Usage

- The function requires to functions and filters in [Goljan](http://dde.binghamton.edu/download/camera_fingerprint/ "Goljan, Fridrich"), and Matlab files PCE2.m and circxcorr2.m in [Hybrid_reference_for_videos.zip](https://lesc.dinfo.unifi.it/it/node/184)
- The function has 3 input values: 
     - visionPath: The path of the VISION folder in the [link](https://drive.google.com/file/d/10c7LEA_W2gm8TuD0GHccdV95epwLdAHd/view?usp=sharing). This folder contains pre-calculated 3-Level HGS transforms, corresponding video frames, and camera fingerprints.
     - isPreCalculated: If it is set 0, the code will calculate transform for frames and reference patterns. It can be used for running 3 levels HGS method for the other 2 datasets or new videos. If it is not set 0, the code will use pre-calculated transforms in the VISION folder.
     - isWeightExist: Due to the size of the weighting masks, they were not provided in the repository. If this parameter is set 0, the code will skip weighting, otherwise, it will read corresponding mat files for weighting. They can be created according to Ref. [2] by using [Weighter](https://github.com/VideoPRNUExtractor/Weighter).

### Example

When code run as `VISIONDataSetConstrained3LevelHGS('VISION/',1,0) `, it will calculate the number of verified videos and PCE threshold by using pre-calculated transforms in the VISION folder for 44 matching cases and 152 non-matching cases.

Getting results takes nearly 10 minutes for the average Linux computer. It will correctly verify 10 videos by using PCE value 52 as a threshold.
All non-matching tests were made by using the reference pattern of the previous camera in the dataset, except videos of the first camera which were tested with the reference pattern of the last camera.

## Datasets
3 different datasets which of detailed information can be found in Ref. [1] can be download from [link](https://drive.google.com/drive/folders/1ZOnvyOKRDQSmOcFADzK4mF4dsHn5sC0B?usp=sharing).

### VISION

VISION dataset contains the first 10 I-frame, except the first one, of the high or low stabilized videos of the VISION dataset in Ref. [3], their camera reference pattern, and pre-calculated 3 levels HGS transform of the matching case of highly stabilized videos and the non-matching case of all videos.

### iPhone SE-XR Dataset

This dataset includes media captured by 8 different phones (2 iPhone SE and 6 iPhone XR models) with a total of 263 photos and 41 videos (11-54 photos and 2-14 videos per camera). These media are collected by searching for public Flickr profiles and
using phones that we had access to.

### Adobe Premiere Pro Stabilized (APS) Video Dataset

This dataset contains 23 videos externally stabilized using Adobe Premiere Pro software suite, their corresponding reference patterns, and their 10 loop filter compensated frames that were used for testing constrained 3 level HGS.

## References

Ref [1]: E. Altinisik and H. T. Sencar, "Source Camera Verification for Strongly Stabilized Videos," in IEEE Transactions on Information Forensics and Security. [Link](https://ieeexplore.ieee.org/document/9169924) 

Ref [2]: E. Altınışık, K. Taşdemir, and H. T. Sencar, "Mitigation of H.264 and H.265 Video Compression for Reliable PRNU Estimation” IEEE Transactions on Information Forensics and Security, 2019. [.PDF](https://ieeexplore.ieee.org/document/8854840 ".PDF")

Ref [3]:Shullani, Dasara, Marco Fontani, Massimo Iuliani, Omar Al Shaya, and Alessandro Piva. "VISION: a video and image dataset for source identification." EURASIP Journal on Information Security 2017, no. 1 (2017): 15. [Link](https://link.springer.com/article/10.1186/s13635-017-0067-2)

