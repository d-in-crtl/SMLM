### Single molecule localization analysis tools

These are `MATLAB` codes developed in the Rothenberg lab used across projects for analyses related to single molecule localization experiments. 

Authors of this code: 
1. Yandong Yin, New York University School of Medicine, New York, NY 10016, USA (present address: Shenzhen Bay Laboratory) [Main author]
2. Eli Rothenberg, New York University School of Medicine, New York, NY 10016 [For correspondence: eli.rothenberg@nyulangone.org]
 
#### Notes

The code is under development and is provided as is, and may contain certain duplicates or older file versions which will be cleaned up later. Codes include default parameters determined based on system configuations in the Rothenberg Lab at New York University School of Medicine. 
The rough organization is as follows:


```bash
2D_SMLM_Reconstruction #
ROIRender: # codes for rendering regions of interest
smPairCorrelation_CoordinateBased # codes for Auto and Cross pair correlation analyses
smTripleCorrelation_CoordinateBased # codes for Triple correlation analyses.
```

Codes were tested and executed in `MATLAB 2020b` and GPU Nvidia GTX 1060, CUDA 8.0 


### Codes used as per publication: 
```bash
Stepwise requirements for Polymerases δ and θ in Theta-mediated end joining. Stroik et al
Image reconstruction: 
1. 2D_SMLM_Reconstruction/main_2D_SCMOS.m

Auto and Cross pair corrlation analyses 
1. smPairCorrelation_CoordinateBased/smPCmain_v22112020.m
2. smPairCorrelation_CoordinateBased/smPCread_v22112020.m
3. smPairCorrelation_CoordinateBased/smPCmain_ForRNDv22112020.m
5. smPCread_ForRND.m

Triple pair correlation analyses
1. smPairCorrelation_CoordinateBased/smPCmain_ForTriv22112020.m
2. smTripleCorrelation_CoordinateBased/CoordTCmain_3CH_v171130.m
3. smTripleCorrelation_CoordinateBased/CoorTCAvg_vDelta.m
```
