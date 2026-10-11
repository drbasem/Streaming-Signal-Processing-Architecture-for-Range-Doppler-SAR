# Streaming Signal-Processing Architecture for Range-Doppler SAR

This repository provides a MATLAB and Xilinx System Generator implementation of a streaming signal-processing architecture for the Range-Doppler Algorithm (RDA) used in Synthetic Aperture Radar (SAR) processing.

The repository includes both the conventional/classical processing flow and the proposed streaming processing flow. The goal is to provide a clear implementation reference showing how the Range-Doppler SAR processing chain can be organized as a streaming architecture while preserving the same signal-processing operations.

## Overview

The Range-Doppler Algorithm is commonly implemented as a sequence of frequency-domain processing stages, including range compression and azimuth compression. In this repository, two equivalent processing organizations are provided:

1. A classical MATLAB implementation of the Range-Doppler processing flow.
2. A streaming MATLAB implementation of the same processing chain.
3. A Xilinx System Generator model containing both the classical and streaming processing paths.

The MATLAB code and the System Generator model are developed in parallel so that their outputs can be compared directly.

## Main Idea

This project does not introduce a new SAR image-formation algorithm.  
Instead, it focuses on the signal-processing architecture and dataflow organization of the Range-Doppler processing chain.

The streaming architecture reorganizes the processing flow so that the SAR signal-processing stages can operate in a continuous streaming manner. The same main processing blocks are preserved:

- Range FFT
- Range matched filtering
- Range IFFT
- Azimuth FFT
- Azimuth matched filtering
- Azimuth IFFT

The emphasis of this repository is on the streaming organization of these processing stages and their implementation in MATLAB and Xilinx System Generator.

## Repository Contents

The shared files include:

- MATLAB code for the classical Range-Doppler processing flow
- MATLAB code for the streaming Range-Doppler processing flow
- Xilinx System Generator model including both the classical and streaming paths
- Supporting scripts and data needed to compare the MATLAB and System Generator outputs

## MATLAB Implementation

The MATLAB code contains both processing versions:

- Classical/reference Range-Doppler SAR processing
- Streaming-oriented Range-Doppler SAR processing

Both implementations are intended to process the same input data and generate comparable output results.

## System Generator Implementation

The Xilinx System Generator model contains hardware-oriented implementations of the classical and streaming processing paths.

The System Generator implementation is aligned with the MATLAB code so that the output of the hardware-oriented model can be compared against the MATLAB reference results.

## Output Comparison

After running the MATLAB code and the System Generator model, the generated outputs can be compared to verify that the streaming organization preserves the intended Range-Doppler processing behavior.

The comparison is intended to show consistency between:

- MATLAB classical processing output
- MATLAB streaming processing output
- System Generator classical processing output
- System Generator streaming processing output

## How to Use

1. Open the MATLAB project files.
2. Run the classical Range-Doppler processing script.
3. Run the streaming Range-Doppler processing script.
4. Open the Xilinx System Generator model.
5. Run the classical and streaming paths in System Generator.
6. Compare the MATLAB and System Generator outputs.

## Requirements

The implementation was developed for a MATLAB/Xilinx System Generator workflow.

Recommended tools:

- MATLAB
- Simulink
- Xilinx System Generator
- Xilinx Vivado

## Notes

This repository is intended as an engineering and research reference for streaming SAR signal-processing architecture.

The focus is on dataflow organization, streaming implementation, and consistency between MATLAB and System Generator models.

## Citation

If you use this code, architecture, figures, or numerical results in an academic work, please cite this repository or the accompanying technical report/preprint.

Suggested citation format:

Abbas Fadavi, "Streaming Signal-Processing Architecture for Range-Doppler SAR," GitHub repository, 2026.

## Author

Abbas Fadavi  
FPGA Hardware Engineer  
GitHub: https://github.com/abbasfadavi
