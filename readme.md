# Vitis HLS Implementation of 802.11a WLAN

This readme provides an overview of the Vitis HLS implementation of 802.11a WLAN. It covers the steps involved in implementing the 802.11a WLAN waveform using Vitis HLS and provides a reference to an example provided by MathWorks Inc.

## Overview
Vitis HLS (High-Level Synthesis) is a tool provided by Xilinx that allows designers to write C, C++, or SystemC code and convert it into RTL (Register Transfer Level) code. This enables faster development and optimization of hardware designs.

802.11a WLAN is a wireless communication standard that operates in the 5 GHz frequency band. It provides high data rates and is commonly used for applications such as video streaming and file transfer.

## Implementation Steps
To implement the 802.11a WLAN waveform using Vitis HLS, the following steps can be followed:

1. Design the waveform: Define the parameters and functionalities of the waveform, such as modulation scheme, coding scheme, and frame structure.

2. Write the C/C++ code: Use Vitis HLS to write the C/C++ code that implements the desired functionality of the waveform. This code will be converted into RTL code.

3. Verify and simulate: Use simulation tools to verify the functionality of the waveform implementation. This step helps identify and fix any issues before proceeding to synthesis.

4. Synthesize the design: Use Vitis HLS to synthesize the C/C++ code into RTL code. This step converts the high-level code into a gate-level representation that can be used for FPGA implementation.

5. Implement on FPGA: Use Xilinx Vivado or other FPGA implementation tools to implement the synthesized RTL code on an FPGA device.

6. Test and validate: Test the implemented design on the FPGA to ensure it meets the desired performance and functionality requirements.

## Example Reference
For a detailed example of implementing the 802.11a WLAN waveform using Vitis HLS, you can refer to the following website provided by MathWorks Inc.: [Image Transmission and Reception Using 802.11 Waveform SDR](https://au.mathworks.com/help/wlan/ug/image-transmission-reception-using-802-11-waveform-sdr.html)

This example provides step-by-step instructions, code snippets, and simulation results to guide you through the implementation process.
