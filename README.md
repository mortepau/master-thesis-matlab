# MATLAB scripts for "Implementation of Frequency-Offset Tolerant DFT-Based BFSK Demodulator"

This repository contains scripts used for verifying and testing a hardware
implementation of the BFSK demodulator described in the paper "On the Low
Complexity Implementation of the DFT-Based BFSK Demodulator for Ultra-Narrowband
Communications"

## Content

### eject.m

Eject computed data to a file, the format is determined based on the stage
we eject data from.

### inject.m

Inject precomputed data from a file, the file content is determined based
on the stage which we specify.

### generate_packet.m

Generates a packet for transmission by appending a preamble to it.

### generate_twiddle_factors.m

Generates a set of twiddle factors for the numerator N and denominator M.

### goertzel.m

Computes the DFT values of the samples given using the goertzel algorithm.

### rsdft.m

Computes the DFT values of the samples given using the rSDFT algorithm.
Possible to specify the r-factor, and the twiddle factors to use.

### tb_demodulator.m

The entrypoint for running simulations.
Contains the full structure/pipeline of the demodulator and allows
for ejecting and/or injecting precomputed data at any point.

### window_alignment_stage.m

The implementation of the window alignment stage used in the implementation
of the demodulator.

### zoom_stage.m

The implementation of the zoom stage used in the implementation of the
demodulator.
