# ddp-shev

A repository for the Differential Dynamic Programming algorithm I developed for my PhD thesis
*Modern Dynamic Programming Algorithms for Optimal Control Problems in Hybrid Electric Vehicles*

## Content
This repository contains:
- Functions to run the differential dynamic programming algorithm described in the thesis.
- Two toy examples.
- An application to energy management strategy design of an hybrid electric vehicle.
- A draft of my PhD thesis. Part II describes the algorithm and the aforementioned application.

## Notes
The DDP algorithm depends on [CasADi](https://web.casadi.org/get/). You will need to download and install it to run the examples.

The HEV application also includes a script to solve the problem with (regular) dynamic programming. This script requires to install the [DynaProg toolbox](https://www.mathworks.com/matlabcentral/fileexchange/84260-dynaprog).

## Licensing
The contents of this repository are released under the [MIT license](LICENSE.md).

## Contact
For all questions or suggestions, feel free to contact me at federico.miretti@polito.it.