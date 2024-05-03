# Unstable-Blowup-NLS
This repository contains the code for the paper "Mechanisms of unstable blowup in a quadratic nonlinear Schrödinger equation". Some parts are written in Mathematica, and other parts are written in MATLAB. 

Most of the MATLAB code uses the Chebfun package, which will need to be downloaded and added to MATLAB's path. See:
https://www.chebfun.org/docs/guide/

To run the numerical simulations described in the paper and store the data, one should run the following programs:

Within folder:  Matlab_Github/Blowup_and_relative_Modes/

<> Run 1st: spinOP_Fourier

<> Run 2nd: process_data

Within folder:  Matlab_Github/Self_Similar_Analysis/

<> Run 1st: spinOP_Fourier_LargerData

<> Run 2nd: spinOP_IntegrableBlowup


After one has created data files with the above method, then each of the figures in the paper were created with the respective files below: 

<> Figure 1 - Matlab_Github/Continuation_Search/spinOP_Fourier_Continuation

<> Figure 2 - Matlab_Github/Blowup_and_relative_Modes/Plotting_Chaotic_Blowup_vs2

<> Figure 3 - Mathematica/parameterization 2 Mode

<> Figure 4 - Matlab_Github/Blowup_and_relative_Modes/FourMode 

<> Figure 5&6 (a)   - Matlab_Github/Self_Similar_Analysis/Blowup_Profile_Simplify

<> Figure 5&6 (b)   - Matlab_Github/Self_Similar_Analysis/Plot_Norms

<> Figure 5&6 (c&d) - Matlab_Github/Self_Similar_Analysis/Self_Similar_2


Note that for making the figures in "Self_Similar_Analysis", you need to adust the code at the beginning of the file depending on whether you want to make FIgure 5 or 6.

