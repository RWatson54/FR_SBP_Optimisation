% -- This is a driving script for generateFRMatrixKernels
% -- Because it's numerical, it relies on the getBasisFunctions script in ../basisFunctions

% -- Tidy up from a previous run
close all
clear all

% -- Set the coordinates of the outer points which define the shape of the element
xOutr = [                    1                  -1
                             1                   1
                            -1                   1
                            -1                  -1];

% -- Set the coordinates of the "shape" points which define the control points for positioning and scaling the element
xShap = [                    1                  -1
                             1                   1
                            -1                   1
                            -1                  -1];

% -- Set the solution point coordinates required
xSoln = [   -0.906179845938664  -0.906179845938664
            -0.906179845938664                   0
            -0.906179845938664   0.538469310105683
            -0.906179845938664   0.906179845938664
            -0.906179845938664  -0.538469310105683
            -0.538469310105683  -0.906179845938664
            -0.538469310105683   0.906179845938664
            -0.538469310105683  -0.538469310105683
            -0.538469310105683                   0
            -0.538469310105683   0.538469310105683
                             0  -0.906179845938664
                             0   0.906179845938663
                             0  -0.538469310105683
                             0   0.538469310105683
                             0                   0
             0.538469310105683  -0.538469310105683
             0.538469310105683                   0
             0.538469310105683   0.538469310105683
             0.538469310105683   0.906179845938663
             0.538469310105683  -0.906179845938664
             0.906179845938663  -0.906179845938663
             0.906179845938663                   0
             0.906179845938663   0.906179845938663
             0.906179845938663   0.538469310105683
             0.906179845938664  -0.538469310105683];

% -- Set up the distribution of points around the element
nFluF = [5 5 5 5];

% -- Set up the computational flux points
xFlux = [    1.000000000000000  -0.906179845938664
             1.000000000000000  -0.538469310105683
             1.000000000000000                   0
             1.000000000000000   0.538469310105683
             1.000000000000000   0.906179845938664
             0.906179845938664   1.000000000000000
             0.538469310105683   1.000000000000000
                             0   1.000000000000000
            -0.538469310105683   1.000000000000000
            -0.906179845938664   1.000000000000000
            -1.000000000000000   0.906179845938664
            -1.000000000000000   0.538469310105683
            -1.000000000000000                   0
            -1.000000000000000  -0.538469310105683
            -1.000000000000000  -0.906179845938664
            -0.906179845938664  -1.000000000000000
            -0.538469310105683  -1.000000000000000
                             0  -1.000000000000000
             0.538469310105683  -1.000000000000000
             0.906179845938664  -1.000000000000000];

% -- Set the basis type for the shape and computational basis
shapBasisType = 'Maximal2D';
compBasisType = 'Maximal2D';

% -- Call the generating code
[S1, C1] = getSBP_Numeric(xOutr, xShap, xSoln, xFlux, nFluF, shapBasisType, compBasisType);