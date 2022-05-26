% -- This is a driving script for generateFRMatrixKernels

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

% -- Set up the vector of symbols we need
Xsym = sym('x',[1 2], 'real');

% -- Generate the basis for the shape positioning and scaling (here, maximal order polynomials)
nOMaxOrder = round(sqrt(size(xOutr,1)));
shapBasis = @(x,d) squarePower(x,d,nOMaxOrder);

% -- Generate the basis for the computational behaviour (here, maximal order polynomials again)
nSMaxOrder = round(sqrt(size(xSoln,1)));
compBasis = @(x,d) squarePower(x,d,nSMaxOrder);

% -- Clean up some of the junk
clear nOMaxOrder nSMaxOrder

% -- Call the generating code
[MOSO, MOSD, MOFO, MOFD, M1, M2, M3, M4, M5, M6, M7] = generateFRMatrixKernels_Numeric(xOutr, xShap, xSoln, xFlux, nFluF, shapBasis, compBasis);


% -- Functions for generating the relevant function handles
% ----------------------------------------------------------
% ----------------------------------------------------------
% -- Functions for maximal order basis functions
function b = squarePower(x,d,n)

    nOMaxOrder = n;
    pow1D = 0:1:nOMaxOrder-1;
    xPow = reshape(repmat(pow1D',1,nOMaxOrder),[],1);
    yPow = reshape(repmat(pow1D ,nOMaxOrder,1),[],1);

    if d == 0
        b = (x(:,1)'.^(xPow) .* x(:,2)'.^(yPow));
    elseif d == 1   
        b = (max(0,xPow) .* x(:,1)'.^(max(0,xPow-1)) .* x(:,2)'.^(yPow));   
    elseif d == 2
        b = (x(:,1)'.^(xPow) .* max(0,yPow) .* x(:,2)'.^(max(0,yPow-1)));
    end

end
