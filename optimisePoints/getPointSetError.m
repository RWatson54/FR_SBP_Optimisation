function [E] = getPointSetError(xOutr, xFlux, xEval, realBasis, testBasis, xIntg, wIntg, xIntF, nFluF, Msq, Msq_d, Mfq, N, wT, p)

%GETPOINTSETERROR Generates the L2 error norm of fitting a basis using a set of points

% -------------------------------------------------------------

% -- The inputs are:

%       xOutr -- the outer coordinates of the polygon [nOutr -x- nDim]
%       xShap -- the shape defining coordinates of the polygon [nShap -x- nDim]
%       xEval -- the location of the points [nSoln -x- nDim]
%       realBasis -- a function handle which returns the evaluated basis to be used for the real fitting basis, as a function of x and d [= @(x,d) x(1) + x(2);]
%       testBasis -- a function handle which returns the evaluated basis to be used for the test basis, as a function of x and d [= @(x,d) x(1) + x(2);]
%       xIntg -- the coordinates of the interior integration points [nII -x- nDim]
%       wIntg -- the weights of the interior integration points [nII -x- 1]
%       xIntF -- the coordinates of the face integration points [nIF -x- nDim]
%       wT   -- a weighting function for each of the test bases, [nTest -x- 1].
%       p -- the degree of the norm to be taken

% -- And the outputs are:
%       E -- the estimated error over the entire domain

% -------------------------------------------------------------

% -- If two of the points are extremely close, return a massive number
if min(pdist(xEval)) < 1e-2
    E = 1e10;
    return
end

% -- Generate the key alternant matrices from the supplied functions
AlternantSS = realBasis(xEval,0,xEval);
AlternantSI = realBasis(xIntg,0,xEval);
AlternantTI = testBasis(xIntg,0,xEval);

% -- If the condition number of the solution alternant is too high, return a huge number
if cond(AlternantSS) > 1e9
    E = 1e9;
    return
end

% -- If the SBP Property is missing, return a huge number
SBP = getSBP_External(xOutr, xEval, xFlux, xIntg, xIntF, nFluF, realBasis, Msq, Msq_d, Mfq, N);
if SBP > 1e-6
    E = 1e8;
    return
end

% -- Get the interpolating matrix from the solution (test) points to the integration points
M1 = AlternantSS \ AlternantSI;

% -- Evaluate the whole shebang of bases at the solution points
approxSolution = testBasis(xEval,0,xEval) * M1;

% -- Get the output error value
% So, (approxSolution-AlternantTI) is the bit inside the ()brackets in the paper (Eq. 19)
% To analyse properly, we multiply each row of this array by the relevant weighting function value, then add up along the columns
% Then we put each term to the p-norm power, and integrate the columns along the single remaining row
E = (sum(wT .* (approxSolution - AlternantTI),1).^p) * wIntg;

end
