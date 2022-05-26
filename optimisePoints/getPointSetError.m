function [E] = getPointSetError(xEval, realBasis, testBasis, xIntg, wIntg, wT, p)

%GETPOINTSETERROR Generates the L2 error norm of fitting a basis using a set of points

% -------------------------------------------------------------

% -- The inputs are:

%       xEval -- the location of the points [nSoln -x- nDim]
%       realBasis -- a function handle which returns the evaluated basis to be used for the real fitting basis, as a function of x and d [= @(x,d) x(1) + x(2);]
%       testBasis -- a function handle which returns the evaluated basis to be used for the test basis, as a function of x and d [= @(x,d) x(1) + x(2);]
%       xIntg -- the coordinates of the integration points [nI -x- nDim]
%       wIntg -- the weights of the integration points [nI -x- 1]
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
AlternantSS = realBasis(xEval,0);
AlternantSI = realBasis(xIntg,0);
AlternantTI = testBasis(xIntg,0);

% -- If the condition number of the solution alternant is too high, return a huge number
if cond(AlternantSS) > 1e10
    E = 1e9;
    return
end

% -- Get the interpolating matrix from the solution (test) points to the integration points
M1 = AlternantSS \ AlternantSI;

% -- Evaluate the whole shebang of bases at the solution points
approxSolution = testBasis(xEval,0) * M1;

% -- Get the output error value
% So, (approxSolution-AlternantTI) is the bit inside the ()brackets in the paper (Eq. 19)
% To analyse properly, we multiply each row of this array by the relevant weighting function value, then add up along the columns
% Then we put each term to the p-norm power, and integrate the columns along the single remaining row
E = (sum(wT .* (approxSolution - AlternantTI),1).^p) * wIntg;

end
