function [MOSO, MOSD, MOFO, MOFD, M1, M2, M3, M4, M5, M6, M7] = generateFRMatrixKernels_Symbolic(xOutr, xShap, xSoln, xFlux, nFluF, shapBasis, compBasis)

% -- This code generates matrix kernels **symbolically** for 2D FR

% -- The inputs are:

%       xOutr -- the coordinates of the outside points, useful for normals and integration area (anticlockwise numbering, must be convex)
%       xShap -- the coordinates of the points which will be used to geometrically fit the shape (so, might have 9 points for a P2 quad)
%       xSoln -- the coordinates of the solution points within the shape
%       xFlux -- the coordinates of the flux points, compatible with the order defined by nFluF
%       nFluF -- the number of flux points on each face (e.g., [5 1 5 1]). Must be compatible with xFlux and xOutr
%       shapBasis -- a symbolic column vector of the basis to be used for geometrically fitting the shape
%       compBasis -- a symbolic column vector of the basis to be used for determining the solution

% -- And the outputs are:
%       Lots of matrices for use as RBF kernels

% -- Extract all the symbolic expressions from the symbolic bases
Xsym = symvar([compBasis; shapBasis]);

% -- Get the number of dimensions from the outer points
nDim = size(xOutr,2);

% -- Extract the computational normals from xOutr
nFace = size(nFluF,2);
for iFace = 1:nFace

    % -- Work out which of the outer points bound this face
    lP = mod(iFace  -1,nFace)+1;
    rP = mod(iFace+1-1,nFace)+1;

    % -- Get the direction of the normal vector
    n(1) = (xOutr(rP,2) - xOutr(lP,2));
    n(2) = (xOutr(lP,1) - xOutr(rP,1));

    % -- And add in the normalised normals
    xNorm(sum(nFluF(1:iFace-1))+1:sum(nFluF(1:iFace)),1) = n(1) / norm(n);
    xNorm(sum(nFluF(1:iFace-1))+1:sum(nFluF(1:iFace)),2) = n(2) / norm(n);
end
fprintf('Computational normals successfully set\n')

% -- Set up the computational domain, useful for plotting easily
compPoly = polyshape(xOutr);
fprintf('Computational domain set\n')

% -- Integrate the "outer product" of the bases to form the Gramian matrix (can do this
%    numerically, but easy here to do this symbolically). Interesting to
%    note the use of the outer product to find the inner product of each!
gramShapO = polyInt(Xsym, shapBasis * shapBasis', xOutr);
gramCompO = polyInt(Xsym, compBasis * compBasis', xOutr);
fprintf('Gramian matrices calculated\n')

% -- Now use the Cholesky decomposition of the Gramian as a quick MGS
basisShapN = inv(chol(gramShapO))' * shapBasis;
basisCompN = inv(chol(gramCompO))' * compBasis;
fprintf('Basis functions orthonormalised over domain\n')

% -- Plot the domain and the flux and solution points, for a check
plot(compPoly); hold on; plot(xSoln(:,1), xSoln(:,2),'x'), plot(xFlux(:,1), xFlux(:,2), 'o')
fprintf('Flux and solution point locations set\n')

% -- Start building the matrices

% -- Compute the computational basis Alternant matrices
AlternantCSO = subs(basisCompN, {Xsym(1), Xsym(2)}, {xSoln(:,1)', xSoln(:,2)'});
AlternantCFO = subs(basisCompN, {Xsym(1), Xsym(2)}, {xFlux(:,1)', xFlux(:,2)'});
for iDim = 1:nDim
    AlternantCSD{iDim} = subs(diff(basisCompN,Xsym(iDim)), {Xsym(1), Xsym(2)}, {xSoln(:,1)', xSoln(:,2)'});
end
fprintf('Alternant matrices formed on the computational basis\n')

% -- Compute the outer basis Alternant matrices
AlternantSOO = subs(basisShapN, {Xsym(1), Xsym(2)}, {xShap(:,1)', xShap(:,2)'});
AlternantSSO = subs(basisShapN, {Xsym(1), Xsym(2)}, {xSoln(:,1)', xSoln(:,2)'});
AlternantSFO = subs(basisShapN, {Xsym(1), Xsym(2)}, {xFlux(:,1)', xFlux(:,2)'});
for iDim = 1:nDim
    AlternantSSD{iDim} = subs(diff(basisShapN,Xsym(iDim)), {Xsym(1), Xsym(2)}, {xSoln(:,1)', xSoln(:,2)'});
end
for iDim = 1:nDim
    AlternantSFD{iDim} = subs(diff(basisShapN,Xsym(iDim)), {Xsym(1), Xsym(2)}, {xFlux(:,1)', xFlux(:,2)'});
end
fprintf('Alternant matrices formed on the outer basis\n')

% -- Build the matrices which depend on the outer points - for setting up
%    solution and flux points, normals, jacobians, determinants...

% -- Solution point projection matrix (MOSO)
MOSO = cleanZeros(double((AlternantSOO) \ (AlternantSSO))');
fprintf('Solution point projection matrix found\n')

% -- Flux point project matrix (MOFO)
MOFO = cleanZeros(double((AlternantSOO) \ (AlternantSFO))');
fprintf('Flux point projection matrix found\n')

% -- Solution point metric projection matrix (MOSD)
for iDim = 1:nDim
    MOSDt{iDim} = cleanZeros(double((AlternantSOO) \ (AlternantSSD{iDim}))');
end
MOSD = interleave(MOSDt{:},'row');
fprintf('Gradient solution point projection matrix found\n')

% -- Flux point projection matrix (MOFD)
for iDim = 1:nDim
    MOFDt{iDim} = cleanZeros(double((AlternantSOO) \ (AlternantSFD{iDim}))');
end
MOFD = interleave(MOFDt{:},'row');
fprintf('Gradient flux point projection matrix found\n')

% -- Now build the matrix kernels which are run every iteration solve

% -- Set up the flux point primitive projection matrix (M1)
M1 = cleanZeros(double(((AlternantCSO) \ (AlternantCFO))'));
fprintf('Solution point to flux point primitive projection matrix found\n')

% -- Set up the flux divergence matrix (M2)
for iDim = 1:nDim
    M2t{:,iDim} = cleanZeros(double((AlternantCSO) \ (AlternantCSD{iDim}))');
end
M2 = interleave(M2t{:}, 'col');
fprintf('Solution point flux divergence matrix found\n')

% -- Set up the flux point flux projection matrix (M3)
for iDim = 1:nDim
    M3t{iDim} = cleanZeros(M1 .* xNorm(:,iDim));
end
M3 = interleave(M3t{:}, 'col');
fprintf('Solution point to flux point flux projection matrix found\n')

% -- Set up the correction function matrix (M4) -- by far the most complicated!
syms t real

% -- Express faces parametrically for line integration
for iFace = 1:nFace

    % -- Work out which of the outer points bound this face
    lP = mod(iFace  -1,nFace)+1;
    rP = mod(iFace+1-1,nFace)+1;

    % -- Put the face (line) in parametric form
    tT(1) = t * (xOutr(rP,1) - xOutr(lP,1)) + xOutr(lP,1);
    tT(2) = t * (xOutr(rP,2) - xOutr(lP,2)) + xOutr(lP,2);

    % -- Work out the dX/dt magnitude for line integral
    tJ = norm((xOutr(rP,:) - xOutr(lP,:)));

    % -- Find the parametric values corresponding to the flux points
    tFP = zeros(1,nFluF(iFace));
    for iFP = 1:nFluF(iFace)
        tFP(iFP) = solve(tT == xFlux(sum(nFluF(1:iFace-1))+iFP,:));
    end

    % -- Find the polynomials (in t) which a 1 at the relevant flux point and zero at the others
    pF = sym(zeros(1,nFluF(iFace)));
    for iFP = 1:nFluF(iFace)
        yFP = zeros(1,nFluF(iFace)); yFP(iFP) = 1;
        pF(iFP) = poly2sym(polyfit(double(tFP),yFP,nFluF(iFace)-1),t);
    end

    % -- Integrate over the system
    sTC(:,sum(nFluF(1:iFace-1))+1:sum(nFluF(1:iFace))) = int(subs(subs(basisCompN,Xsym(1),tT(1)),Xsym(2),tT(2))*pF,t,0,1)*tJ;
end
M4 = cleanZeros(double(AlternantCSO' * sTC));
fprintf('Correction function matrix found\n')

% -- Set up the solution point primitive differentiation matrix (M5)
M5 = interleave(M2t{:},'row');
fprintf('Solution point primitive differentiation matrix found\n')

% -- Set up the gradient projection matrix (M6)
for iDim = 1:nDim
    M6t{iDim} = M1;
end
M6t = blkdiag(M6t{:});
% -- Build a pair of fairly complex permutation vectors
for iDim = 1:nDim
    prVt{iDim} = 1:size(xFlux,1) + (iDim-1)*size(xFlux,1);
    pcVt{iDim} = 1:size(xSoln,1) + (iDim-1)*size(xSoln,1);
end
prV = interleave(prVt{:});
pcV = interleave(pcVt{:});
M6 = M6t(prV,pcV);
M6 = cleanZeros(M6);
fprintf('Solution point to flux point primitive gradient projection matrix found\n')

% -- Set up the gradient correction matrix (M7)
for iDim = 1:nDim
    M7((iDim-1)*size(xSoln,1)+1:iDim*size(xSoln,1),:) = repmat(xNorm(:,iDim)',size(M4,1),1) .* M4;
end
% -- Build a permutation vector
for iDim = 1:nDim
    prVt{iDim} = 1:size(xSoln,1) + (iDim-1)*size(xSoln,1);
end
prV = interleave(prVt{:});
M7 = M7(prV,:);
M7 = cleanZeros(M7);
fprintf('Primitive gradient correction matrix found\n')

end

% -- HELPER FUNCTIONS
% -- polyInt: symbolically integrate over any polynomial
function fI = polyInt(Xsym, f, xOutr)
%POLYINT This function integrates f symbolically over a polygon defined by xOutr

%     Xsym is a (row) vector of the symbols that appear in f
%     f is a symbolic array of any size, in terms of the variables in Xsym
%     xOutr is a set of points which define the domain of integration

% -- Start by triangulating xOutr
TR = delaunayTriangulation(xOutr);

% -- Zero the response matix
fI = sym(zeros(size(f,1), size(f,2)));

% -- Loop over all the triangles in the list, accumulating the integrals
for iT = 1:size(TR.ConnectivityList,1)
    fI = fI + triInt(Xsym, f, TR.Points(TR.ConnectivityList(iT,:),:));
end

end
% -- triInt: symbolically integrate over a triangle
function fI = triInt(Xsym, f, xTri)
%TRIINT This function integrates f symbolically over a triangle defined by xTri

% -- Set up some symbols
syms r s

% -- Calculate the change of variables
tX(1) = xTri(1,1) + (xTri(2,1) - xTri(1,1)) * r + (xTri(3,1) - xTri(1,1)) * s;
tX(2) = xTri(1,2) + (xTri(2,2) - xTri(1,2)) * r + (xTri(3,2) - xTri(1,2)) * s;

% -- Do substitution, replace x1 and x2 with s and t
fS = subs(f, {Xsym(1), Xsym(2)}, {tX(1), tX(2)});

% Get the determinant of the Jacobian of the transformation
dJ = diff(tX(1),r) * diff(tX(2),s) - diff(tX(1),s) * diff(tX(2),r);

% -- Integrate over the s and t triangle and rescale
fI = int(int(fS, s, 0, 1-r), r, 0, 1) .* dJ;

end
% -- CleanZeros: tidy up very small numbers from matrices
function M = cleanZeros(M)

    % -- This function just cleans up any very small values

    % -- Set the tolerance
    tol = 0.0000000001;

    % Zero the values which fall below this
    M(abs(M)<tol) = 0; 

end
% -- Interleave: thread matrices together row by row / column by column
function output = interleave(varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function interleaves the rows or columns of the specified vectors or
% matrices A, B, C, .... The last argument indicates whether to interleave
% the rows or columns by passing 'row' or 'col', respectively. If
% interleaving the rows, the number of columns of A and B must be the same;
% if interleaving the columns, the number of rows must be the same. 
%
% The order of input vectors/matrices indicate the order of interleaving
% (e.g. first row of A will be the first row, first row of B will be the
% second row of the interleaved matrix).
% 
% Extra rows or columns will be appended to the end of the final vector or
% matrix.
%
% If all inputs are vectors, their elements will be interleaved without
% the vectors needing to be all rows or all columns if the last argument
% ('row' or 'col') is not specified. If the orientation of the interleaving
% is specified, then the program will follow that orientation.
% 
% Examples:
% 1) Interleaving rows of matrices
% A = [1 2           B = [5 6
%      3 4]               7 8]     
% C = interleave(A, B, 'row')
% C = [1 2
%      5 6
%      3 4
%      7 8]
%
% 2) Interleaving columns of matrices
% C = interleave(A, B, 'col')
% C = [1 5 2 6
%      3 7 4 8]
%
% 3) Interleaving vectors
% A = [1 2 3 4]      B = [5 6 7 8 9]'
% C = interleave(A, B)
% C = [1 5 2 6 3 7 4 8 9]
%
% 4) Interleaving >2 matrices
% A = [1 2           B = [5 6         C = [9 10      D = [13 14
%      3 4]               7 8]            11 12]          15 16]
% E = interleave(A, B, C, D, 'col')
% E = [1 5  9 13 2 6 10 14
%      3 7 11 15 4 8 12 16]
% 
% 5) Interleaving columns of 2 matrices with unequal columns
% A = [1 2           B = [5 6  7  8
%      3 4]               9 10 11 12] 
% C = interleave(A, B, 'col')
% C = [1 5 2 6  7  8
%      3 9 4 10 11 12]
%     
% 6) Interleaving >2 vectors of unequal lengths
% A = [1 2 3 4]    B = [5 6 7]    C = [8 9 10 11 12 13]
% D = interleave(A, B, C, D)
% D = [1 5 8 2 6 9 3 7 10 4 11 12 13]
% 
% Written by Steven Shimizu
% March 4, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if interleaving orientation is specified
if ischar(varargin{end}) %check if last argument is a char array
    orientation = varargin{end};
    
    %Then check if interleaving orientation is anything but 'row' or 'col'
    if ~(isequal(orientation,'row') || isequal(orientation,'col'))
        %if not equal to 'row' or 'col', throw an error
        error('Last argument should specify ''row'' or ''col'' unless the inputs only consist of vectors.');
    else
        %The orientation is specified, then there are n-1 input numeric
        %arrays (single numbers, vectors or 2d matrices).
        numInputs = nargin - 1;
        
        for i = 1:1:numInputs
            if ~isnumeric(varargin{i})
                error('Inputs must be numeric arrays.');
            end
        end
    end
    
%check if last input argument is a numeric array, then check that all
%of the input arguments are vectors only, otherwise throw an error.
elseif isnumeric(varargin{end})
    orientation = 'none';
    numInputs = nargin;
    for i = 1:1:numInputs
        if ~isvector(varargin{i})
            error('The interleaving orientation ''row'' or ''col'' was not specified, which is only valid if all inputs are vectors.');
        end
    end
    
end
%Now that all inputs are known to be numeric arrays and interleaving
%orientation is specified, check that they are the right size
sizeArray = zeros(numInputs,2);
lengthList = zeros(numInputs,1);
for i = 1:1:numInputs
    sizeArray(i,:) = size(varargin{i});
    if isequal(orientation,'none')
        lengthList(i) = length(varargin{i});
    end   
end
if isequal(orientation,'row') && ~(sum(sizeArray(:,2)==sizeArray(1,2)) == numInputs)
    %check if number of columns are the same
    error('The number of columns of input matrices must be the same for interleaving the rows');
elseif isequal(orientation,'col') && ~(sum(sizeArray(:,1)==sizeArray(1,1)) == numInputs)
    error('The number of rows of input matrices must be the same for interleaving the columns');
end
output = [];        
%If only a single numeric array is passed, then return the same array
if numInputs == 1
    output = varargin{1};
elseif isequal(orientation,'none')
    %interleave vectors
    output = zeros(sum(lengthList),1);
    maxLength = max(lengthList);
    i_output = 1;
    for i = 1:1:maxLength
        for n = 1:1:numInputs
            if i <= lengthList(n)
                output(i_output) = varargin{n}(i);
                i_output = i_output + 1;
            end
        end
    end
elseif isequal(orientation,'row')
    %interleave by rows
    output = zeros(sum(sizeArray(:,1)),sizeArray(1,2));
    maxRows = max(sizeArray(:,1));
    i_output = 1;
    for i = 1:1:maxRows
        for n = 1:1:numInputs
            if i <= sizeArray(n,1)
                output(i_output,:) = varargin{n}(i,:);
                i_output = i_output + 1;
            end
        end
    end
elseif isequal(orientation,'col')
    %interleave by columns
    output = zeros(sizeArray(1,1),sum(sizeArray(:,2)));
    maxCols = max(sizeArray(:,2));
    i_output = 1;
    for i = 1:1:maxCols
        for n = 1:1:numInputs
            if i <= sizeArray(n,2)
                output(:,i_output) = varargin{n}(:,i);
                i_output = i_output + 1;
            end
        end
    end
end
end