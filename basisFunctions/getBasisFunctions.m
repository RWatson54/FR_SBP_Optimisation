function [basis] = getBasisFunctions(basisType, nBasis, compPoly)

%GETBASISFUNCTIONS generates a function handle which evaluates the requested nBasis bases as a column vector at each of the input coordinates

% -------------------------------------------------------------

% -- The inputs are:

%       basisType -- The type of the basis that is being requested
%       nBasis -- how many bases of this type are being asked for
%       compPoly (optional) -- the polygon that we're using here, for the orthonormalisation

% -- The outputs are:

%       basis -- the function handle of the orthonormalised basis

% -------------------------------------------------------------

% -- Test to see if the basis number requested is reasonable
switch basisType

    case 'Maximal2D'

        % -- For the 2D maximal order basis, the requested number of bases should be a square number
        isSquare = mod(sqrt(nBasis),1)==0;
        if ~isSquare

            error('Incompatible number of bases requested for a Maximal2D basis')

        end

    case 'Total2D'

        % -- For the 2D total order basis, the requested number of bases should be a triangular number
        isTriangular = mod(sqrt(8*nBasis+1),1)==0;
        if ~isTriangular

            error('Incompatible number of bases requested for a Total2D basis')

        end

    case 'Euclid2D'

        % -- For the Euclidean basis, there's no easy closed form for the test, so list the numbers (OEIS A000603)
        isTriangular = ismember(nBasis,[1 3 6 11 17 26 35 45 58 73 90 106 123 146 168 193 216 243]);
        if ~isTriangular

            error('Incompatible number of bases requested for a Total2D basis')

        end
    otherwise

        error('No corresponding trial basis functions set in that size. Goodbye.\n')

end

% -- Set up the bases for the requested type
switch basisType

    case 'Maximal2D'

        % -- Set the 2D Maximal Order basis polynomial function handle
        basisO = @(x,d) maximalOrder2D(x,d,nBasis);

    case 'Total2D'

        % -- Set the 2D Total Order basis polynomial function handle
        basisO = @(x,d) totalOrder2D(x,d,nBasis);

    case 'Euclid2D'

        % -- Set the 2D Euclidean Order basis polynomial function handle
        basisO = @(x,d) euclidOrder2D(x,d,nBasis);
            
    otherwise

        error('No corresponding trial basis functions set. Goodbye.\n')

end

% -- Test to see if the polygon shape has been passed in for orthonormalisation
if exist('compPoly','var')
    % -- Generate the integration weights
    [xInt, wInt] = getIntegrationPoints(compPoly, 20);
    % -- Integrate the basis over the domain to get the Gramian
    Gm = (basisO(xInt,0) .* wInt') * basisO(xInt,0)';
    % -- And calculate the matrix for orthonormalisation with the inverse of the Cholesky
    oMat = inv(chol(Gm))';
else
    % -- Otherwise, set the orthonormalisation matrix to an identity
    oMat = eye(nBasis);
end

% -- Othonormalise the bases if required
basis = @(x,d) oMat * basisO(x,d);

end

% ---------------------------------------------
% ---------------------------------------------
%  Basis generating functions
% ---------------------------------------------
% ---------------------------------------------

% Feel free to add your own basis functions in here in the style of maximalOrder2D!

% -- Function for two dimensional maximal order basis
function [b] = maximalOrder2D(x,d,n)

    % -------------------------------------------------------------

    % -- The inputs are:

    %       x -- the coordinates to evaluate the basis at [nCoord -x- nDim]
    %       d -- the direction in which to differentiate, if any
    %       n -- the number of bases to compute

    % -- The outputs are:

    %       b -- the basis evaluated at the coordinate points

    % -------------------------------------------------------------

    % -- Get the powers involved
    nOMaxOrder = round(sqrt(n));
    pow1D = 0:1:nOMaxOrder-1;
    xPow = reshape(repmat(pow1D',1,nOMaxOrder),[],1);
    yPow = reshape(repmat(pow1D ,nOMaxOrder,1),[],1);

    % -- Compute the basis
    if d == 0
        b = (x(:,1)'.^(xPow) .* x(:,2)'.^(yPow));
    elseif d == 1   
        b = (max(0,xPow) .* x(:,1)'.^(max(0,xPow-1)) .* x(:,2)'.^(yPow));   
    elseif d == 2
        b = (x(:,1)'.^(xPow) .* max(0,yPow) .* x(:,2)'.^(max(0,yPow-1)));
    end

end

% -- Function for two dimensional total order basis
function [b] = totalOrder2D(x,d,n)

    % -------------------------------------------------------------

    % -- The inputs are:

    %       x -- the coordinates to evaluate the basis at [nCoord -x- nDim]
    %       d -- the direction in which to differentiate, if any
    %       n -- the number of bases to compute

    % -- The outputs are:

    %       b -- the basis evaluated at the coordinate points

    % -------------------------------------------------------------

    % -- Get the corresponding triangular number
    nOTotalOrder = (round(sqrt(8*n+1)) - 1) / 2;

    % -- Get the maximal powers involved
    nOMaxOrder = nOTotalOrder;
    pow1D = 0:1:nOMaxOrder-1;
    xPowM = reshape(repmat(pow1D',1,nOMaxOrder),[],1);
    yPowM = reshape(repmat(pow1D ,nOMaxOrder,1),[],1);

    % -- Truncate the maximal powers to get just the total order basis
    xPow = xPowM(xPowM + yPowM <= (nOTotalOrder-1));
    yPow = yPowM(xPowM + yPowM <= (nOTotalOrder-1));   

    % -- Compute the basis
    if d == 0
        b = (x(:,1)'.^(xPow) .* x(:,2)'.^(yPow));
    elseif d == 1   
        b = (max(0,xPow) .* x(:,1)'.^(max(0,xPow-1)) .* x(:,2)'.^(yPow));   
    elseif d == 2
        b = (x(:,1)'.^(xPow) .* max(0,yPow) .* x(:,2)'.^(max(0,yPow-1)));
    end

end

% -- Function for two dimensional Euclidean order basis
function [b] = euclidOrder2D(x,d,n)

    % -------------------------------------------------------------

    % -- The inputs are:

    %       x -- the coordinates to evaluate the basis at [nCoord -x- nDim]
    %       d -- the direction in which to differentiate, if any
    %       n -- the number of bases to compute

    % -- The outputs are:

    %       b -- the basis evaluated at the coordinate points

    % -------------------------------------------------------------

    % -- Get the corresponding triangular number

    % -- No easy closed form for this, so use find
    nOEuclidOrder = find([1 3 6 11 17 26 35 45 58 73 90 106 123 146 168 193 216 243] == n);
    
    % -- Get the maximal powers involved
    nOMaxOrder = nOEuclidOrder;
    pow1D = 0:1:nOMaxOrder-1;
    xPowM = reshape(repmat(pow1D',1,nOMaxOrder),[],1);
    yPowM = reshape(repmat(pow1D ,nOMaxOrder,1),[],1);

    % -- Truncate the maximal powers to get just the total order basis
    xPow = xPowM(xPowM.^2 + yPowM.^2 <= (nOEuclidOrder-1)^2);
    yPow = yPowM(xPowM.^2 + yPowM.^2 <= (nOEuclidOrder-1)^2);   

    % -- Compute the basis
    if d == 0
        b = (x(:,1)'.^(xPow) .* x(:,2)'.^(yPow));
    elseif d == 1   
        b = (max(0,xPow) .* x(:,1)'.^(max(0,xPow-1)) .* x(:,2)'.^(yPow));   
    elseif d == 2
        b = (x(:,1)'.^(xPow) .* max(0,yPow) .* x(:,2)'.^(max(0,yPow-1)));
    end

end