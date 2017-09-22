function coeffs = FieldToCoeffs( b0, pos, sphharm, amp, polyd, v )

% FUNCTION:     FieldToCoeffs
% DESCRIPTION:  Uses the "pos" and "b0" variables to calculate the
%               spherical harmonic decomposition coefficients.
%               Multiple sets of data can be decomposed. If b0 is a matrix
%               then each column c is a new set of b0 data corresponding to
%               the position matrix pos(:,:,c) with amplitude amp(c).
% INPUTS:       b0      - A [m x n] matrix where each column c is a vector
%                         of b0 values corresponding to the positions in
%                         pos(:,:,c) and amplitude amp(c). Therefore, there
%                         are n-sets of data to decompose.
%               pos     - A [m x 3 x n] matrix. n-sets of position data
%                         stored in each pos(:,:,c) OR if n=1 then the same
%                         set of position data is used to calculate the
%                         coeffs for each b0 vector.
%               sphharm - (optional) Indicies of spherical harmonic
%                         functions to be used for the decomposition.
%               amp     - (optional) A vector of length n that stores the
%                         applied amplitude.
%               polyd   - (optional) Simultaneously decomposes and fits a
%                         polynomial of degree "polyd - 1" to the data.
%               v       - (optional) verbosity flag:
%                          none (0), text (1), plots (2).
% OUTPUTS:      coeffs  - spherical harmonic coefficients.
% DEPENDENCIES: textprogressbar.m
%               spha.mexw64

% addpath( '../../../common' );

coeffs = [];
if nargin < 2, display( 'ERROR: Not enough inputs.' ); return; end
if size( pos,2 ) ~= 3, display( 'ERROR: Incorrect size of "pos" variable.' ); return; end

%% Set default input arguments

% Truncate "b0" and "pos" to have the same number of sets otherwise if
% there is only one set of "pos" then use this for each "b0" set
nsets = size( b0,2 );
if size( pos, 3 ) ~= 1
    nsets = min( size( b0,2 ), size( pos,3 ) );
    if length( amp ) > 1 && length( amp ) ~= nsets
        b0  = b0 ( :,  1:nsets );
        pos = pos( :,:,1:nsets );
        amp = amp( 1:nsets );
    end
end

% Check input arguments
if nargin < 3 || isempty( sphharm ), sphharm = 1:49; end                    % Full 6th order decomp
if nargin < 4 || isempty( amp ), amp = 1; end                               % Amplitudes = 1
if nargin < 5 || isempty( polyd ), polyd = 0; end
if nargin < 6 || isempty( v ), v = 1; end
if length( amp ) == 1, amp = amp*ones( 1, nsets ); end

%% Generate posvec and b based on the specified method

len = size( b0,1 );
numSH = length( sphharm );
b0 = b0./( ones( len,1 )*amp );

if polyd < 1 || nsets == 1
    %% Decompose separately
    
    if v > 0, textprogressbar( 'Decomposing separately: ',1 ); end
    coeffs = zeros( numSH, nsets );
    
    K = spha( sphharm, pos );
    for k = 1 : nsets
        if v > 0, textprogressbar( k/nsets ); end
        pset = k; if size( K, 3 ) == 1, pset = 1; end
        coeffs(:,k) = pinv( K(:,:,pset)'*K(:,:,pset) )*K(:,:,pset)'*b0(:,k);
    end
    if v > 0, textprogressbar( ' done.' ); end
    
elseif polyd == 1
    %% Stack "b0" and "pos" sets
    
    if size( pos, 3 ) ~= 1
        posvec = permute( pos, [1,3,2] );
        posvec = reshape( posvec, [ len*nsets, 3 ] );
    else
        posvec = zeros( len*nsets, 3 );
        for k = 1:nsets, posvec( ( k - 1 )*len + 1 : k*len,: ) = pos; end
    end
    A = spha( sphharm, posvec );
    coeffs = pinv( A'*A )*A'*b0(:);
    if v > 0, display( 'Decomposing stacked: done.' ); end

else
    %% Decompose by fitting a "polyd - 1" degree polynomial.
    
    % Partitions the data if it is too big.
    if v > 0, textprogressbar( 'Decomposing with interpolation: ',1 ); end
    nparts = ceil( numel( b0 )/1e6 );
    for p = 1 : nparts
        bp   = b0 ( p : nparts : end, : );
        posp = pos( p : nparts : end, :, : );
        len  = size( bp, 1 );

        A = zeros( len*nsets, numSH*polyd );
        K = spha( sphharm, pos );
        for k = 1 : nsets
            
            if v > 0, textprogressbar( ( (p-1)*nparts+k )/( nparts*nsets ) ); end

            % Concatenate "x^0", "x^1", ... "x^(polyd-1)"
            pset = k; if size( K, 3 ) == 1, pset = 1; end
            for j = 0 : polyd - 1
                A( (k-1)*len+1:k*len, j*numSH+1:numSH*(j+1) ) = K(:,:,pset)*( amp(k)^j );
            end
            
        end
        
        coeffs(:,:,p) = reshape( pinv(A'*A)*A'*bp(:), [ numSH, polyd ] );
    end
    coeffs = mean( coeffs, 3 );
    if v > 0, textprogressbar( ' done.' ); end
    
end

% Plot coefficients
if v > 1, figure; plot( coeffs ); end

end