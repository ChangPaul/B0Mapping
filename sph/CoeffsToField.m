function [ recon, meas ] = CoeffsToField( coeffs, sphharm, pathOrPos, voiOrB0, units, v )

% FUNCTION:     CoeffsToField
% DESCRIPTION:  Generates the reconstructed b0 maps from the given
%               coefficients. These maps are compared to the actual
%               measured b0 field (if "fieldPath" is provided).
%               The default setting reconstructs on a 200x200x200mm FOV
%               with a spacing of 1mm.
% INPUTS:       coeffs  - A mxn size matrix where each column is a set of
%                         spherical harmonic coefficients.
%               sphharm - (optional) Indices of the spherical harmonic
%                         functions corresponding to each coefficient.
%               pathOrPos - (optional) file name of the b0 field that is
%                           to be compared
%                        OR this is a Nx3 array with the position of the
%                           physical coordinates of the b0 value to be
%                           reconstructed.
%                           If the value is a string then it's assumed to
%                           be the file name, otherwise if it is a numeric
%                           array then it's assumed to be positions.
%               voiOrB0 - (optional) A Nx3 matrix of indices of the volume
%                         of interest (pixels outside are masked)
%                       OR a vector of the b0 values corresponding to the
%                         positions given in "pathOrPos". This can also be
%                         a matrix with the same number of columns as the
%                         "coeffs" matrix. Reconstructed b0 values will be
%                         compared column-by-column.
%               units   - (optional) Units of measurement of coeffs.
%                          0 - m; 1 - dm; 2 - cm; 3 - mm.
%               v       - (optional) verbosity flag:
%                          none (0), text (1), plots (2).
% OUTPUTS:      recon   - Set of reconstructed b0 maps from each set of
%                         coefficients. Column c of coeffs corresponds to
%                         3d matrix recon(:,:,:,c).
%               meas    - Measures b0 field. If the b0 field is not provide
%                         (i.e. "fieldPath") then diff is filled with
%                         zeros.
% DEPENDENCIES: textprogressbar.m
%               spha.mexw64
%               CheckMatVars.m

%% Check input arguments and initialise variables

recon = []; meas = [];
if nargin < 1, return; end
if nargin < 2 || isempty( sphharm ), sphharm = 1 : size( coeffs, 1 ); end
if nargin < 3 || isempty( pathOrPos ), pathOrPos = ''; end
if nargin < 4, voiOrB0 = []; end
if nargin < 5 || isempty( units ), units = 0; end
if nargin < 6 || isempty( v ), v = 1; end
numSH = length( sphharm );

if( isnumeric( pathOrPos ) && ~isempty( voiOrB0 ) && ...
        length( voiOrB0 ) ~= size( pathOrPos, 1 ) ) || ...
        ( ~isnumeric( pathOrPos ) && ~isempty( voiOrB0 ) && ...
        size( voiOrB0, 2 ) ~= 3 )
    display( 'ERROR: Inconsistent inputs to CoeffsToField.' );
    return;
end

%% Calculate Nx3 position array

if v > 0
    textprogressbar( 'Reconstructing fields: ', 1 );
    textprogressbar( 0.2 );
end
if isnumeric( pathOrPos )
    %% Pos matrix from "pathOrPos" variable
    sz = size( pathOrPos );
    pos = zeros( sz(1), 3 );
    pos( :, 1:min( sz(2), 3 ) ) = pathOrPos( :, 1:min( sz(2), 3 ) );
    
    %% "b0" vector from "voiOrB0" variable (transpose if necessary)
    if ~isempty( voiOrB0 )
        if size( voiOrB0, 1 ) == sz(1), b0 = voiOrB0;
        elseif size( voiOrB0, 2 ) == sz(1), b0 = voiOrB0';
            display( 'ERROR: "pos" and "B0" are different lengths' );
            return;
        end
    end
    num = sz( 1 );
    
else
    %% Generate index meshgrid
    [ valid, fname ] = CheckMatVars( pathOrPos, 'b0' );
    if valid                                                                % Read size from file
        m = matfile( pathOrPos );
        eval( ['b0 = m.', fname, 'b0;'] );
        sz = size( b0 ); sz = sz( 1:3 );
    elseif ~isempty( voiOrB0 ), sz = max( voiOrB0 );                        % Obtain size from "voi"
    else sz = [ 201, 201, 201 ];                                            % Default size: 201x201x201
    end
    [ X, Y, Z ] = meshgrid( 1:sz(1), 1:sz(2), 1:sz(3) );
    
    %% Convert indices to physical coordinates
    
    % Read scale and offset from file
    minfo = [];
    if CheckMatVars( pathOrPos, 'info' )
        if ~exist( 'm', 'var' ), m = matfile( pathOrPos ); end
        eval( ['minfo   = m.', fname, 'info;'] );
        eval( ['scale  = minfo.scale;'] );
        eval( ['offset = minfo.pos;'] );
    else scale = [ 1 1 1 ]; offset = [ -100 -100 -100 ];                    % Default: "scale" 1mm
%     else scale = [ 1 1 1 ]; offset = [ -50 -50 -50 ];                       % Default: "scale" 1mm
    end                                                                     %          "offset" -100mm
    
    % Calculate pos
    num = numel( X );
    pos = CorrectOrientation( [Y(:) X(:) Z(:)], sz, minfo );
    pos = ( pos.*( ones( num, 1 )*scale ) + ...                             % Scale pos
        ( ones( num, 1 )*offset ) )*10^( units - 3 );                       % Add offset; convert to [units]
    clear scale offset;
    
    %% Initialise the voi mask
    if isempty( voiOrB0 )
        if CheckMatVars( pathOrPos, 'mask' )                                % Read voi from file
            if ~exist( 'm', 'var' ), m = matfile( pathOrPos ); end
            eval( ['voiOrB0 = m.', fname, 'mask;'] );
        else voiOrB0 = [ X(:), Y(:), Z(:) ];                                % Default voi: whole region
        end
    end
    mask = zeros( sz );
    for k = 1 : 3, voiOrB0( voiOrB0(:,k) > sz(k),: ) = []; end              % Remove out-of-bound voi points
    mask( sub2ind( sz, voiOrB0(:,1), voiOrB0(:,2), voiOrB0(:,3) ) ) = 1;
    clear X Y Z;
    
end

%% Calculate reconstructed b0 fields

if v > 0, textprogressbar( 0.4 ); end

% Initialise output variables
if isnumeric( pathOrPos )
    recon = zeros( sz(1), size( coeffs, 2 ) );
    meas  = zeros( sz(1), size( coeffs, 2 ) );
else
    recon = zeros( sz(1), sz(2), sz(3), size( coeffs, 2 ) );
    meas  = zeros( sz(1), sz(2), sz(3), size( coeffs, 2 ) );
end

% Reconstruct data for each set of coefficients or set of b0 maps
for k = 1 : size( coeffs, 2 )
    if v > 0, textprogressbar( 0.4 + k/size( coeffs, 2)*0.6 ); end
    
    % Reconstruct field
    if isnumeric( pathOrPos )
        recon(:,k) = spha( sphharm, pos )*coeffs(:,k);                      % Reconstructed
        if exist( 'b0', 'var' )
            meas(:,k)  = b0(:,k);                                           % Difference
            if k > size( b0, 2 ), break; end
        end
    else
        
        recon(:,:,:,k) = reshape( spha( sphharm, pos )*coeffs(:,k), sz ).*mask;
        if exist( 'b0', 'var' )
            deltaTE = 1;
            % Sets of coeffs == Sets of b0 maps (one-to-one comparison)
            if size( b0, 4 ) == size( coeffs, 2 )
                if exist( 'minfo', 'var' )
                    deltaTE = ( minfo.TE( k + 1 ) - minfo.TE(k) )*1e-3;
                end
                meas(:,:,:,k)  = b0(:,:,:,k)/deltaTE;
                
                % Compare nsets of recon maps to first b0 map
            elseif size( b0, 4 ) < size( coeffs, 2 )
                if exist( 'minfo', 'var' )
                    deltaTE = ( minfo.TE(2) - minfo.TE(1) )*1e-3;
                end
                meas(:,:,:,k)  = b0(:,:,:,1)/deltaTE;
                if size( b0, 4 ) ~= 1
                    display( 'WARNING: Using first b0 map for comparison' );
                end
                
                % Compare first recon map to nsets of b0 maps
            else
                for j = 1 : size( b0, 4 )
                    if exist( 'minfo', 'var' )
                        deltaTE = ( minfo.TE( j + 1 ) - minfo.TE(j) )*1e-3;
                    end
                    meas(:,:,:,j)  = b0(:,:,:,j)/deltaTE;
                end
                if size( coeffs,2 ) ~= 1
                    display( 'WARNING: Using first set of coefficients' );
                end
                break;
            end
        end
        
    end
end

if v > 0, textprogressbar( ' done.' ); end

%% Plot output
if v > 1
    if isnumeric( pathOrPos )
        if exist( 'b0', 'var' )
            figure; plot( [ recon(:,1), meas(:,1) ] );
        else
            figure; plot( recon(:,1) );
        end
    else
        PlotSlices( recon(:,:,:,1) );
        if exist( 'b0', 'var' ), PlotSlices( meas(:,:,:,1) ); end
    end
end

end