function [ b0, pos ] = ProcessImage( fieldPath, refPath, def_voi, dTE, v )

% FUNCTION:     ProcessImage
% DESCRIPTION:  Calculates the b0 target vector (in Hz) and the position
%               matrix (Nx3) [in metres] for each of the files specifed
%               in "fieldPath". This file is expected to be .mat files
%               and requires the "b0" and "info" variables. If the
%               "def_voi" variable is not specified then the .mat files
%               should also contain the "mask" variable.
%               A reference b0 map can also be used if the "refPath" file
%               name is specified. This should be the file name of a .mat
%               file with the "b0" variable. This ref b0 map will only be
%               used if it is the same size as the "b0" variable of the b0
%               maps.
% INPUTS:       fieldPath- file name of the b0 field .mat file.
%               refPath  - (optional) file name of the ref .mat file.
%               def_voi  - (optional) User specified volume of interest.
%                          A Nx3xM size array with indices of the
%                          voi. For M=1, the columns are the indices of
%                          the field map. for M=2, the columns (optional)
%                          are the actual positions of the point.
%                          Overrides the automatically segmented voi
%                          ("mask" variable).
%               dTE      - (optional) Override the deltaTE in the header.
%               v        - (optional) verbosity flag:
%                          none (0), text (1), plots (2).
% OUTPUTS:      b0       - b0 vector (Hz)
%               pos      - Nx3 pos array (metres)
% DEPENDENCIES: textprogressbar.m
%               CorrectOrientation.m
%               CheckMatVars.m

%% Initialise variables and check files

% Set default input arguments
voi = [];
if nargin < 1, display( 'Error: Not enough inputs.' ); fieldPath = []; end
if nargin < 2, refPath = []; end
if nargin < 3, def_voi = []; end
if nargin < 4, dTE = 0.0;    end
if nargin < 5, v = 1;        end
if isempty( fieldPath ) | ( fieldPath == 0 ), return; end

if v > 0, textprogressbar( 'Calculating b0 and pos: ',1 ); end

% Read the reference b0 map from the specified path
ref = []; refmask = [];
if ~isempty( refPath )
    m = matfile( refPath );
    [ valid, fname ] = CheckMatVars( refPath, 'b0' );
    if valid, eval( ['ref = m.', fname, 'b0;'] ); end
    [ valid, fname ] = CheckMatVars( refPath, 'mask' );
    if valid, eval( ['refmask = m.', fname, 'mask;'] ); end
end

% Check for required variables in the file
[ valid, fname ] = CheckMatVars( fieldPath, { 'b0', 'info' } );
if ~valid, display( 'ERROR: Could not find "b0" or "info" in .mat file.' ); return; end

m = matfile( fieldPath );
eval( ['b0map = m.', fname, 'b0;' ] );
eval( ['info = m.', fname, 'info;'] );

if v > 0, textprogressbar( 0.25 ); end

%% Generate the "b0" vec

% Use default mask if the voi mask is not specified
if ~isempty( def_voi ), voi = def_voi;                                      % User voi
else
    if ~isempty( refmask )                                                  % Reference voi (supercedes field voi)
        nsize = size( b0map );
        if max( refmask ) <= nsize( 1:3 ), voi = refmask; end
    end
    clear refmask nsize;
end
if isempty( voi )
    if CheckMatVars( fieldPath, 'mask' )                                    % Field voi
        % Default mask was found
        eval( ['voi = m.', fname, 'mask;' ] );
    else                                                                    % Entire volume as voi
        % Otherwise use entire volume (not recommended)
        display( 'WARNING: Voi mask could not be found. Using entire volume.' );
        [ m, n, p ] = size( b0map );
        voi = ones( m, n, p );
        clear m n p;
    end
end
clear valid vars;

if v > 0, textprogressbar( 0.5 ); end

b0 = zeros( size( voi, 1 ), size( b0map, 4 ) );
for k = 1 : size( b0map, 4 )
    b0mapk = b0map(:,:,:,k);
    
    % Mask the b0
    if size( voi, 3 ) > 2 || size( voi, 2 )  > 3
        idx = sub2ind( size( b0mapk ), voi(:,1,1), voi(:,2,1), voi(:,3,1) );
        b0vec = b0mapk( idx );
    else
        voi = round( voi );
        b0vec = zeros( size( voi,1 ), 1 );
        for p = 1 : size( voi,1 )
            b0vec(p) = b0mapk( voi(p,1), voi(p,2), voi(p,3) );
        end
    end
    
    % Mask reference b0
    if ndims( ref ) == ndims( b0map ) && ndims( ref ) > 2
        if size( ref,4 ) == size( b0map, 4 )
            refk = ref(:,:,:,k);
        else
            refk = ref(:,:,:,1);
            display('WARNING: Using first ref map.');
        end
        % Use correct indices for refmap (if available)
        if size( voi, 3 ) < 3 || size( voi, 2 )  == 3
            tmp = zeros( size( voi,1 ), 1 );
            for p = 1 : size( voi,1 )
                tmp(p) = refk( voi(p,1), voi(p,2), voi(p,3) );
            end
            refk = tmp; clear tmp;
        else
            refk = refk( idx );
        end
        b0vec = b0vec - refk;
    end
    
    % Calculate resultant b0 in Hz
    if dTE > 0, deltaTE = dTE;
    elseif length( info.TE ) == 2, deltaTE = info.TE(2) - info.TE(1);
    elseif length( info.TE ) == size( b0map, 4 ) + 1, deltaTE = info.TE(k+1) - info.TE(k);
    elseif ~isempty( info.TE ), deltaTE = info.TE(1);
    else deltaTE = 1;
    end
    b0(:,k) = b0vec/( deltaTE*1e-3 );
end

if v > 0, textprogressbar( 0.75 ); end

%% Generate the "pos" matrix

% Calculate positions from positions indices and header info
num = size( b0vec, 1 );
pos = ( CorrectOrientation( voi( :,:,min( size( voi, 3 ), 2 ) ), ...        % Correct orientation
                            size( b0mapk ), info ).* ...
        ( ones( num, 1 )*info.scale ) + ...                                 % Scale pos
        ( ones( num, 1 )*info.pos ) )/1000;                                 % Add offset; convert to [m]

if v > 0, textprogressbar( 1.0 ); end
if v > 0, textprogressbar( ' done.' ); end

end
