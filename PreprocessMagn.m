function [ outPath, mask ] = PreprocessMagn( files, pathName, stages, outfile, v )

% FUNCTION:     PreprocessMagn
% DESCRIPTION:  Reads a series of dicom files (expected to be the magnitude
%               plots) and automatically creates the mask (indices).
%               The "magn" and "mask" variables are saved in a .mat file
%               in the parent directory (if the output file name is not
%               specified). For example, if pathName = "/home/b0/" then
%               outPath = "/home/b0.mat".
%               The "files" and "pathName" variables are required and can
%               be generated from the "uigetfile" Matlab function.
% INPUTS:       files    - list of dicom file names.
%               pathName - directory name of the files.
%               stages   - (optional) vector of flags selecting which
%                          stages of the preprocessing to run.
%                        - 1: Read in dicom files and save "magn" maps.
%                        - 2: Load magn maps and generate "masks".
%               outfile  - (optional) name of the output file. Cannot
%                          start with a number (if it does then the number
%                          is appended to the end of the name.
%               v        - (optional) verbosity flag:
%                          none (0), text (1), plots (2).
% OUTPUTS:      outPath  - Name of the output file.
%               mask     - Indices of the automatically generated voi mask.
% DEPENDENCIES: textprogressbar.m
%               read_dicom_header_local.m
%               CheckMatVars.m
%               PlotSlices.m
%               region_seg.m (bwdistsc.m)

%% Set default input arguments

outPath = []; mask = [];
if nargin < 5, v = 1;          end
if nargin < 4, outfile = [];   end
if nargin < 3 || isempty( stages ), stages = [1 1]; end
if nargin < 2, display( 'Error: Not enough inputs.' ); return; end
if ~iscell( files ) && isempty( outfile ) && stages(1) == 1
    if files == 0, return; else files = { files }; end
end

%% Generate output path (mat file)

if ~exist( 'outfile', 'var' ) || isempty( outfile )
    tokens = regexp( pathName, '\', 'split' );
    if length( tokens ) == 1, tokens = regexp( pathName, '/', 'split' ); end
    if length( tokens ) == 1, display( 'Warning: Could not save .mat file.' ); end
    pathstr = [pathName, '..'];
    fname   = tokens{ length( tokens ) - 1 };
    ext     = '.mat';
    clear tokens;
else
    [ pathstr, fname, ext ] = fileparts( outfile );
    if isempty( pathstr ), pathstr = '.'; end
end
fname = regexprep( fname, '\W', '' );
[~, idx] = regexp( fname, '^\d*W*' );
if ~isempty( idx ), fname = strcat( [ fname(idx+1:end), 'V' ], fname(1:idx) ); end
outPath = [ pathstr, '\', fname, ext ];

clear pathstr ext idx;

%% Read data from dicom files

% Run stage 1: "magn"
if stages(1) || ~CheckMatVars( outPath, 'magn' )
    
    % Iterate through files
    if v > 0, textprogressbar( 'Reading magnitude files: ',1 ); end
    for k = 1 : length( files )
        
        % Read dicom file
        if v > 0, textprogressbar( (k)/length(files) ); end
        obj = read_dicom_header_local( strcat( pathName, files{k} ) );
        
        % Store pixel data (overwrites with latest data)
        if isempty(obj.nslice )
            slice = obj.ima;
        elseif obj.nslice == 1
            slice = mod( obj.ima - 1, obj.nislab ) + 1;                     % 3D scans
        else
            slice = mod( obj.ima - 1, obj.nslice ) + 1;                     % 2D scans
        end
        magn( :,:,slice ) = obj.pix; clear slice;
        
    end
    if v > 0, textprogressbar( ' done.' ); end
    
    % Correct the orientation based on the header information
    m = matfile( outPath, 'Writable', true );
    eval( ['m.', fname, 'magn = magn;'] );
    if v > 1, PlotSlices( magn ); end

end

%% Generate mask from magnitude plots

% Run stage 2: "mask"
if stages(2) || ~CheckMatVars( outPath, 'mask' )
    
    if ~exist( 'm', 'var' ), m = matfile( outPath, 'Writable', true ); end
    eval( ['magn = m.', fname, 'magn;'] );
    ns = size( magn,3 );
    
    % Create initial mask for centre slice
    cslice = round( ns/2 );
    width  = round( 0.05*size( magn ) );
    mask   = ones( size( magn ) );
    mask( width(1):end - width(1), width(2):end - width(2), cslice ) = 0;
    
    % Mask the background
    if v > 0, textprogressbar( 'Creating masks: ',1 ); end
    k = cslice;
    mask(:,:,k) = region_seg( magn(:,:,k), mask(:,:,k), 1000, 3, false );   % Mask centre slice
    for k = ( cslice + 1 ) : ns                                             % Mask: Iterate forward
        if cslice == 1, continue; end
        if v > 0, textprogressbar( ( k - cslice )/ns ); end
        mask(:,:,k) = region_seg( magn(:,:,k), mask(:,:,k-1), 100, 3, false );
    end
    for k = cslice : - 1 : 1                                                % Mask: Iterate backwards
        if cslice == 1, continue; end
        if v > 0, textprogressbar( ( cslice - k )/ns + 0.5 ); end
        mask(:,:,k) = region_seg( magn(:,:,k), mask(:,:,k+1), 100, 3, false );
    end
    if v > 0, textprogressbar(1.0); textprogressbar( ' done.' ); end
    
    % Save the mask indices (inverted)
    [ ix, iy, iz ] = ind2sub( size( mask ), find( mask == 0 ) );
    mask = [ ix, iy, iz ];
    eval( ['m.', fname, 'mask = mask;'] );
    clear cslice width k ix iy iz;
    
    % Plot "mask"
    if v > 1
        maskplot = zeros( size( magn ) );
        idx = sub2ind( size( magn ), mask(:,1), mask(:,2), mask(:,3) );
        maskplot( idx ) = magn( idx ); 
        PlotSlices( maskplot );
    end

end

end