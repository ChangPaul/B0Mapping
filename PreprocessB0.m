function [ outPath, b0 ] = PreprocessB0( files, pathName, stages, outfile, v )

% FUNCTION:     PreprocessB0
% DESCRIPTION:  Reads a series of dicom files (expected to be the phase
%               maps) and automatically creates the b0 maps.
%               The "phases", "info" and "b0" variables are saved in a
%               .mat file in the parent directory (if the output file name
%               is not specified). For example, if pathName = "/home/b0/"
%               then outPath = "/home/b0.mat".
%               The "files" and "pathName" variables are required and can
%               be generated from the "uigetfile" Matlab function.
% INPUTS:       files    - list of dicom file names.
%               pathName - directory name of the files.
%               stages   - (optional) vector of flags selecting which
%                          stages of the preprocessing to run.
%                        - (1): Read in dicom files and save "phases" maps
%                               and "info" header variable.
%                               skip (0), phase (1), phase unwrap (2)
%                        - (2): Load "phases" maps and generates "b0" maps.
%                               skip (0), b0 (1), b0 unwrap (2)
%               outfile  - (optional) name of the output file. Cannot
%                          start with a number (if it does then the number
%                          is appended to the end of the name.
%               v        - (optional) verbosity flag:
%                          none (0), text (1), plots (2).
% OUTPUTS:      outPath  - name of the output file.
% DEPENDENCIES: textprogressbar.m
%               read_dicom_header_local.m
%               CheckMatVars.m
%               PlotSlices.m
%               puma_ho.m (mf2.mexw64)

%% Set default input arguments

outPath = []; b0 = [];
if nargin < 5, v = 1; end
if nargin < 4, outfile = []; end
if nargin < 3 || isempty( stages ), stages = [2 2]; end
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

%% Read phase maps from dicom files

% Run stage 1: "phases" and "info"
if stages(1) || ~CheckMatVars( outPath, { 'phases', 'info' } )
    
    % Initialise header
    info.TE = [];    info.pos = [0 0 0];    info.scale = [1 1 1];
    
    % Read dicom files
    if v > 0, textprogressbar( 'Reading phase files: ',1 ); end
    for k = 1 : length( files )
        
        if v > 0, textprogressbar( (k)/length(files) ); end
        if ~iscell( files ), if v > 0, textprogressbar( ' done.' ); end, return; end
        obj = read_dicom_header_local( strcat( pathName, files{k} ) );
        if obj.ima == 1, info.pos = obj.pos; end
        
        % Store pixel data (overwrites with latest data)
        if isempty(obj.nslice)
            slice = obj.ima;
        elseif obj.nslice == 1
            slice = mod( obj.ima - 1, obj.nislab ) + 1;                     % 3D scans
        else
            slice = mod( obj.ima - 1, obj.nslice ) + 1;                     % 2D scans
        end
        phases( :,:,slice,obj.necho ) = obj.pix;
        info.TE( obj.necho ) = obj.TE;
        
    end
    if v > 0, textprogressbar( ' done.' ); end
    
    % Save header info
    if obj.ima > 1                                                          % Calculates slice dist
        slicedist = ( obj.pos(3) - info.pos(3) )/( slice - 1 );
        if slicedist < 0.1*obj.thickness, slicedist = obj.thickness; end
    else slicedist = obj.thickness;
    end
    if info.TE(1) == 0
        info.TE = info.TE( 2:end );
        phases( :,:,:,1 ) = [];
    end
    info.scale = [ obj.psize, slicedist ];
    info.MainOrientation = obj.MainOrientation;
    info.snorm = obj.snorm;
    info.band  = obj.band;
    info.freq  = obj.freq;
    clear slice slicedist obj files k;
    
    % Unwrap phase maps
    phases = phases*2*pi/4096;                                       % Convert phases to rad
    if stages(1) > 1
        [~,~,ns,np] = size( phases );
        if v > 0, textprogressbar( 'Unwrapping phase maps: ' ); end
        for p = 1 : np
            for s = 1 : ns
                if v > 0, textprogressbar( ( (p-1)*ns + s )/( np*ns ) ); end
                phases(:,:,s,p) = puma_ho( squeeze( phases(:,:,s,p) ), 2 );
            end
        end
        if v > 0, textprogressbar( ' done.' ); end
    end
    
    % Save "phases" and "info"
    m = matfile( outPath, 'Writable', true );
    eval( ['m.', fname, 'phases = phases;'] );
    eval( ['m.', fname, 'info   = info;'] );
    if v > 1, PlotSlices( phases ); end
end

%% Generate b0 maps from phase maps

% Run stage 2: "b0"
if stages(2) || ~CheckMatVars( outPath, 'b0' )
    
    if ~exist( 'm', 'var' ), m = matfile( outPath, 'Writable', true ); end
    eval( ['phases = m.', fname, 'phases;'] );
    eval( ['info   = m.', fname, 'info;'] );
    
    % Calculate b0 maps from phase maps
    b0 = zeros( size( phases ) ); b0 = b0( :,:,:,1 );
    switch( size( phases,4 ) )
        case 1
            b0 = phases; info.deltaTE(1) = info.TE(1);
        case 2
            b0 = phases(:,:,:,2) - phases(:,:,:,1);
            info.deltaTE(1) = info.TE(2) - info.TE(1);
        otherwise
            for k = 2 : size( phases,4 )
                b0(:,:,:,k-1) = phases(:,:,:,k) - phases(:,:,:,k-1);
                info.deltaTE(k-1) = info.TE(k) - info.TE(k-1);
            end
%             % Linearly interpolate
%             x = [ ones( size( phases,4), 1 ), info.TE' ];
%             info.deltaTE = info.TE( end ) - info.TE( 1 );
%             if v > 0, textprogressbar( 'Interpolating b0: ' ); end
%             for r = 1 : size( phases,1 )
%                 if v > 0, textprogressbar( r/size( phases,1 ) ); end
%                 for c = 1 : size( phases,2 )
%                     for s = 1 : size( phases,3 )
%                         y = unwrap( squeeze( phases( r,c,s,: ) ) );
%                         k = pinv( x'*x )*x'*y;
%                         b0( r,c,s ) = k(2);
%                     end
%                 end
%             end
%             b0 = b0*info.deltaTE;
%             if v > 0, textprogressbar( ' done.' ); end
%             clear r c s x y k;
    end
    
    % Unwrap the b0 maps
    if stages(2) > 1
        if v > 0, textprogressbar( 'Unwrapping b0 maps: ',1 ); end
        [ nr,nc,ns,np ] = size( b0 );
        for p = 1 : np
            
            % Unwrap the b0 maps
            for s = 1 : ns
                if v > 0, textprogressbar( ( (p-1)*ns + s )/( np*ns ) ); end
                b0(:,:,s,p) = puma_ho( squeeze( b0(:,:,s,p) ), 2 );
            end
            b0 = b0/2/pi;
            
            % Correct phase jumps
            ridx = floor( nr/3 ) : floor( nr/3*2);
            cidx = floor( nc/3 ) : floor( nc/3*2);
            for s = round( ns/2 ) : ns                                      % Phase jumps: iterate forwards
                meandiff = mean( mean( b0(ridx,cidx,s,p) - ...
                    b0(ridx,cidx,s-1,p) ) );
                if abs( meandiff ) > 0.5
                    b0(:,:,s,p) = b0(:,:,s,p) - round( meandiff );
                end
            end
            for s = round( ns/2 ) - 1 : -1 : 2                              % Phase jumps: iterate backwards
                meandiff = mean( mean( b0(ridx,cidx,s+1,p) - ...
                    b0(ridx,cidx,s,p) ) );
                if abs( meandiff ) > 0.5
                    b0(:,:,s,p) = b0(:,:,s,p) + round( meandiff );
                end
            end
            
        end
        if v > 0, textprogressbar( ' done.' ); end
    else
        b0 = b0/2/pi;
    end
    
    % Save "b0" and "info"
    eval( ['m.', fname, 'b0 = b0;'] );
    eval( ['m.', fname, 'info = info;'] );
    clear meandiff ridx cidx nr nc np s p;
    if v > 1 && ~isempty( b0 ), PlotSlices( b0(:,:,:,1) ); end
end

end