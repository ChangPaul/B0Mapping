function mask = ReadVoi( filename, arrsize )

% FUNCTION:     ReadVoi
% DESCRIPTION:  Reads the vertices of a volume of interest (voi) in .obj
%               format. A 3d mask array is returned where the points
%               within (and on) the boundaries of the voi have value = 1.
%               The size of the mask is obtained from the header
%               information in the .obj file. If this is not avaible then
%               the maximum number of slices (and resolution) of the
%               vertices is used.
% INPUTS:       filename - full file name of the .obj file.
%               arrsize  - (optional) Size of the data (bounds on the
%                          maximum indices of the roi). Otherwise this is
%                          automatically determined from the header of the
%                          .obj file or if this is not available then the
%                          maximum indices of the roi are used.
% OUTPUTS:      mask     - 3d mask array of the points which lie inside the
%                          specified voi.

%% Check that the file exists
if nargin < 1, mask = []; return; end
if nargin < 2, arrsize = []; end
if ~exist( filename, 'file' ), display( 'ERROR: Roi file not found.' ); return; end

vertices = {};
ns = 0; res = [ 0, 0 ];

%% Read file
fid = fopen( filename );
tline = fgetl( fid );
while ischar( tline )
    
    % Separate comment from the line
    idx = find( tline == '#' );
    if ~isempty( idx )
        comment = tline( idx + 1 : end );
        tline = tline( 1 : idx - 1 );
    end
    
    % Check the comments for header information
    if length( comment ) > 0
        [ ftype comment ] = strtok( comment );
        if strcmpi( ftype, 'Slices:' )
            ns = str2num( comment );
        elseif strcmpi( ftype, 'Resolution:' )
            res = str2num( comment );
        end
    end
    
    % Read the roi properties
    if length( tline ) > 0
        [ ftype tline ] = strtok( tline, ' \t' );
        % Read vertices
        if ftype == 'v'
            c = str2num( tline );
            if length( vertices ) < c(3) + 1, vertices{ c(3)+1 } = []; end
            vertices{ c(3)+1 } = [vertices{ c(3)+1 }; c(1:2) ];
        end
    end
    tline = fgetl( fid );
end
fclose( fid );

%% Calculate the size of the mask array
cidx   = cellfun( @(x) ~isempty(x), vertices );
maxx   = max( cellfun( @(x) max( x(:,1) ), { vertices{ cidx } } ) );
maxy   = max( cellfun( @(x) max( x(:,2) ), { vertices{ cidx } } ) );
ns     = max( length( vertices ), ns );
res(1) = max( maxx, res(1) );
res(2) = max( maxy, res(2) );

%% Generate mask from voi vertices
voi = zeros( res(1), res(2), ns );
[ xq, yq ] = meshgrid( 1:res(1), 1:res(2) );
for s = 1 : min( ns, length( vertices ) )
    
    % Skip slices with no vertices
    if isempty( vertices{ s } ), continue; end
    
    % Find points within the polygon
    [ in, on ] = inpolygon( xq, yq, vertices{ s }(:,1), vertices{ s }(:,2) );
    voi(:,:,s) = in | on;
    
end

%% Resize the mask if the "arrsize" variable is specified

if length( arrsize ) > 2
    maskidx = min( arrsize(1:3), size(voi) );
    ix = 1:maskidx(1); iy = 1:maskidx(2); iz = 1:maskidx(3);
    temp = zeros( arrsize(1:3) );
    temp( ix,iy,iz ) = voi( ix,iy,iz );
    voi = temp; clear temp;
end

%% Return indices of mask

[ ix, iy, iz ] = ind2sub( size(voi), find( voi ) );
mask = [ ix, iy, iz ];

end