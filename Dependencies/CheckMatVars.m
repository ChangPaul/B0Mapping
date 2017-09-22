function [ valid, fname ] = CheckMatVars( pathName, varList )

% FUNCTION:     CheckMatVars
% DESCRIPTION:  Check if the specified .mat file in "pathName" contains the
%               variables listed in "varList".
%               If the name of the file is "ref.mat" and the varList
%               contains a variable "b0" then it checks that ref.mat has
%               the variable "refb0" and returns [fname] = "ref".
% INPUTS:       pathName - name of the .mat file.
%               varList  - cell list of variables to check for. This
%                          variable can also be a single string if only one
%                          variable needs to be checked.
% OUTPUTS:      valid    - "true" if all the variables exist in the file.
%               fname    - name of the file.

% Check input arguments
valid = false; fname = [];
if nargin < 2, return; end
if ~iscell( varList ), varList = { varList }; end

% Check the file exists
if exist( pathName, 'file' )
    
    % Extract the file name (and remove punctuation)
    [ ~, fname, ~ ] = fileparts( pathName );
    fname = regexprep( fname, '\W', '' );
    [~, idx] = regexp( fname, '^\d*W*' );
    if ~isempty( idx ), fname = strcat( [ fname(idx+1:end), 'V' ], fname(1:idx) ); end

    % Check for the existence of the variables
    vars = whos( '-file', pathName );
    valid = true;
    for k = 1 : length( varList )
        varname = strcat( fname, varList{k} );
        if sum( arrayfun(@(x) strcmp( x.name, varname ), vars ) ) == 0
            valid = false;
            break;
        end
    end
end

end