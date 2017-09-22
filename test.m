addpath 'Dependencies';
voipath = '';                                % Insert filepath of VOI

% Process magnitude
[fileName,pathName]  = uigetfile({'*.ima;*.IMA','IMA files (*.IMA)'},'Please select SIEMENS .IMA files','MultiSelect','on');
outPath = PreprocessMagn( fileName1, pathName1 );

% Process B0
[fileName,pathName]  = uigetfile({'*.ima;*.IMA','IMA files (*.IMA)'},'Please select SIEMENS .IMA files','MultiSelect','on');
[outpath, b0 ] = PreprocessB0( fileName, pathName );

voi = ReadVoi( voiPath, size( b0 ) );
[ b0, pos ] = ProcessImage( outpath, [], voi, 0.76e-3 );

coeffs = FieldToCoeffs( b0, pos, 1:49, 0.5*42.576e6 );
