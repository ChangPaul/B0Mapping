addpath 'Dependencies';
addpath 'D:\MediaFire\PhD\ProjectManagement\Common';

refpath = 'C:\Users\changp\Desktop\shimfields\ref.mat';
voiPath = 'C:\Users\changp\Desktop\shimfields\roi.obj';

[fileName,pathName]  = uigetfile({'*.ima;*.IMA','IMA files (*.IMA)'},'Please select SIEMENS .IMA files','MultiSelect','on');
[outpath, b0 ] = PreprocessB0( fileName, pathName, [1 1], [] );

% voi = ReadVoi( voiPath, size( b0 ) );
% [ b0, pos ] = ProcessImage( outpath, refpath, voi, 0.76e-3 );
% 
% coeffs = FieldToCoeffs( b0, pos, 1:49, 0.5*42.576e6 );

% amp = -1514.38e-3;
% refPath = 'D:\MediaFire\PhdData\2015-11-24_2Dvs3D\shim terms\ref.mat';
% matPath = 'D:\MediaFire\PhdData\2015-11-24_2Dvs3D\shim terms\z2_N1514p38.mat';
% [ b0, pos ] = ProcessImage( matPath, refPath, [], 0.76e-3 );
% coeffs = FieldToCoeffs( b0, pos, 1:49, amp*42.576e6 );
% 
