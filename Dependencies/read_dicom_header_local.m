function [m, off] = read_dicom_header(fname, out_level)
% ================================================================================  
%function m = read_dicom_header(fname, out_level)
%	fname:	filename of raw data file
%	out_level:	output informations
%			-3 data only, required tags (to do)
%			-2 no data, main tags
%			-1 no data
%			0 nothing
%			1 data of known tags
%			2 all tags
if nargin < 1,
   fname       = spm_get(1,'*','select dicom image');
end
if nargin < 2 
% out_flag=-1;	% read no data
 out_flag=0;	% nothing
% out_flag=1;	% data of known tags
% out_flag=2;	% all tags
else
	out_flag=out_level;
end
name = fname;
% disp(['out_flag: ' num2str(out_flag) num2str(out_flag < 0) num2str(out_flag==0)]);

%struct tag ( ...
%					'group', 0, ...
%					'member', 0, ...
%					'length', 0, ...			
%				};
%fid = fopen(fname, 'r', 'b');
fid = fopen(fname, 'r', 'l');
if fid < 0
        disp(['can''t open file ' fname]);
        return;
     end
%status=fstatus(fid);    
% fseek(fid,0,'eof'); flength=ftell(fid);fseek(fid,0,'bof');
fseek(fid,128,'bof');
kennung=char(fread(fid,4,'char')');
if isempty(kennung), m=[]; return, end
style_new=1;

if ~strcmp(kennung,'DICM')
	fseek(fid,0,'bof');
	style_new=0;
end
   
while ~feof(fid)
   group=fread(fid,1,'int16'); 
   member=fread(fid,1,'int16');
   if style_new
	   kennung=char(fread(fid,2,'char')');
	   if strcmp(kennung,'OB') | strcmp(kennung,'OW') | ...
	         strcmp(kennung,'SQ') | strcmp(kennung,'UN') | ...
	         strcmp(kennung,'UT')
	      length=fread(fid,1,'int16'); % dummy
	      length=fread(fid,1,'int32');
	   else
	      length=fread(fid,1,'int16');
	   end
   else
   	   kennung='XX';
	   length=fread(fid,1,'int32');
   end
%   fseek(fid,length,'cof');
   if (out_flag > 0), fprintf('\n');   end
   if (out_flag > 1), fprintf('  group: %4.4X member: %4.4X Kennung: %s length: %d ', ...
      								group, member, kennung, length);   end
%   if ((group == hex2dec('7FE0')) & (member==hex2dec('0010')))
%
%  evaluate different tags
%
   switch group 
   % =================================================================
     case hex2dec('0008'),		% Acquisition Group
		switch member
		 case hex2dec('0002'),	% Samples per Pixel
		   if ~style_new, kennung='IS'; end
	   	m.samples=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Samples per Pixel = %d',m.samples); end
	    case hex2dec('0008'),	% Acquisition Name
		   if ~style_new, kennung='CS'; end
	      m.acqname=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Acquisition Name = %s',m.acqname); end
      	    case hex2dec('0018') 
		   if ~style_new, kennung='UI'; end
              m.imgnumbr=get_dicom_value(fid,kennung,length);
	    case hex2dec('0033'),	% Image Date
		   if ~style_new, kennung='TM'; end
	      m.date=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Image Date = %s',m.date); end
	    case hex2dec('1030'),	% Study Description
		   if ~style_new, kennung='LO'; end
	      m.prog=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Study Description = %s',m.prog); end
	    case hex2dec('1090') 
		   if ~style_new, kennung='LO'; end
          	m.scanner=get_dicom_value(fid,kennung,length);    
	      if (out_flag > 0), fprintf('   -> Scanner Description = %s',m.scanner); end
	    case hex2dec('103E'),	% Series Description
		   if ~style_new, kennung='LO'; end
	      m.descript=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Series Description = %s',m.descript); end
	    case hex2dec('1140'),	% Seq Var
	      m.seqvar=get_dicom_value(fid,kennung,length);
	      %if (out_flag > 0), fprintf('   -> Seq Var = %s',m.seqvar); end
	      if (out_flag > 0), fprintf('   -> Seq Var'); end
	    otherwise,
		if (out_flag < -1)
		  fseek(fid,length,'cof');
		else
	      	  m.test=get_dicom_value(fid,kennung,length);
		  switch kennung
		     case {'CS','LO','SH','LT','IS','DS','UI','PN','AS','DA','DT','TM','ST','OB'}
		      	if (out_flag > 1), fprintf('   -> ??? = %s',m.test); end
		     otherwise
		      	if (out_flag > 1), fprintf('   -> ??? = %d',m.test); end
		  end
		end
	end
   % =================================================================
     case hex2dec('0010'),		% Patient Group
	   switch member
       case hex2dec('0010'),	% Patient's Name
	 off.pname=ftell(fid);
		   if ~style_new, kennung='PN'; end
         m.pname=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Patient''s Name = %s',m.pname); end
       case hex2dec('0020'),	% Patient ID
		   if ~style_new, kennung='LO'; end
	      m.pID=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Patient ID = %s',m.pID); end
       case hex2dec('0030'),	% Patient's Birth Date
		   if ~style_new, kennung='DA'; end
	      m.birth=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Patient''s Birth Date = %s',m.birth); end
       case hex2dec('0040'),	% Patient's Sex
		   if ~style_new, kennung='CS'; end
	      m.sex=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Patient''s Sex = %s',m.sex); end
	case hex2dec('1030'),	% Patient's Weight
		   if ~style_new, kennung='DS'; end
	   	m.weight=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Patient''s Weight = %d',m.weight); end
	case hex2dec('4000'),	% Additional Comment
	   	m.comment=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Additional Comment = %s',m.comment); end
	    otherwise,
		if (out_flag < -1)
			fseek(fid,length,'cof');
		else
	      	  m.test=get_dicom_value(fid,kennung,length);
		  switch kennung
		     case {'CS','LO','SH','LT','IS','DS','UI','PN','AS','DA','DT','TM','ST','OB'}
		      	if (out_flag > 1), fprintf('   -> ??? = %s',m.test); end
		     otherwise
		      	if (out_flag > 1), fprintf('   -> ??? = %d',m.test); end
		  end
		end
	end
   % =================================================================
     case hex2dec('0018'),		% Sequence Group
		switch member
		 case hex2dec('0020'),	% Scanning Sequence
		   if ~style_new, kennung='CS'; end
	   	m.scan=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Scanning Sequence = %s',m.scan); end
	    case hex2dec('0024'),	% Sequence Name
		   if ~style_new, kennung='CS'; end
	      m.seqname=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Sequence Name = %s',m.seqname); end
	    case hex2dec('0050'),	% Slice Thickness
		   if ~style_new, kennung='DS'; end
	      m.thickness=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Slice Thickness = %d',m.thickness); end
	    case hex2dec('0080'),	% Repetition Time
		   if ~style_new, kennung='DS'; end
	      m.TR=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Repetition Time = %d',m.TR); end
	    case hex2dec('0081'),	% Echo Time
		   if ~style_new, kennung='DS'; end
	      m.TE=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Echo Time = %d',m.TE); end
	    case hex2dec('0082'),	% Inversion Time
		   if ~style_new, kennung='DS'; end
	      m.TI=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Inversion Time = %d',m.TI); end
	    case hex2dec('0083'),	% Number of Averages
		   if ~style_new, kennung='DS'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Number of Averages = %d',m.test); end
	    case hex2dec('0084'),	% Imaging Frequency
		   if ~style_new, kennung='DS'; end
	      m.freq=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Imaging Frequency = %s',m.freq); end
	    case hex2dec('0086'),	% Echo Number
		   if ~style_new, kennung='IS'; end
	      m.necho=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Echo Number = %d',m.necho); end
	    case hex2dec('0088'),	% Slice Spacing
		   if ~style_new, kennung='DS'; end
	      m.spac=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Slice Spacing = %d',m.spac); end
	    case hex2dec('0089'),	% Number of Phase Encoding Steps
		   if ~style_new, kennung='IS'; end
	      m.phase=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Number of Phase Encoding Steps = %d',m.phase); end
	    case hex2dec('0093'),	% Percent Sampling
		   if ~style_new, kennung='DS'; end
	      m.psamp=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Percent Sampling = %d',m.psamp); end
	    case hex2dec('0094'),	% Percent Phase FoV
		   if ~style_new, kennung='DS'; end
	      m.pFoV=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Percent Phase FoV = %d',m.pFoV); end
	    case hex2dec('0095'),	% Pixel Bandwidth
		   if ~style_new, kennung='DS'; end
	      m.band=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Pixel Bandwidth = %d',m.band); end
	    case hex2dec('1020'),	% Software Version
		   if ~style_new, kennung='LO'; end
	      m.version=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Software Version = %s',m.version); end
	    case hex2dec('1310'),	% Acquisition Matrix
		   if ~style_new, kennung='IS'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Acquisition Matrix = %d',m.test); end
	    case hex2dec('1312'),	% Phase Encoding Direction
		   if ~style_new, kennung='CS'; end
	      m.pdir=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Phase Encoding Direction = %s',m.pdir); end
	    case hex2dec('1314'),	% Flip Angle
		   if ~style_new, kennung='DS'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Flip Angle = %d',m.test); end
	    case hex2dec('1318'),	% dB/dt
		   if ~style_new, kennung='DS'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> dB/dt = %d',m.test); end
	    otherwise,
		if (out_flag < -1)
		  fseek(fid,length,'cof');
		else
	      	  m.test=get_dicom_value(fid,kennung,length);
		  switch kennung
		     case {'CS','LO','SH','LT','IS','DS','UI','PN','AS','DA','DT','TM','ST','OB'}
		      	if (out_flag > 1), fprintf('   -> ??? = %s',m.test); end
		     otherwise
		      	if (out_flag > 1), fprintf('   -> ??? = %d',m.test); end
		  end
		end
	end
   % =================================================================
     case hex2dec('0020'),		% Data Group
		switch member
		 case hex2dec('0010'),	% Study ID
	        off.study=ftell(fid);
		   if ~style_new, kennung='SH'; end
	   	m.study=str2num(get_dicom_value(fid,kennung,length));
	      if (out_flag > 0), fprintf('   -> Study ID = %d',m.study); end
		 case hex2dec('0011'),	% Series Number
	        off.series=ftell(fid);
		   if ~style_new, kennung='IS'; end
	   	m.series=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Series Number = %d',m.series); end
		 case hex2dec('0012'),	% Acquisition Number
	        off.acq=ftell(fid);
		   if ~style_new, kennung='IS'; end
	   	m.acq=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Acquisition Number = %d',m.acq); end
		 case hex2dec('0013'),	% Image Number
	        off.ima=ftell(fid);
		   if ~style_new, kennung='IS'; end
	   	m.ima=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Image Number = %d',m.ima); end
	    case hex2dec('0032'),	% Image Position
		   if ~style_new, kennung='DS'; end
	      m.pos=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Image Position = %d',m.pos); end
	    case hex2dec('0037'),	% Image Orientation
		   if ~style_new, kennung='DS'; end
	      m.norm=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Image Orientation = %d',m.norm); end
	    case hex2dec('4000'),	% Image Comments
		   if ~style_new, kennung='CS'; end
	      m.comment=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Image Comments = %s',m.comment); end
	    otherwise,
		if (out_flag < -1)
		  fseek(fid,length,'cof');
		else
	      	  m.test=get_dicom_value(fid,kennung,length);
		  switch kennung
		    case {'CS','LO','SH','LT','IS','DS','UI','PN','AS','DA','DT','TM','ST','OB'}
		      	if (out_flag > 1), fprintf('   -> ??? = %s',m.test); end
		    otherwise
		      	if (out_flag > 1), fprintf('   -> ??? = %d',m.test); end
		  end
		end
	end
   % =================================================================
     case hex2dec('0028'),		% Image Presentation Group
	switch member
	   case hex2dec('0002'),	% Samples per Pixel
		   if ~style_new, kennung='US'; end
	      m.samples=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Samples per Pixel = %d',m.samples); end
	   case hex2dec('0004'),	% Photometric Interpretation
		   if ~style_new, kennung='CS'; end
	      m.photo=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Photometric Interpretation = %s',m.photo); end
	   case hex2dec('0008'),	% Number of Frames
		   if ~style_new, kennung='US'; end
	      m.frames=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Number of Frames = %d',m.frames); end
	   case hex2dec('0010'),	% Rows
		   if ~style_new, kennung='US'; end
	      m.rows=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Number of Rows = %d',m.rows); end
%	   case hex2dec('0029'),	% Cols
	   case hex2dec('0011'),	% Cols
		   if ~style_new, kennung='US'; end
	      m.cols=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Number of Cols = %d',m.cols); end
	   case hex2dec('0030'),	% Pixel Spacing
		   if ~style_new, kennung='DS'; end
	      m.psize=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Pixel Spacing = %d',m.psize); end
	   case hex2dec('0100'),	% Bits Allocated
		   if ~style_new, kennung='US'; end
	      m.nbits=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Bits Allocated = %d',m.nbits); end
	   case hex2dec('0101'),	% Bits Stored
		   if ~style_new, kennung='US'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Bits Stored = %d',m.test); end
	   case hex2dec('0102'),	% High Bit
		   if ~style_new, kennung='US'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> High Bit = %d',m.test); end
	   case hex2dec('0103'),	% Pixel Representation
		   if ~style_new, kennung='US'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> ??? = %d',m.test); end
	   case hex2dec('0106'),	% Smallest Image Pixel Value
		   if ~style_new, kennung='US'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> ??? = %d',m.test); end
	   case hex2dec('0107'),	% Largest Image Pixel Value
		   if ~style_new, kennung='US'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> ??? = %d',m.test); end
	   case hex2dec('1050'),	% Window Center
		   if ~style_new, kennung='DS'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> ??? = %s',m.test); end
	   case hex2dec('1051'),	% Window Width
		   if ~style_new, kennung='DS'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> ??? = %s',m.test); end
	   case hex2dec('1055'),	% Window Center & Width Explanation
		   if ~style_new, kennung='LO'; end
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> ??? = %s',m.test); end
	    otherwise,
		fseek(fid,length,'cof');
	end
    % =================================================================
     case hex2dec('0029'),		% Siemens Data Group
	switch member
	 case hex2dec('1010'),	% MEDCOM HEADER (Protokoll)
	 	if (length==5448), length=5444; end
	   	m.medcom=char(fread(fid,length,'uchar')');
	   	% medcom=char(fread(fid,length,'uchar')');
	        if (out_flag > 0), fprintf('   -> MEDCOM???'); end
	 case hex2dec('1020'),	% CSA HEADER
	   	m.csa=char(fread(fid,length,'uchar')');
		apos=findstr(m.csa,'### ASCCONV');
		if ~isempty(apos)
			m.ahead=m.csa(apos(1)+21:apos(2)-1);
		        if (out_flag > 0), fprintf('   -> ahead'); end
		end
	 case hex2dec('1110'),	% MEDCOM HEADER [csi]
	   	m.medcom=char(fread(fid,length,'uchar')');
	   	% medcom=char(fread(fid,length,'uchar')');
	        if (out_flag > 0), fprintf('   -> MEDCOM CSI'); end
	 case hex2dec('1120'),	% CSA HEADER [csi] (Protokoll)
	   	m.csa=char(fread(fid,length,'uchar')');
		apos=findstr(m.csa,'### ASCCONV');
		if ~isempty(apos)
			m.ahead=m.csa(apos(1)+21:apos(2)-1);
		        if (out_flag > 0), fprintf('   -> ahead'); end
		end
	 case hex2dec('1210'),	% MEDCOM HEADER [spec]
	   	m.medcom=char(fread(fid,length,'uchar')');
	   	% medcom=char(fread(fid,length,'uchar')');
	        if (out_flag > 0), fprintf('   -> MEDCOM SPEC'); end
	 case hex2dec('1220'),	% CSA HEADER [spec] (Protokoll)
	   	m.csa=char(fread(fid,length,'uchar')');
		apos=findstr(m.csa,'### ASCCONV');
		if ~isempty(apos)
			m.ahead=m.csa(apos(1)+21:apos(2)-1);
		        if (out_flag > 0), fprintf('   -> ahead'); end
		end
	    otherwise,
		if (out_flag < -1)
		  fseek(fid,length,'cof');
		else
	      	  m.test=get_dicom_value(fid,kennung,length);
		  switch kennung
		    case {'CS','LO','SH','LT','IS','DS','UI','PN','AS','DA','DT','TM','ST','OB'}
		      	if (out_flag > 1), fprintf('   -> ??? = %s',m.test); end
		    otherwise
		      	if (out_flag > 1), fprintf('   -> ??? = %d',m.test); end
		  end
		end
	end
  % =================================================================
     case hex2dec('6000'),		% Overlay Group
	switch member
	   case hex2dec('0010'),	% Overlay Rows
	      m.orows=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Number of Overlay Rows = %d',m.orows); end
	   case hex2dec('0011'),	% Overlay Cols
	      m.ocols=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Number of Overlay Cols = %d',m.ocols); end
	   case hex2dec('0015'),	% Number of Frames in Overlay
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   ->  Number of Frames in Overlay = %d',m.test); end
	   case hex2dec('0022'),	% Overlay Description
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   ->  Overlay Description = %s',m.test); end
	   case hex2dec('0040'),	% Overlay Type
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   ->  Overlay Type = %s',m.test); end
	   case hex2dec('0050'),	% Overlay Origin
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   ->  Overlay Origin = %d',m.test); end
	   case hex2dec('0051'),	% Image Frame Origin
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   ->  Image Frame Origin = %d',m.test); end
	   case hex2dec('0100'),	% Overlay Bits Allocated
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Overlay Bits Allocated = %d',m.test); end
	   case hex2dec('0102'),	% Overlay Bit Position
	      m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> Overlay Bit Position = %d',m.test); end
	   case hex2dec('3000'),	% Pixel Data
	      if (out_flag > 0), fprintf('Datenbereich (Overlay)!!! feof= %d',feof(fid)); end
	      if (out_flag >= 0)	      
		   %fseek(fid,-length,'eof');
		   %  m.Cols=sqrt(length/2);
		   %  m.Rows=m.Cols;
		      m.over=fread(fid,length*8,'ubit1');
		   %  m.pix=reshape(m.pix(:),m.Cols,m.Rows);
		      m.over=reshape(m.over(:),m.ocols,m.orows);
		else
		      fseek(fid,length,'cof');
		end
	      if (out_flag > 0), fprintf(' --> Datenbereich (Overlay)!!! feof= %d',feof(fid)); end
	      %break;
	    otherwise,
		if (out_flag < -1)
		  fseek(fid,length,'cof');
		else
	      	  m.test=get_dicom_value(fid,kennung,length);
		  switch kennung
		    case {'CS','LO','SH','LT','IS','DS','UI','PN','AS','DA','DT','TM','ST','OB'}
		      	if (out_flag > 1), fprintf('   -> ??? = %s',m.test); end
		    otherwise
		      	if (out_flag > 1), fprintf('   -> ??? = %d',m.test); end
		  end
		end
	end
  % =================================================================
     case hex2dec('7FE0'),		% Pixel Data Group
	switch member
	   case hex2dec('0010'),	% Pixel Data
	      m.offset=ftell(fid);
	      off.pix=ftell(fid);
	      if (out_flag > 0), fprintf('Datenbereich!!! feof= %d',feof(fid)); end
	      if (out_flag >= 0)	      
		   %fseek(fid,-length,'eof');
		   %  m.Cols=sqrt(length/2);
		   %  m.Rows=m.Cols;
		   if (m.nbits == 8) 
		      m.pix=fread(fid,length,'uint8');
		   else
		      m.pix=fread(fid,length/2,'int16');
		   end
		   %  m.pix=reshape(m.pix(:),m.Cols,m.Rows);
		   %  m.pix=reshape(m.pix(:),m.cols,m.rows);
		      m.pix=squeeze(reshape(m.pix,m.samples,m.cols,m.rows));
		      if m.samples > 1, m.pix=shiftdim(m.pix,1); end
		else
		      fseek(fid,length,'cof');
		end
	      if (out_flag > 0), fprintf(' --> Datenbereich!!! feof= %d',feof(fid)); end
	      break;
	    otherwise,
		fseek(fid,length,'cof');
	end
   % =================================================================
     case hex2dec('7FE1'),		% Spec Data Group
	switch member
           case hex2dec('0010'),	% Group Header
              m.test=get_dicom_value(fid,kennung,length);
	      if (out_flag > 0), fprintf('   -> ??? = %s',m.test); end
	   case hex2dec('1010'),	% Spec Data
	      if (out_flag > 0), fprintf('Datenbereich!!! feof= %d',feof(fid)); end
	      if (out_flag >= 0)	      
	 	  %m.spec=fread(fid,length/4,'float32');
	 	  m.spec=fread(fid,fix(length/8)*2,'float32');
		  m.spec=complex(m.spec(1:2:end),m.spec(2:2:end));
	      end
	      if (out_flag > 0), fprintf(' --> Datenbereich!!! feof= %d\n',feof(fid)); end
	      break;
	    otherwise,
		fseek(fid,length,'cof');
	end
   % =================================================================
     otherwise,
%    else
	if (out_flag < -1)
		fseek(fid,length,'cof');
	else
	     	m.test=get_dicom_value(fid,kennung,length);
		switch kennung
		   case {'CS','LO','SH','LT','IS','DS','UI','PN','AS','DA','DT','TM','ST','OB'}
		      	if (out_flag > 1), fprintf('   -> ??? = %s',m.test); end
		   otherwise
		      	if (out_flag > 1), fprintf('   -> ??? = %d',m.test); end
		end
	end
   end
   % if ((group == 2) & (member==1)), fseek(fid,6,'cof'); end
   % if ((group == 8) & (member==hex2dec('1140'))), fseek(fid,12,'cof'); end

%fseek(fid, 0, 'bof');  m.h_G08_Ide_StudyDate_Year = fread(fid, 1, 'long');
%fseek(fid, 820, 'bof');  m.h_G10_Pat_PatientSex = fread(fid, 1, 'int');
%fseek(fid, 856, 'bof');  m.h_G10_Pat_PatientSize = fread(fid, 1, 'double');
%fseek(fid, 1952, 'bof');  m.h_G19_Acq1_CM_ACScaleVector = fread(fid, 1, 'float');
%fseek(fid, 6058, 'bof');  t = fread(fid, 27, 'char')';  t(find(t > 127 | t < 0)) = 0;  m.h_G51_Txt_PatientName = char(t);

end
%gap
if isfield(m,'spac') & isfield(m,'thickness')
	m.gap=m.spac-m.thickness;
end

% ================================================================================  
handarbeit=1;
% ================================================================================  
%if (out_flag < -1), handarbeit=0; end;
if isfield(m, 'medcom')
	offset=findstr(m.medcom,'ICE_Dims');
	if ~isempty(offset),m.ICEDims=deblank(m.medcom(offset+100:offset+131)); end
end
if isfield(m, 'csa')
	test=get_ahead_value(m.csa,'sFastImaging.lEchoSpacing', out_flag);
	if ~isempty(test),m.EchoSpacing=test; end
end
if ~isfield(m, 'ahead'), handarbeit=0; end;
if handarbeit
%   fseek(fid, 0, 'bof');
%%   header=239914-320*320*2;
%   header=flength-m.rows*m.cols*2;
%	head=char(fread(fid,header,'char')');
%	apos=findstr(head,'### ASCCONV');
%	ahead=head(apos(1)+21:apos(2)-1);

%	npos=findstr(m.ahead,'sSliceArray.lSize');
%	tmp_str=m.ahead(npos+42:npos+42+5));
%	tmp_str=tmp_str(1:findstr(tmp_str,10)-1);
%	m.nslice=str2num(m.ahead(npos+42:npos+45));
	m.dwell=get_ahead_value(m.ahead,'sRXSPEC.alDwellTime[0]', out_flag);
	m.npart=get_ahead_value(m.ahead,'sKSpace.lPartitions', out_flag);
	m.nslice=get_ahead_value(m.ahead,'sSliceArray.lSize', out_flag);
	m.EPIFactor=get_ahead_value(m.ahead,'sFastImaging.lEPIFactor', out_flag);

	m.AcqMode=get_ahead_value(m.ahead,'sSliceArray.ucMode', out_flag);
	if ~isempty(m.AcqMode),
	 switch m.AcqMode
	  case '0x1',
		m.AcqMode='ASCENDING';
	  case '0x2',
		m.AcqMode='DESCENDING';
	  case '0x4',
		m.AcqMode='INTERLEAVED';
	  otherwise,
	 end
	end
	test=get_ahead_value(m.ahead,'sSliceArray.ucImageNumbSag', out_flag);
	if ~isempty(test),
	 switch test
	  case '0x0',
		m.NumMode='R>>L';
	  case '0x1',
		m.NumMode='L>>R';
	  otherwise,
	 end
	end
	test=get_ahead_value(m.ahead,'sSliceArray.ucImageNumbCor', out_flag);
	if ~isempty(test),
	 switch test
	  case '0x0',
		m.NumMode='A>>P';
	  case '0x1',
		m.NumMode='P>>A';
	  otherwise,
	 end
	end
	test=get_ahead_value(m.ahead,'sSliceArray.ucImageNumbTra', out_flag);
	if ~isempty(test),
	 switch test
	  case '0x0',
		m.NumMode='F>>H';
	  case '0x1',
		m.NumMode='H>>F';
	  otherwise,
	 end
	end
%	npos=findstr(m.ahead,'Thickness');
%	m.thickness=str2num(m.ahead(npos(1)+18:npos(1)+20));
	if ~isfield(m, 'thickness'), m.thickness=get_ahead_value(m.ahead,'Thickness', out_flag); end;

    m.specTE= get_ahead_value(m.ahead,'alTE[0]', out_flag);
    m.specTR= get_ahead_value(m.ahead,'alTR[0]', out_flag);

	m.mtrx       = get_ahead_value(m.ahead,'sKSpace.lBaseResolution', out_flag);
	m.mtry       = get_ahead_value(m.ahead,'sKSpace.lPhaseEncodingLines', out_flag);
	m.nislab     = get_ahead_value(m.ahead,'sKSpace.lImagesPerSlab', out_flag);
        %disp('line 634')
	m.ReadoutFoV = get_ahead_value(m.ahead,'sSliceArray.asSlice[0].dReadoutFOV', out_flag);
	m.PhaseFoV   = get_ahead_value(m.ahead,'sSliceArray.asSlice[0].dPhaseFOV', out_flag);
        %disp('line 637')
	%if ~isfield(m,'mtrx')
		if (isfield(m,'ReadoutFoV') & isfield(m,'psize'))
			m.mtrx=round(m.ReadoutFoV/m.psize(1));
		end
	%end
	%if ~isfield(m,'mtry')
		if (isfield(m,'PhaseFoV') & isfield(m,'psize'))
			m.mtry=round(m.PhaseFoV/m.psize(2));
		end
	%end
	test=get_ahead_value(m.ahead,'sSliceArray.asSlice[0].dInPlaneRot', out_flag);
	if ~isempty(test),m.SliceRot=test; end
	if isfield(m,'pdir')
	  if strcmp(m.pdir ,'ROW ')
		tmp=m.mtrx;
		m.mtrx=m.mtry;
		m.mtry=tmp;
	   end
	end
	test=get_ahead_value(m.ahead,'sTXSPEC.flReferenceAmplitude', out_flag);
	if ~isempty(test),m.tref=test; end
	test=get_ahead_value(m.ahead,'sWiPMemBlock.alFree[4]', out_flag);
	if ~isempty(test),m.supermosaic=test; end
	for p=1:14
		test=get_ahead_value(m.ahead,['sWiPMemBlock.alFree[' num2str(p) ']'], out_flag);
		if ~isempty(test),m.WiP_l(p)=test; end
	end
	for p=1:14
		test=get_ahead_value(m.ahead,['sWiPMemBlock.adFree[' num2str(p) ']'], out_flag);
		if ~isempty(test),m.WiP_d(p)=test; end
	end

    
    
    %vectorsize
    
    test=get_ahead_value(m.ahead,'sSpecPara.lVectorSize', out_flag);
    if ~isempty(test),m.vectorsize=test; end
	%
    %finalni velikost matice
    test=get_ahead_value(m.ahead,'sSpecPara.lFinalMatrixSizePhase', out_flag);
	if ~isempty(test),m.finalmatrixPhase=test; end
    test=get_ahead_value(m.ahead,'sSpecPara.lFinalMatrixSizeRead', out_flag);
	if ~isempty(test),m.finalmatrixRead=test; end
    test=get_ahead_value(m.ahead,'sSpecPara.lFinalMatrixSizeSlice', out_flag);
	if ~isempty(test),m.finalmatrixSlice=test; end



	%
	% Spectroscopy Voxel Position
	%
	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.sPosition.dSag', out_flag);
	if ~isempty(test),VposSag=test; else VposSag=0; end
	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.sPosition.dCor', out_flag);
	if ~isempty(test),VposCor=test; else VposCor=0; end
	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.sPosition.dTra', out_flag);
	if ~isempty(test),VposTra=test; else VposTra=0; end
	if max([VposSag VposCor VposTra]), m.VOIpos=[VposSag VposCor VposTra]; end
	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.sNormal.dSag', out_flag);
	if ~isempty(test),VnormSag=test; else VnormSag=0; end
	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.sNormal.dCor', out_flag);
	if ~isempty(test),VnormCor=test; else VnormCor=0; end
	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.sNormal.dTra', out_flag);
	if ~isempty(test),VnormTra=test; else VnormTra=0; end
	if max([VnormSag VnormCor VnormTra]), m.VOInorm=[VnormSag VnormCor VnormTra]; end
%	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.dThickness');
%	if ~isempty(test),m.VOIs=test; end
%	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.dPhaseFOV');
%	if ~isempty(test),m.VOIp=test; end
%	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.dReadoutFOV');
%	if ~isempty(test),m.VOIr=test; end
	VOIsizS=get_ahead_value(m.ahead,'sSpecPara.sVoI.dThickness', out_flag);
	VOIsizP=get_ahead_value(m.ahead,'sSpecPara.sVoI.dPhaseFOV', out_flag);
	VOIsizR=get_ahead_value(m.ahead,'sSpecPara.sVoI.dReadoutFOV', out_flag);
	if ~isempty(VOIsizS) & ~isempty(VOIsizP) & ~isempty(VOIsizR)
		%m.VOIsiz=[VOIsizS VOIsizP VOIsizR];
		m.VOIsiz=[VOIsizR VOIsizP VOIsizS];
	end
	test=get_ahead_value(m.ahead,'sSpecPara.sVoI.dInPlaneRot', out_flag);
	if ~isempty(test),m.VOIrot=test; end
	%
	% Diffusion Parameters
	%
	test=get_ahead_value(m.ahead,'sDiffusion.lDiffWeightings', out_flag);
	if ~isempty(test),m.DiffWeights=test; end
	test=get_ahead_value(m.ahead,'sDiffusion.lDiffDirections', out_flag);
	if ~isempty(test),m.DiffDirs=test; end
	% test=get_ahead_value(m.ahead,'sDiffusion.alBValue[1]', out_flag);
	if isfield(m,'seqname')
		test=strfind(m.seqname,'ep_b');
		if ~isempty(test),test=sscanf(m.seqname(test:end),'ep_b%d');end
		if ~isempty(test),m.b_value=test; end
		test=findstr(m.seqname,'#');
		if ~isempty(test), 
			test=sscanf(m.seqname(test:end),'#%d');
			if ~isempty(test), m.b_dir = test; end
		end;
	end;
	if isfield(m,'b_dir')
	   if isfield(m,'DiffDirs')
		switch m.DiffDirs
		case 6
			switch m.b_dir
				case 0, m.b_grad= [ 1  0  1];
				case 1, m.b_grad= [-1  0  1];
				case 2, m.b_grad= [ 0  1  1];
				case 3, m.b_grad= [ 0  1 -1];
				case 4, m.b_grad= [ 1  1  0];
				case 5, m.b_grad= [-1  1  0];
			end
		case 12
			if strcmp(m.version,'syngo MR B13 4VB13A ')
				switch m.b_dir
					case  1, m.b_grad= [ 1.000000,  0.414250, -0.414250 ];
					case  2, m.b_grad= [ 1.000000, -0.414250, -0.414250 ];
					case  3, m.b_grad= [ 1.000000, -0.414250,  0.414250 ];
					case  4, m.b_grad= [ 1.000000,  0.414250,  0.414250 ];
					case  5, m.b_grad= [ 0.414250,  0.414250,  1.000000 ];
					case  6, m.b_grad= [ 0.414250,  1.000000,  0.414250 ];
					case  7, m.b_grad= [ 0.414250,  1.000000, -0.414250 ];
					case  8, m.b_grad= [ 0.414250,  0.414250, -1.000000 ];
					case  9, m.b_grad= [ 0.414250, -0.414250, -1.000000 ];
					case 10, m.b_grad= [ 0.414250, -1.000000, -0.414250 ];
					case 11, m.b_grad= [ 0.414250, -1.000000,  0.414250 ];
					case 12, m.b_grad= [ 0.414250, -0.414250,  1.000000 ];
				end
			else
				switch m.b_dir
					case  0, m.b_grad= [  1   0  0.5];
					case  1, m.b_grad= [  0  0.5  1 ];
					case  2, m.b_grad= [ 0.5  1   0 ];
					case  3, m.b_grad= [  1  0.5  0 ];
					case  4, m.b_grad= [  0   1  0.5];
					case  5, m.b_grad= [ 0.5  0   1 ];
					case  6, m.b_grad= [  1   0 -0.5];
					case  7, m.b_grad= [  0 -0.5  1 ];
					case  8, m.b_grad= [-0.5  1   0 ];
					case  9, m.b_grad= [  1 -0.5  0 ];
					case 10, m.b_grad= [  0   1 -0.5];
					case 11, m.b_grad= [-0.5  0   1 ];
				end
			end
		case 30
			switch m.b_dir
				case  1, m.b_grad= [ -0.208098,  0.525514,  0.850005 ];
				case  2, m.b_grad= [  0.202387,  0.526131,  0.851002 ];
				case  3, m.b_grad= [  0.409956,  0.175267,  0.918257 ];
				case  4, m.b_grad= [ -0.412630,  0.742620,  0.565889 ];
				case  5, m.b_grad= [ -0.207127,  0.959492,  0.280092 ];
				case  6, m.b_grad= [ -0.872713,  0.525505,  0.064764 ];
				case  7, m.b_grad= [ -0.746815,  0.526129,  0.455449 ];
				case  8, m.b_grad= [ -0.415238,  0.175473,  0.915841 ];
				case  9, m.b_grad= [ -0.746636,  0.175268,  0.673642 ];
				case 10, m.b_grad= [ -0.665701,  0.742619, -0.217574 ];
				case 11, m.b_grad= [ -0.330391,  0.959489, -0.110458 ];
				case 12, m.b_grad= [ -0.331275,  0.525513, -0.809983 ];
				case 13, m.b_grad= [ -0.663936,  0.526130, -0.569521 ];
				case 14, m.b_grad= [ -0.999332,  0.175474, -0.111904 ];
				case 15, m.b_grad= [ -0.871398,  0.175267, -0.501922 ];
				case 16, m.b_grad= [  0.001214,  0.742616, -0.700356 ];
				case 17, m.b_grad= [  0.002949,  0.959483, -0.348370 ];
				case 18, m.b_grad= [  0.667975,  0.525509, -0.565356 ];
				case 19, m.b_grad= [  0.336490,  0.526126, -0.807431 ];
				case 20, m.b_grad= [  0.202383, -0.175470,  0.985002 ];
				case 21, m.b_grad= [  0.208094,  0.175265, -0.983848 ];
				case 22, m.b_grad= [  0.666452,  0.742619, -0.215262 ];
				case 23, m.b_grad= [  0.332212,  0.959489, -0.104850 ];
				case 24, m.b_grad= [  0.205064,  0.958364,  0.285421 ];
				case 25, m.b_grad= [  0.412630,  0.742620,  0.565889 ];
				case 26, m.b_grad= [  0.746093,  0.175315,  0.674232 ];
				case 27, m.b_grad= [  0.744110,  0.525505,  0.460568 ];
				case 28, m.b_grad= [  0.871894,  0.526125,  0.070507 ];
				case 29, m.b_grad= [  0.874264,  0.175471, -0.496841 ];
				case 30, m.b_grad= [  1.000000,  0.175267, -0.106112 ];
			end
		case 126
			load('/nrad/sekt/matlab5/126_dirs.mat');
			m.b_grad=test(:,m.b_dir)';
			%m.b_grad=grads2_icos(126);
			%m.b_grad=m.b_grad(:,m.b_dir)';
		case 252
			m.b_grad=grads2_icos(252);
			m.b_grad=m.b_grad(:,m.b_dir+1)';
		otherwise
		end
	   end
	end
	%
	% Slice Position
	%
	mitte1=0; mitte2=0;
	if (m.nslice > 1)
		mitte1=round(m.nslice/2)-1;
		mitte2=round(m.nslice/2)+1-mod(m.nslice,2)-1;
	end
	test=get_ahead_value(m.ahead,['sSliceArray.asSlice[' num2str(mitte1) '].sPosition.dSag'], out_flag);
	if ~isempty(test),posSag=test; else posSag=0; end
	test=get_ahead_value(m.ahead,['sSliceArray.asSlice[' num2str(mitte1) '].sPosition.dCor'], out_flag);
	if ~isempty(test),posCor=test; else posCor=0; end
	test=get_ahead_value(m.ahead,['sSliceArray.asSlice[' num2str(mitte1) '].sPosition.dTra'], out_flag);
	if ~isempty(test),posTra=test; else posTra=0; end
	if (mitte1 ~= mitte2)
		test=get_ahead_value(m.ahead,['sSliceArray.asSlice[' num2str(mitte2) '].sPosition.dSag'], out_flag);
		if ~isempty(test),posSag=(posSag+test)/2; end
		test=get_ahead_value(m.ahead,['sSliceArray.asSlice[' num2str(mitte2) '].sPosition.dCor'], out_flag);
		if ~isempty(test),posCor=(posCor+test)/2; end
		test=get_ahead_value(m.ahead,['sSliceArray.asSlice[' num2str(mitte2) '].sPosition.dTra'], out_flag);
		if ~isempty(test),posTra=(posTra+test)/2; end
	end
	m.spos=[posSag posCor posTra];
	%
	% Main Orientation
	%
	test=get_ahead_value(m.ahead,'sSliceArray.asSlice[0].sNormal.dSag', out_flag);
	if ~isempty(test),normSag=test; else normSag=0; end
	test=get_ahead_value(m.ahead,'sSliceArray.asSlice[0].sNormal.dCor', out_flag);
	if ~isempty(test),normCor=test; else normCor=0; end
	test=get_ahead_value(m.ahead,'sSliceArray.asSlice[0].sNormal.dTra', out_flag);
	if ~isempty(test),normTra=test; else normTra=0; end
	[normMax, ind]=max([normSag normCor normTra]);
	if (normMax == 0), ind=0; end
	switch ind
		case 1	
	 		m.MainOrientation='SAGITTAL';   ori='dSag'; m.snorm=[normSag normCor normTra];
		case 2	
	 		m.MainOrientation='CORONAL';    ori='dCor'; m.snorm=[normSag normCor normTra];
		case 3	
	 		m.MainOrientation='TRANSVERSE'; ori='dTra'; m.snorm=[normSag normCor normTra];
		otherwise
	 		m.MainOrientation='UNKNOWN';
	end
	if (m.nslice > 1)
	   if ~isfield(m,'gap')
% gap ist irgend wie immer falsch
% gap value geht so nur für transversaler Stapel !!!
%		npos=findstr(m.ahead,['Position.dTra']);
%		npos=findstr(m.ahead,['Position.' ori]);
%		if ~isempty(npos)
%		m.gap=m.thickness +str2num(m.ahead(npos(1)+18:npos(1)+30))...
%   					-str2num(m.ahead(npos(2)+18:npos(2)+30));
%		m.gap=abs(m.gap);
		%test1=get_ahead_value(m.ahead,['sSliceArray.asSlice[0].sPosition.' ori]);
		%test2=get_ahead_value(m.ahead,['sSliceArray.asSlice[1].sPosition.' ori]);
		test1=get_ahead_value(m.ahead,'sSliceArray.asSlice[0].sPosition.dSag', out_flag);
		test2=get_ahead_value(m.ahead,'sSliceArray.asSlice[0].sPosition.dCor', out_flag);
		test3=get_ahead_value(m.ahead,'sSliceArray.asSlice[0].sPosition.dTra', out_flag);
		test4=get_ahead_value(m.ahead,'sSliceArray.asSlice[1].sPosition.dSag', out_flag);
		test5=get_ahead_value(m.ahead,'sSliceArray.asSlice[1].sPosition.dCor', out_flag);
		test6=get_ahead_value(m.ahead,'sSliceArray.asSlice[1].sPosition.dTra', out_flag);
		if isempty(test1), test1=0; end
		if isempty(test2), test2=0; end
		if isempty(test3), test3=0; end
		if isempty(test4), test4=0; end
		if isempty(test5), test5=0; end
		if isempty(test6), test6=0; end
%		m.gap=str2num(m.ahead(npos(1)+18:npos(1)+30))...
%   					-str2num(m.ahead(npos(2)+18:npos(2)+30))
%		m.gap=abs(m.gap)-m.thickness;
		%if ~isempty(test1) & ~isempty(test2)
		%	m.gap=abs(test2-test1)-m.thickness;
		%if ~isempty(test1) & ~isempty(test2) & ~isempty(test3) & ~isempty(test4) & ~isempty(test5) & ~isempty(test6)
		if any([test1,test2,test3,test4,test5,test6])
			m.gap=sqrt(sum([test1-test4, test2-test5, test3-test6].^2))-m.thickness;
		else
			m.gap=0;
		end
		test1=get_ahead_value(m.ahead,'sGroupArray.asGroup[0].dDistFact', out_flag);
		if ~isempty(test1) & isfield(m,'thickness')
			m.gap=test1*m.thickness;
		end
		if abs(m.gap) < 1.0e-8, m.gap=0; end
	   end
	end;
end
fclose(fid);

% eof
	if ~isfield(m,'mtrx') & isfield(m, 'cols'), m.mtrx=m.cols; end
	if ~isfield(m,'mtry') & isfield(m, 'rows'), m.mtry=m.rows; end

% ================================================================================  
function v=get_dicom_value(fid,VR, length)
% ================================================================================  
	%fprintf('VR: %s, length: %d\n', VR,length);
%keyboard
  switch VR
%    case {'LT','UI','AS','DT','TM','ST'}
    case {'UI','AS','DT','ST'}
    	%ftell(fid)
    	%v=char(fread(fid,length,'char')');
    	v=fread(fid,length,'uchar');
    case {'CS','LO','SH','PN','DA','LT','TM'}
    	v=char(fread(fid,length,'uchar')');
    case {'IS','DS'}
%    	v=str2num(fread(fid,length,'uchar'));
    	v=char(fread(fid,length,'uchar')');
	%disp(['DS: ' v])
	v=str2num(strrep(v, '\',' '));
%    	v=str2num(char(fread(fid,length,'uchar')'));
    case 'SS',
%    	v=fread(fid,1,'int16');
    	v=fread(fid,length/2,'int16');
	if length ~= 2
		%fprintf('  ----> !!! unsigned short mit %d byte', length);
	    	%fseek(fid, length-2,'cof');
	end
    case 'US',
%    	v=fread(fid,1,'uint16');
    	v=fread(fid,length/2,'uint16');
	if length ~= 2
		%fprintf('  ----> !!! unsigned short mit %d byte', length);
	    	%fseek(fid, length-2,'cof');
	end
    case 'UL',
%    	v=fread(fid,1,'uint32');
    	v=fread(fid,length/4,'uint32');
	if length ~= 4
		%fprintf('  ----> !!! unsigned long mit %d byte', length);
	    	%fseek(fid, length-4,'cof');
	end
    case 'FL',
%    	v=fread(fid,1,'float32');
    	v=fread(fid,length/4,'float32');
    case 'FD',
%    	v=fread(fid,1,'float64');
    	v=fread(fid,length/8,'float64');
    case 'OB',			% other byte string
    	fseek(fid, length,'cof');
    	v=['OB Field (length ' num2str(length) ' Bytes)'];
    case {'SQ'}
    	%v=char(fread(fid,length,'uchar')');
    	v=fread(fid,length,'uchar')';
    	%fseek(fid, length,'cof');
    	%v=['SQ Field (length ' num2str(length) ' Bytes)'];
    otherwise,
      fseek(fid, length,'cof');
      v=-1;
  end

% ================================================================================  
function v=get_ahead_value(ahead, field, out_flag)
% ================================================================================  
	%disp(field);
	nposit=findstr(ahead,field);
	v=[];
	if ~isempty(nposit)
	   for p=1:length(nposit)
		npos=nposit(p);
		lpos=findstr(ahead(npos:end),10);
		lpos=lpos(1);
		tmp_str=ahead(npos:npos+lpos-2);
%		v=str2num(tmp_str(findstr(tmp_str,'=')+1:end));
		tmp_str=tmp_str(findstr(tmp_str,'=')+2:end);
		val=str2num(tmp_str);
		if isempty(val), v=tmp_str; else v(p)=val; end
		%if isempty(val), v{p}=tmp_str; else v{p}=val; end
		% switch field
		% case 'sSliceArray.lSize'
		% whos val v tmp_str
		% end
	   end
	end
       if (out_flag > 0), fprintf('ahead: %s = %s\n',field,v); end
	
