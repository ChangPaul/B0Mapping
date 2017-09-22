function [unwph,iter,erglist] = puma_ho(psi,p,varargin)
%puma_ho   Graph optimization based phase unwrapping algorithm.
%   [unwph,iter,erglist] = puma_ho(psi,p,'potential',potential,'cliques',cliques, 'qualitymaps', qualitymaps,
%   'schedule',schedule)
%   Unwrapps the observed modulo-2pi phase, casting the problem as an energy minimization via graph mincut
%   calculation. The algorithm is described in "Phase Unwrapping via Graph Cuts" submitted to IEEE IP, October, 2005.
%   Herein we make a generalization, by allowing cliques of any order (not just 1st order).
%
%   Authors: Jose Bioucas-Dias and Gon�alo Valad�o
%
%   Last change: Goncalo Valadao (19/9/2012 23h46m)
%
% ======================== REQUIRED INPUT PARAMETERS ======================
% Parameter             Values
% name                  and description
% =========================================================================
%
% psi                   (double) The wrapped phase image.
% p                     (double) It defines the clique potential exponent (>0).
%
% ======================== OPTIONAL INPUT PARAMETERS ======================
% Parameter             Values
% name                  and description
% =========================================================================
%
% potential             (1x1 struct array) This struct array has 2 fields:
% potential.quantized   ('yes','no') Default: 'no'.
% potential.threshold   (double) it defines a region over which the
%                        potential grows quadratically. By default is pi.
%
% cliques               (nx2 double matrix) Each row defines the
%                       "displacement vector" corresponding to each clique.
%                       The first and the second columns correspond to
%                       cliques along rows and columns in the image, respecively.
%                       By default is [1 0;0 1] (first order cliques).
%
% qualitymaps           (size(psi) x n (nocliques) double array). The quality matrices
%                       may take values between 0 and 1 (value 1: discontinuity presence;
%                       value 0: discontinuity absence).
%                       There is one quality matrix per clique type. By default there is
%                       discontinuity absence. A quality value corresponding to a certain
%                       clique must be signalled in the pixel corresponding to the end of
%                       the clique displacement vector (for each pair of pixels).
%
% schedule              (double vector) This vector contains a schedule of jump sizes.
%
% verbose               ('yes', 'no') -> display the unwrapped  phase along  
%                        the iterations.
%                        Default = 'yes'
%
%
% ========================= OUTPUT PARAMETERS =============================
% Parameter             Values
% name                  and description
% =========================================================================
% unwph                 (double array) This is the unwrapped phase image.
% iter                  (double) This is the number of iterations the algorithm runs through
% erglist               (double vector) This is a sequence of the energies of the unwrapped
%                       phases along the algorithm run.
%
% =========================== EXAMPLES ====================================
%       Note: the optional arguments must be provided in pairs string+value;
%       those pairs may be entered in any order.
%
%       potential.quantized = 'no'; potential.threshold = 0.5;
%       [unwph,iter,erglist] = puma_ho(psi,2,'potential',potential)
%
%       potential.quantized = 'yes'; potential.threshold = 2;
%       cliques = [1 0; 0 1; -1 1];
%       [unwph,iter,erglist] = puma_ho(psi,1,'cliques',cliques,'potential',potential)
%
%       potential.quantized = 'no';
%       potential.threshold = 0.1; cliques = [1 1];
%       qualitymaps = ones(size(psi,1),size(psi,2))
%       [unwph,iter,erglist] = puma_ho(psi,p,'potential',potential,'cliques',cliques,'qualitymaps',qualitymaps)

% ========================== REFERENCES ===================================
%   For reference see:
%   J. Bioucas-Dias and G. Valad�o, "Phase Unwrapping via Graph Cuts"
%   IEEE Transactions Image Processing, 2007 (to appear).
%   The algorithm here coded corresponds to a generalization for any
%   cliques set (not only vertical and horizontal).
%
%   J. Bioucas-Dias and J. Leit�o, "The ZpiM Algorithm for Interferometric Image Reconstruction
%   in SAR/SAS", IEEE Transactions Image Processing, vol. 20, no. Y, 2001.
%
%   The technique here employed is also based on the article:
%   V. kolmogorov and R. Zabih, "What Energy Functions can be Minimized via Graph Cuts?",
%   European Conference on Computer Vision, May 2002.
% =========================================================================
%
%
%
% Modification: 
%     
%    1 - Fix a bug in the way discontinties were deal with 
%        (Gon�alo Valad�o, Sep.,2012)
%
%    2 - Change in the potential default parameters:
%        potential.quantized = 'no';
%        potential.threshold = pi;
%        (J Bioucas-Dias, Sep.,2012)
%
%    3 - Introdution of the verbose input parameter:
%        verbose = 'yes' -> display (iter, enerrg_actual, jump_size)
%                           and display the unwrapped  phase along
%                           the iterations.
%        Default = 'yes'
%        
%    
%   
%
% ------------------------------------------------------------------
% Author: Gon�alo Valad�o & Jose Bioucas-Dias, 2007
%
%

%
%% -------------------------------------------------------------------------
%
% Copyright (July, 2007):        Gon�alo Valad�o (gvaladao@lx.it.pt) 
%                                Jos� Bioucas-Dias (bioucas@lx.it.pt)
%
% PUMA_HO is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Default values
potential.quantized     = 'no';
potential.threshold     = pi;
cliques                 = [1 0;0 1];
qualitymaps             = repmat(zeros(size(psi,1),size(psi,2)),[1,1,2]); qual=0;
schedule                = 1;
verbose                 = 'no';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test for number of required parameters

% Error out if there are not at least the two required input arguments

if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end

% Read the optional parameters
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
elseif length(varargin)~=0
    for i=1:2:(length(varargin)-1)
        % change the value of parameter
        switch varargin{i}
            case 'potential'        % potential definition
                potential = varargin{i+1};
            case 'cliques'          % cliques to consider
                cliques = varargin{i+1};
            case 'qualitymaps'      % the quality maps
                qualitymaps = varargin{i+1};
                qual = 1;
            case 'schedule'         % jump size schedule
                schedule = varargin{i+1};
           case 'verbose'           % display partial unwrappings
                verbose = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end;
    end
end;

if (qual==1)&&(size(qualitymaps,3)~=size(cliques,1))
    error('qualitymaps must be a 3D matrix whos 3D size is equal to no. cliques. Each plane on qualitymaps corresponds to a clique.');
end

% INPUT AND INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%
th    = getfield(potential,'threshold');
quant = getfield(potential,'quantized');


[m,n] = size(psi); % Size of input
kappa = zeros(m,n); % Initial labeling
%kappa = round(rand(m,n)*40);
kappa_aux = kappa;
iter = 0;
erglist = [];

[cliquesm,cliquesn] = size(cliques); % Size of input cliques
if qual ==0
    qualitymaps             = repmat(zeros(size(psi,1),size(psi,2)),[1,1,size(cliques,1)]);
end
disc_bar = 1 - qualitymaps;
% "maxdesl" is the maximum clique length used.
maxdesl = max(max(abs(cliques)));
% We define "base" which is a mask having ones in the region of interest(psi) and zeros upon a passe-partout
% having a constant length maxdesl+1.
base = zeros(2*maxdesl+2+m,2*maxdesl+2+n); base(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = ones(m,n);


      %% Temporary movie maker
%        mov = avifile('example5.avi');
        %% Temporary movie maker

% PROCESSING   %%%%%%%%%%%%%%%%%%%%%%%%
for jump_size = schedule

  
    
    
    possible_improvment = 1;
    erg_previous = energy_ho(kappa,psi,base,p,cliques,disc_bar,th,quant);

    while possible_improvment
        iter = iter + 1;
        erglist = [erglist erg_previous];
        remain = [];
        % Here we put a passe-partout (constant length = maxdesl+1) in the images kappa and psi
        base_kappa = zeros(2*maxdesl+2+m,2*maxdesl+2+n); base_kappa(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = kappa;
        psi_base = zeros(2*maxdesl+2+m,2*maxdesl+2+n); psi_base(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = psi;

        
        %%%%% Added by Goncalo Valadao on 19/11/2012 22:15 %%%%%
        z = size(disc_bar,3);
        base_disc_bar  = repmat(zeros(2*maxdesl+2+m,2*maxdesl+2+n),[1 1 z]); base_disc_bar(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,:) = disc_bar;
        %%%%% Added by Goncalo Valadao on 19/11/2012 22:15 %%%%%
        
        for t = 1:cliquesm
            % The allowed start and end pixels of the "interpixel" directed edge
            base_start(:,:,t) = circshift(base,[-cliques(t,1),-cliques(t,2)]).*base;
            base_end(:,:,t) = circshift(base,[cliques(t,1),cliques(t,2)]).*base;

            % By convention the difference images have the same size as the
            % original ones; the difference information is retrieved in the
            % pixel of the image that is subtracted (end of the diff vector)
            auxili = circshift(base_kappa,[cliques(t,1),cliques(t,2)]);
            t_dkappa(:,:,t) = (base_kappa-auxili);
            auxili2 = circshift(psi_base,[cliques(t,1),cliques(t,2)]);
            dpsi = auxili2 - psi_base;
            
            % Beyond base, we must multiply by
            % circshift(base,[cliques(t,1),cliques(t,2)]) in order to
            % account for frontier pixels that can't have links outside ROI
            
            %%%%%%%%  Changed on 19/11/2012 18:30 by Goncalo Valadao  %%%%%%%%

%             a(:,:,t) = (2*pi*t_dkappa(:,:,t)-dpsi).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);
%             A(:,:,t) = clique_energy_ho(abs(a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);
%             D(:,:,t) = A(:,:,t);
%             C(:,:,t) = clique_energy_ho(abs(2*pi*jump_size + a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);
%             B(:,:,t) = clique_energy_ho(abs(-2*pi*jump_size + a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);
            
            a(:,:,t) = (2*pi*t_dkappa(:,:,t)-dpsi).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);
            A(:,:,t) = clique_energy_ho(abs(a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t);
            D(:,:,t) = A(:,:,t);
            C(:,:,t) = clique_energy_ho(abs(2*pi*jump_size + a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t);
            B(:,:,t) = clique_energy_ho(abs(-2*pi*jump_size + a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t);
            
            %%%%%%%%  End of changes on 19/11/2012 18:30 by Goncalo Valadao  %%%%%%%%  
            
            
            
            
            
            % The circshift by [-cliques(t,1),-cliques(t,2)] is due to the fact that differences are retrieved in the
            % "second=end" pixel. Both "start" and "end" pixels can have source and sink connections.
            source(:,:,t) = circshift((C(:,:,t)-A(:,:,t)).*((C(:,:,t)-A(:,:,t))>0),[-cliques(t,1),-cliques(t,2)]).*base_start(:,:,t);
            sink(:,:,t) = circshift((A(:,:,t)-C(:,:,t)).*((A(:,:,t)-C(:,:,t))>0),[-cliques(t,1),-cliques(t,2)]).*base_start(:,:,t);

            source(:,:,t) = source(:,:,t) + ((D(:,:,t)-C(:,:,t)).*((D(:,:,t)-C(:,:,t))>0)).*base_end(:,:,t);
            sink(:,:,t) = sink(:,:,t) + ((C(:,:,t)-D(:,:,t)).*((C(:,:,t)-D(:,:,t))>0)).*base_end(:,:,t);
        end


        % We get rid of the "pass-partous"
        source(1:maxdesl+1,:,:)=[]; source(m+1:m+maxdesl+1,:,:)=[]; source(:,1:maxdesl+1,:)=[]; source(:,n+1:n+maxdesl+1,:)=[];
        sink(1:maxdesl+1,:,:)=[]; sink(m+1:m+maxdesl+1,:,:)=[];sink(:,1:maxdesl+1,:)=[]; sink(:,n+1:n+maxdesl+1,:)=[];
        auxiliar1 = B + C - A - D;
        auxiliar1(1:maxdesl+1,:,:)=[]; auxiliar1(m+1:m+maxdesl+1,:,:)=[]; auxiliar1(:,1:maxdesl+1,:)=[]; auxiliar1(:,n+1:n+maxdesl+1,:)=[];
        base_start(1:maxdesl+1,:,:)=[]; base_start(m+1:m+maxdesl+1,:,:)=[]; base_start(:,1:maxdesl+1,:)=[]; base_start(:,n+1:n+maxdesl+1,:)=[];
        base_end(1:maxdesl+1,:,:)=[]; base_end(m+1:m+maxdesl+1,:,:)=[]; base_end(:,1:maxdesl+1,:)=[]; base_end(:,n+1:n+maxdesl+1,:)=[];

        % We construct the "remain" and the "sourcesink" matrices
        for t=1:cliquesm
            start = find(base_start(:,:,t)~=0); endd = find(base_end(:,:,t)~=0);
            auxiliar2 = auxiliar1(:,:,t);
            auxiliar3 = [start endd  auxiliar2(endd).*(auxiliar2(endd)>0) zeros(size(endd,1),1)];
            %auxiliar3 = [start endd  auxiliar2(endd) zeros(size(endd,1),1)];
            remain = [remain; auxiliar3];
        end
        %remain = sortrows(remain,[1 2]);
        sourcefinal = sum(source,3);
        sinkfinal = sum(sink,3);
        sourcesink = [(1:m*n)' sourcefinal(:) sinkfinal(:)];

        % KAPPA RELABELING
        [flow,cutside] = mincut(sourcesink,remain);
        % The relabeling of each pixel relies on the cutside of that pixel:
        %       if the cutside = 0 (source) then we increment the label by jump_size
        %       if the cutside = 1 (sink) then the label remains unchanged
        kappa_aux(cutside(:,1)) = kappa(cutside(:,1)) + (1 - cutside(:,2))*jump_size;

        % CHECK ENERGY IMPROVEMENT
        erg_actual = energy_ho(kappa_aux,psi,base,p,cliques,disc_bar,th,quant);
        %
        if (erg_actual < erg_previous)
            erg_previous = erg_actual;
            kappa = kappa_aux;
        else
            possible_improvment = 0;
            unwph = 2*pi*kappa + psi;
        end
%         mesh( 2*pi*kappa + psi)
%         view(-30,30);
%         surfl(2*pi*kappa + psi);shading interp; colormap(gray);
        
        if strcmp(verbose,'yes')
           imagesc(2*pi*kappa + psi); colormap(gray);
           % display to the current figure
           figure(gcf);
           drawnow;
        end
           

        %surfl(2*pi*kappa + psi);shading interp; 
        %colormap(gray);
        %imagesc( 2*pi*kappa + psi);
        

        
        
        %% Temporary movie maker
%        F=getframe(gca);
%        mov=addframe(mov,F);
        %% Temporary movie maker        
        
        %drawnow;
        %iter %#ok<NOPRT>
        %erg_actual %#ok<NOPRT>
        %jump_size %#ok<NOPRT>
        clear base_start base_end source sink auxiliar1 auxiliar2 A B C D;
    end % while
    if strcmp(verbose,'yes')
        title('Puma solution');
    end    
end %for

% Temporary movie maker
%mov=close(mov);
% Temporary movie maker
end

function e=clique_energy_ho(d,p,th,quant)
%clique_energy_ho  Computes clique energy: e=th^(p-2)*d.^2.*mask + d.^p.*(1-mask)
%        e=clique_energy_ho(d,p,th,quant)
%
%  Input arguments --------------------
%  d             -> clique difference
%  p             -> power law exponent
%  th            -> it defines a region over which the potential grows quadratically
%  quant         -> it defines whether or not the potential is quantized

switch quant
    case 'no'        % non quantized potential
        d=abs(d);
    case 'yes'       % quantized potential (2pi Quantization of phase difference)
        d=abs(round(d/2/pi)*2*pi);
end;
%

if th~=0
    mask = (d<=th);
    e = th^(p-2)*d.^2.*mask + d.^p.*(1-mask);
else
    e = d.^p;
end
return
end

function erg = energy_ho(kappa,psi,base,p,cliques,disc_bar,th,quant)
%energy_ho   Energy from kappa labeling and psi phase measurements.
%   erg = energy_ho(kappa,psi,base,p,cliques,disc_bar,p,th,quant) returns the energy of kappa labeling given the 
%   psi measurements image, the base ROI image (having ones in the region of interest (psi) and a passe-partout
%   made of zeros), the exponent p, the cliques matrix (each row indicating a displacement vector corresponding
%   to each clique), the disc_bar (complement to one of the quality maps), a threshold th defining a region for
%   which the potential (before a possible quantization) is quadratic, and quant which is a flag defining whether
%   the potential is or is not quantized.
%   (see J. Bioucas-Dias and G. Valad�o, "Phase Unwrapping via Graph Cuts"
%   submitted to IEEE Transactions Image Processing, October, 2005).
%   SITE: www.lx.it.pt/~bioucas/ 

[m,n] = size(psi);
[cliquesm,cliquesn] = size(cliques); % Size of input cliques
maxdesl = max(max(abs(cliques))); % This is the maximum clique length used
% Here we put a passe-partout (constant length = maxdesl+1) in the images kappa and psi
base_kappa    = zeros(2*maxdesl+2+m,2*maxdesl+2+n); base_kappa(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = kappa;
psi_base      = zeros(2*maxdesl+2+m,2*maxdesl+2+n); psi_base(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = psi;
z = size(disc_bar,3);
base_disc_bar  = repmat(zeros(2*maxdesl+2+m,2*maxdesl+2+n),[1 1 z]); base_disc_bar(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,:) = disc_bar;

for t = 1:cliquesm
    % The allowed start and end pixels of the "interpixel" directed edge
    base_start(:,:,t) = circshift(base,[-cliques(t,1),-cliques(t,2)]).*base;
    base_end(:,:,t) = circshift(base,[cliques(t,1),cliques(t,2)]).*base;
    
    % By convention the difference images have the same size as the
    % original ones; the difference information is retrieved in the
    % pixel of the image that is subtracted (end of the diff vector)
    auxili = circshift(base_kappa,[cliques(t,1),cliques(t,2)]);
    t_dkappa(:,:,t) = (base_kappa-auxili);
    auxili2 = circshift(psi_base,[cliques(t,1),cliques(t,2)]);
    dpsi = auxili2 - psi_base;
    % Beyond base, we must multiply by
    % circshift(base,[cliques(t,1),cliques(t,2)]) in order to
    % account for frontier pixels that can't have links outside ROI
    a(:,:,t) = (2*pi*t_dkappa(:,:,t)-dpsi).*base.*circshift(base,[cliques(t,1),cliques(t,2)])...
               .*base_disc_bar(:,:,t);
end

erg = sum(sum(sum((clique_energy_ho(a,p,th,quant)))));

end

function [flow,cutside] = mincut(sourcesink,remain)
%   MINCUT   Graph Max-Flow/Min-Cut.
%   [FLOW,CUTSIDE] = MINCUT(SOURCESINK,REMAIN) returns a graph's maximum flow and minimum cut.
%   The graph is defined by two matrices given as the two function arguments:
%   
%   SOURCESINK: nx3 matrix where the first column has a node number and the other
%   two the correspondent weights of the arcs linking that node to the source and
%   sink respectively. There is one line for each node (n = total amount of nodes).
%   
%   REMAIN: this is a mx4 matrix containing the information about the links between
%   nodes and respective weights. Column 1 goes for starting node, column 2 for final 
%   node and the final two columns have the direct and inverse link weights respectively.     
%
%   The outputs are given by:
%   
%   FLOW: a scalar.
%   
%   CUTSIDE: nx2 matrix where the first column has a node number and the second column
%   has one of the following values: 
%
%   0 - the node stays at the source side of the CUT.
%   1 - the node stays at the sink side of the CUT.


%   This m-file calls a mex file that contains a software library that is a modification 
%   of the maxflow algorithm described in:
% 
% 	An Experimental Comparison of Min-Cut/Max-Flow Algorithms
% 	for Energy Minimization in Computer Vision.
% 	Yuri Boykov and Vladimir Kolmogorov.
% 	In Third International Workshop on Energy Minimization
% 	Methods in Computer Vision and Pattern Recognition, September 2001

sourcesinkg = single(sourcesink);
remaing = single(remain);

[flowg,cutsideg] = mf2(sourcesinkg,remaing);

flow = double(flowg); 
cutside = double(cutsideg);

% End
end
