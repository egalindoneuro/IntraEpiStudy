function CSDTime = CSD_2D(dirname,newFs)
% INPUT:  directory and name WITHOUT channel number in the form:
%         '/Volumes/F/Data/LxS02_spont_EU96L-M32x1_01_main/RAW_1_EU96L_0';
%     dirname     =  '/Volumes/F/Data/LxS02_spont_EU96L-M32x1_01_main/RAW_1_EU96L_0';
%
% OUTPUT: 3D array, X position, Y position and time. For the analysis I use
%         an interpolation and that is why there are more elements in the
%         grid.
%
% e.g.

filtorder   =   4;
%newFs       =   500;
%
xc   = [  0     1.2    2.4    3.6    4.8    6.0    7.2    8.4    9.6   0.6    1.8   3.0   4.2   5.4   6.6  7.8 ...
        9.0    10.2      0    1.2    2.4    3.6    4.8    6.0    7.2   8.4    9.6  10.8   0.6   1.8   3.0  ...
        4.2     5.4    6.6    7.8    9.0   10.2   11.4    1.2    2.4   3.6    4.8   6.0   7.2   8.4   9.6  10.8 ...
        12.0    0.6    1.8    3.0    4.2    5.4    6.6    7.8    9.0  10.2   11.4     0   1.2   2.4   3.6  ...
        4.8     6.0    7.2    8.4    9.6   10.8   12.0   -0.6    0.6   1.8    3.0   4.2   5.4   6.6   7.8  ...
        9.0    10.2   11.4    2.4    3.6    4.8    6.0    7.2    8.4   9.6   10.8   1.8   3.0   4.2   5.4  ...
        6.6     7.8    9.0   10.2];

yc   = [   0      0      0      0      0      0      0      0      0  1.04   1.04   1.04  1.04  1.04  1.04  1.04 ...
        1.04   1.04   2.08   2.08   2.08   2.08   2.08   2.08   2.08  2.08   2.08   2.08  3.12  3.12  3.12  ...
        3.12   3.12   3.12   3.12   3.12   3.12   3.12   4.16   4.16  4.16   4.16   4.16  4.16  4.16  4.16  4.16 ...
        4.16    5.2    5.2    5.2    5.2    5.2    5.2    5.2    5.2   5.2    5.2   6.34  6.34  6.34  6.34  ...
        6.34   6.34   6.34   6.34   6.34   6.34   6.34   7.48   7.48  7.48   7.48   7.48  7.48  7.48  7.48  ...
        7.48   7.48   7.48   8.52   8.52   8.52   8.52   8.52   8.52  8.52   8.52   9.66  9.66  9.66  9.66  ...
        9.66   9.66   9.66   9.66   ];
%
uniX  =  unique(xc);
uniY  =  unique(yc);
%
% remapping into square tiling
xy     =  0;
for i1 = 1 : numel(uniX)
    for j1 = 1 : numel(uniY)
        chnltmp = find(xc == uniX(i1) & yc == uniY(j1));
        if ~isempty(chnltmp)
            xy             = xy + 1;
            chn2in(xy)     = chnltmp;
            if chnltmp < 10
                chnl2 = ['0',num2str(chnltmp)];
            else
                chnl2 = num2str(chnltmp);
            end
            %
            namefile          =  ([dirname chnl2]);
            namesig           =  [namefile, '.bin'];
            namehdr           =  [namefile, '.mat'];
            load(namehdr)
%            NyFre             =  hdr.fs/2;
            NyFre             =  newFs/2;
            % filter Hz  ==============================
            [b1,a1]           =  butter(filtorder,[1 4]/NyFre);
            fid               =  fopen(namesig);
            sigCh             =  fread(fid,'int16=>double');
            fclose(fid);       
            [n,d]             =  rat(hdr.fs/newFs);
            sigCh2tmp         =  zscore(resample(sigCh,d,n));
            %
            if i1 == 1 & j1 == 1
                VFiel  =  NaN(numel(uniY),numel(uniX),numel(sigCh2tmp));
                VShort =  NaN(numel(uniY),numel(sigCh2tmp));
            end
            %
            sigCh2            =  filtfilt(b1,a1,sigCh2tmp);
            VFiel(j1,i1,:)    =  sigCh2;
            VShort(chnltmp,:) =  sigCh2;
%        else
        end
    end
end
%
% filling empty tiles
[XX,YY] = meshgrid(uniX,uniY);
clear VF
smpLFP  = NaN(10,22,numel(sigCh2tmp),1);
VFShort = NaN(11,13);
%
for ii = 1 : numel(sigCh2tmp)
    VF               =  squeeze(VFiel(:,:,ii));
%   VF               =  flipud(VF);
%   VF(end+1,end+1)  =  NaN;
    for r2 = 2 : 2 : size(VF,1)
        VFShort(r2,2:11) = VF(r2,2:2:end-2) + (diff(VF(r2,2:2:end))/2);
    end
    for r2 = 1 : 2 : size(VF,1)
        VFShort(r2,1:11) = VF(r2,1:2:end-1);
    end
     sig = squeeze(VShort(:,ii));
% %    sig = flipud(sig);
     zi  = griddata(xc',yc',sig,XX,YY);
     zi  = flipud(zi);
     smpLFP(:,:,ii,1) = zi;
end
%% CSD
Laplc = NaN(10,22,numel(sigCh2tmp));
for t3 = 1 : numel(sigCh2tmp)
    % ================================
    %
    data1          =  smpLFP(:,:,t3);
    [grax,gray]    =  gradient(data1);
    [graxx,graxy]  =  gradient(grax);
    [grayx,grayy]  =  gradient(gray);
    %
    % ================================
    %
    Laplc(:,:,t3) = graxx + grayy;    
end
%
CSDTime   = Laplc;
%CSDTime   = zscore(Laplc,[],3);

% figure
% for t4 = 1 : 1000   
%     pcolor(squeeze(Laplc3(:,:,t4)));
%     set(gca,'clim',[-8 8])
%     colormap('jet')
%     pause(0.001)
% end
% 
% LaplcGT = Laplc3;
% LaplcLT = Laplc3;
% %
% % detect sink/sources in ECoG (3 std)
% LaplcGT(Laplc3<3)   = 0;
% LaplcGT(Laplc3>=3)  = 1;
% %
% LaplcLT(Laplc3>-3)  = 0;
% LaplcLT(Laplc3<=-3) = -1;
% %
% % finding coordinates and time for sources and sinks
% XY = [];
% T3 = [];
% for x1 = 1 : size(LaplcGT,1)
%     [y2,t3]                   = find(squeeze(LaplcGT(x1,:,:))==1);    
%     XY(:,end+1:end+numel(y2)) = [x1*ones(numel(y2),1) y2]';
%     T3(end+1:end+numel(y2))   = t3;    
% end
% [TAll,indxOrd] = sort(T3);
% XYCoor = XY(:,indxOrd);
% 
% uniX(4)
% uniY(4)
% %
% 
% 
% 
