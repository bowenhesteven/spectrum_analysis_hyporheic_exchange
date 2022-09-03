clear all;
%cd '/scratch/heb1'; % cd to the image directory;
I = imread('hanlidar_nad88.tif'); % read the index image which is the first image of the file;
imshow(I);

% [J,rect]=imcrop(I); % Crop the image;
xmin=1.43e+04;
height=299;
ymin=10200+height;
width=999;
% RGB = insertMarker(I,[xmin ymin-height;xmin ymin;xmin+width ymin-height;xmin+width ymin]);
rect=[xmin ymin width height];
J = imcrop(I,rect);
topo_piece=J; % The elevation of the cropped area;
[nr,nc]=size(topo_piece); % Get the x coordinate and y coordinate size of the cropped area of image;
x = 1:nc;
y = 1:nr; y = fliplr(y);
[X,Y] = meshgrid(x,y);
topopiece_vector = struct('x',X(1,:),'y',Y(:,1),'z',topo_piece); % Reconstruct the index of the cropped area of the image;
%% Set up the Matlab desktop

% First, take a moment to dock the Editor so you can see this script and
% the command line at the same time.

% Make new figures docked and white
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureColor','w')

%% Import gridded data

% The LiDAR data are stored as an ArcGIS-formatted ASCII grid (*.asc)
%[Z, dim] = ReadArcGrid('merced');

% The struct variable dim contains fields describing the coordinates and
% dimensions of Z. These are based on the information extracted from the
% grid file header.
dx = abs(topopiece_vector.x(2) - topopiece_vector.x(1)); % grid spacing in the x-direction
dy = abs(topopiece_vector.y(2) - topopiece_vector.y(1)); % grid spacing in the y-direction
[Ny,Nx] = size(topopiece_vector.z); % grid dimensions
%% View data

% Display a shaded relief map (Fig. 2)
figure('Name','Fig. 2: Shaded relief','NumberTitle','off')
ShadePlot(topopiece_vector.x,topopiece_vector.y,topopiece_vector.z)

% Use the magnifying glass tool to zoom in for a closer look at the mounds.
% What is the typical inter-mound spacing? Note that the grid spacing is 
% 1m. They're about 10 pixels apart, so we'll be on the lookout for a
% feature at a wavelength of ~10m in the power spectrum.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. CALCULATE POWER SPECTRUM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing

% Before calculating the power spectrum, we'll make some decisions about
% the preprocessing steps we'll take. If the DEM has a non-zero mean or a
% background slope across the entire grid, this will contaminate the
% spectrum with long-wavelength signals. So our first step will be to 
% detrend the DEM by fitting a least-squares plane to the elevations and 
% then subtracting this fit.
Zo = topopiece_vector.z; % Save the original elevations for later
Z = Detrend(topopiece_vector.z);
plane = Zo - Z; % Save the least-squares plane for re-trending later

% Second, the fast Fourier transform proceeds fastest if the dimensions of 
% the input matrix are integer powers of two, which we can achieve by 
% padding the DEM with zeros. 
pad = 1; % 1 means pad the data, 0 no padding.

% Third, because the edges of our DEM are not perfectly periodic, the
% spectrum can become contaminated by frequencies used to "fit" the edge
% discontinuity. We can mitigate this effect by multiplying the DEM by a
% function that tapers to zero at the edges.
window = 1; % 1 means window the data, 0 no window
%% 2D FFT

% Calculate the power spectrum using the 2D Fast Fourier Transform (FFT).
% Note that the zero padding and windowing happen inside fft2D.
[Pm, fm, Pv, fv] = fft2D(Z,dx,dy,pad,window); 

% A few words about the output from the previous step:
% 
% We have calculated the Discrete Fourier Transform (DFT) periodogram. This
% is different from the power spectral density, which is also commonly used
% as an estimate of the power spectrum. 
%
% Pm is the 2D power spectrum matrix, and fm is the corresponding matrix of 
% radial frequencies. Pv is a vector version of the power spectrum, and fv 
% is the corresponding vector of frequencies. Units of Pm and Pv are 
% [units of Z]^2, which in our case is length^2 since Z is a matrix of 
% elevations. Units of fm and fv are [units of x]^-1, which for us is 
% 1/length since x is distance. Figs. 3a-e show examples of a simple 1D
% signal and its power spectrum, and Fig. 3f shows an example of a simple
% surface and its 2D power spectrum.
%
% What do the units mean? The periodogram is a measure of how much of the
% original elevation field's variance falls within a given frequency range.
% You can check that the sum of the periodogram is roughly equal to the
% variance in Z. (It will be somewhat less due to the zero padding.)
%
% What about the frequency limits? Wavelength is equal to the inverse of 
% frequency. The shortest wavelength that can be resolved by 2D data is 
% oriented diagonally, and is equal to sqrt(dx^2 + dy^2), or sqrt(2) for 
% our data (Fig. 4). The longest wavelength that can be resolved is 
% infinite in principle, but in practice we wouldn't feel comfortable 
% trying to detect any signal with a wavelength longer than our dataset. To 
% be even more conservative, the longest wavelength we trust should be a 
% few times shorter than the width of the DEM -- a few hundred meters in 
% this case.

% Clean up workspace
clear pad window
%% %%%%%%%%%%%%%%%%%%%%%%%%
% 3. ANALYZE THE SPECTRUM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D Power Spectrum

% Plot the 1D version of the spectrum (Fig. 5). We'll plot the 2D spectrum 
% a few steps later, but for now it's easier to visualize in 1D.
s1d = figure('Name','Fig. 5: 1D Spectrum','NumberTitle','off');
% subplot(2,1,1)
[axf,axw] = SpecPlot1D(fv,Pv);

% Do we see anything at a 10m wavelength? Yes, there is a broad peak. The 
% shorter-wavelength peaks are just pixel-scale noise.
hold(axw,'on')
plot([10 10],get(axw,'ylim'),'k')
%% Background Spectrum

% Now let's see what the spectrum looks like in 2 dimensions. If we just
% looked at it as is, we wouldn't see much: as the 1D spectrum shows,
% longer-wavelength signals have much higher power, and would swamp any
% shorter-wavelength detail. To make the 2D spectrum easier to view, we'll
% remove the power-law background trend. Because there are so many more
% points at higher frequencies, we'll bin the 1D spectrum and fit the trend
% to the binned values (Fig. 5).

nbin = 20; % Number of bins
B = bin(log10(fv),log10(Pv),nbin,0); % Bin the log-transformed data

% Plot the binned values
hold(axf,'on')
plot(axf,10.^B(:,1),10.^B(:,2),'ok','markerfacecolor','w')

% Fit a trend with the form P ~ 1/f^n, and plot it
fit = robustfit(B(:,1),B(:,2));
plot(axf,10.^B1(:,1),10^fit_Background(1)*(10.^B1(:,1)).^fit_Background(2),'b'); % 10^(A1+A2logx)=y;
plot(axf,10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k')
%% Background spectral (Artifial randomly rough surfaces generation)
% std = std2(Zo); % standard deviation of the original topography;
% Variance=var(Z(:)); % total variance of the original topography;
% Mean=mean(Zo(:)); % mean value of the original topography;
Variance=sum(Pm(:)); % the sum of the total Pdft array;
std=sqrt(Variance);
%[N,R]=boxcount(Z,'slope'); % also shows the semi-log plot of the local slope
% DF = - diff(log(N))/diff(log(R));
H=0.5; % Hurst exponent (H), which is related to the fractal dimension of a surface topography D=3-H. 
       % For example, a Brownian surface roughness has H= 0.5. 
       % The parameter H can take values between 0 and 1. 
size_Z=size(fm); % get size of the frequency matrix;
n=size_Z(1,1);
m=size_Z(1,2);
Lx=size_Z(1,2);
Pv_t=zeros(524800,1);
% fv_size=size(fv);
% fv_t=zeros(fv_size(1,1),1);
Pm_t=zeros(1024,1024);
loop=5;% std=zeros(loop,1)+std2(Zo);
% H=zeros(loop,1)+0.5;
% Lx=zeros(loop,1)+cols;
% m=zeros(loop,1)+cols;
% n=zeros(loop,1)+rows;
Pm_art=NaN(1024,1024);
Pv_art=NaN(524800,1);
z_art_mat = zeros(1024, 1024, loop);
for i=1:loop
[z_art, PixelWidth, PSD] = artificial_surf(std, H, Lx, m, n); 
z_art_mat(:,:,i)=z_art;
%Zo_BG = z_art_mat(:,:,i); % Save the original elevations for later
% Z_BG= Detrend(z_art_mat(:,:,i));
% plane = Zo_BG - Z_BG; % Save the least-squares plane for re-trending later

% Second, the fast Fourier transform proceeds fastest if the dimensions of 
% the input matrix are integer powers of two, which we can achieve by 
% padding the DEM with zeros. 
 pad = 0; % 1 means pad the data, 0 no padding.

% Third, because the edges of our DEM are not perfectly periodic, the
% spectrum can become contaminated by frequencies used to "fit" the edge
% discontinuity. We can mitigate this effect by multiplying the DEM by a
% function that tapers to zero at the edges.
 window = 0; % 1 means window the data, 0 no window
 
[Pm_BG, fm_BG, Pv_BG, fv_BG] = fft2D(z_art_mat(:,:,i),dx,dy,pad,window); 
 Pm_art(:,:,i)=Pm_BG;
 Pm_t=Pm_t+Pm_art(:,:,i);
 Pv_art(:,:,i)=Pv_BG;
 Pv_t=Pv_t+Pv_art(:,:,i);
end
Pm_average=Pm_t./loop;
Pv_average=Pv_t./loop;


% Plot the 1D version of the spectrum (Fig. 5). We'll plot the 2D spectrum 
% a few steps later, but for now it's easier to visualize in 1D.
s2d = figure('Name','Fig. 6: 1D Spectrum of background','NumberTitle','off');
%subplot(2,1,1)
[axf_BG, axw_BG] = SpecPlot1D(fv_BG,Pv_average);
% axf and axw are two verisons of axis;
% Do we see anything at a 10m wavelength? Yes, there is a broad peak. The 
% shorter-wavelength peaks are just pixel-scale noise.
hold(axw_BG,'on')
%plot([10 10],get(axw,'ylim'),'k') % plot straightline from (10,0) to (10,1)

nbin_BG = 20; % Number of bins
B1 = bin(log10(fv),log10(Pv_average),nbin_BG,0); % bin is a written function! Bin the log-transformed data
% output      is matrix consisting of nbin rows and 8 columns: 
%               [x (center of bin), mean, standard deviation, standard 
%               error, n, max, min, median]

% Plot the binned values
hold(axf_BG,'on')
plot(axf_BG,10.^B1(:,1),10.^B1(:,2),'ok','markerfacecolor','w')

% Fit a trend with the form P ~ 1/f^n, and plot it
fit_Background = robustfit(B1(:,1),B1(:,2)); % returns the 2 coefficient. A1+A2logx=logy;
plot(axf_BG,10.^B1(:,1),10^fit_Background(1)*(10.^B1(:,1)).^fit_Background(2),'k'); % 10^(A1+A2logx)=y;

%% Quantify the good of fitness
s3d = figure('Name','Fig. 7: Goodness of fit','NumberTitle','off');
r2=zeros(20,1); % create a vector to save the r2 value;
RSS=0; % set the initial value 0 to RSS;
y_t=0; % set the initial value 0 to the sum of original data point;
mean_y=0; % set the initial value 0 to the mean of original data;
TSS=0; % set the initial value 0 to TSS;
for i=1:20 % 20 data points;
RSS=RSS+(10.^B(i,2)-10^fit_Background(1)*(10.^B1(i,1)).^fit_Background(2))^2;
y_t=y_t+10.^B(i,2);
mean_y=mean(y_t);
TSS=TSS+((10.^B(i,2)-mean_y))^2;
r2(i,1)=1-RSS/TSS;
end
% [r2,r2adj]=rsquared(10.^B(:,2),10^fit_Background(1)*(10.^B1(:,1)).^fit_Background(2),2);
 PlotR2(10.^B(:,1),r2(:));
 
 %% Normalized spectra (normalized by generated surfaces)
Pv_Nor=Pv./Pv_average;
s4d = figure('Name','Fig. 8: Normalized Power','NumberTitle','off');
%subplot(2,1,1)fv-fv_BG
[axf_Nor , axw_Nor] = NormaSpecPlot1D(fv,Pv_Nor);
% Use this fit to normalize the 2D spectrum
% Pmn = Pm./(10^fit(1)*fm.^fit(2)); % Pm is the DFT periodogram; fm is the frequency matrix;
%% Normalized spectra (normalized by generated surfaces)
Pv_Nor=Pv./Pv_average;
s4d = figure('Name','Fig. 8: Normalized Power','NumberTitle','off');
%subplot(2,1,1)fv-fv_BG
[axf_Nor , axw_Nor] = NormaSpecPlot1D(fv,Pv_Nor);
 %% 2D Power Spectrum

% Use this fit to normalize the 2D spectrum
 %Pmn = Pm./(10^fit(1)*fm.^fit(2)); % Pm is the DFT periodogram; fm is the frequency matrix;
Pm_average(Pm_average<10^-15)=NaN;
Pmn = Pm./Pm_average;
Pmn(Pmn<0.5*th)=0;
% (10^fit(1)*fm.^fit(2)) is the result of DFT periodogram calculated by
% fit;
% Plot the normalized 2D spectrum (Fig. 6)
figure('Name','Fig. 6: 2D Spectrum','NumberTitle','off')
SpecPlot2D(fm,log10(Pmn));
% The donut of elevated values at radial frequencies of ~0.1 is the
% spectral signature of the mima mounds. It's a donut because the mounds
% don't have a strongly preferred orientation over the entire landscape. 
% The two peaks very close to the center of the square correspond to the 
% NW-SE-trending valleys, and the sharp features aligned with these peaks 
% at higher frequencies are probably harmonics.

% Clean up workspace
clear Pm fm Pv nbin B axf axw fit Pmn