clear all;
%% Read the iamge and cut othe piece and plot it
I = imread('hanlidar_nad88.tif'); % read the index image which is the first image of the file;
% imshow(I); % show the image;
imagesc(I);
colormap jet;% Thanks to Bro BAO;
hold on

% The coordinate of the rectangle of the image;
xmin=2.549E+04; % The min x-direction coordinate;
height=5999; % The height of the section;
ymin=1; % The min y-direction coordinate;
width=11899; % The width of the section;

% cut the piece from the image;
topopiece_vector=cut_piece(I,xmin,ymin,width,height);

% Plot the piece of the topo, marker the coodinates of the rectangle/the cutted section;
plot(xmin,ymin+height,'r+', 'MarkerSize', 30)
hold on
plot(xmin,ymin,'r+', 'MarkerSize', 30)
hold on
plot(xmin+width,ymin,'r+', 'MarkerSize', 30)
hold on
plot(xmin+width,ymin+height,'r+', 'MarkerSize', 30)
%% Set up the Matlab desktop
% Make new figures docked and white
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureColor','w')
%% Import gridded data
dx = abs(topopiece_vector.x(2) - topopiece_vector.x(1)); % grid spacing in the x-direction
dy = abs(topopiece_vector.y(2) - topopiece_vector.y(1)); % grid spacing in the y-direction
[Ny,Nx] = size(topopiece_vector.z); % grid dimensions
%% View data
% Display a shaded relief map (Fig. 2)
figure('Name','Fig. 2: Shaded relief','NumberTitle','off')
ShadePlot(topopiece_vector.x,topopiece_vector.y,topopiece_vector.z) % ShadePlot of the topo/surface;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. CALCULATE POWER SPECTRUM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing
Zo = topopiece_vector.z; % Save the original elevations for later;
Z = Detrend(topopiece_vector.z); % Detrend the data;
Size_Section=size(Z); % get the size of the section;
width=Size_Section(1,2); % get the size of the width of the section;
length=Size_Section(1,1); % get the size of the length of the section;
plane = Zo - Z; % Save the least-squares plane for re-trending later;
pad = 0; % 1 means pad the data, 0 no padding.
window = 1; % 1 means window the data, 0 no window
%% 2D FFT
[P, Pm, fm, Pv, fv] = fft2D(Z,dx,dy,pad,window); % Apply the 2D fft;
% P is the 2D fft coefficient;
% Pm is the 2D fft power spectrum matrix (periodogram);
% fm is the 2D radial frequencies matrix;
% Pv is the 1D power spetrum vector;
% fv is the 1D radial frequencies vector;
[M , Wss] = Hann2D(Z); % Output the detrended, windowed topo; 
% M is the topo/surface after the detrend and windowed process;
Origin_Output=P; % save the original output, complex value;
%% %%%%%%%%%%%%%%%%%%%%%%%%
% 3. ANALYZE THE SPECTRUM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D Power Spectrum
% Plot the 1D version of the spectrum (Fig. 5). We'll plot the 2D spectrum 
% a few steps later, but for now it's easier to visualize in 1D.
s1d = figure('Name','Fig. 3: 1D Spectrum','NumberTitle','off');
% subplot(2,1,1)
[axf,axw] = SpecPlot1D(fv,Pv); % plot the 1D power spectrum;
nbin = 40; % Number of bins;
B = bin(log10(fv),log10(Pv),nbin,0); % Bin the log-transformed data
% Plot the binned values
hold(axf,'on')
plot(axf,10.^B(:,1),10.^B(:,2),'ok','markerfacecolor','w') % Plot the background power spectrum;
hold on;
plot(axf,10.^B1(:,1),10^fit_Background(1)*(10.^B1(:,1)).^fit_Background(2),'b'); % 10^(A1+A2logx)=y, the background power spectrum line;
%% Quantify the good of fitness (find out the roll-off frequency)

% Use r^2 to find the frequency roll-off;
s3d = figure('Name','Fig. 4: Goodness of fit','NumberTitle','off');
for i=6:nbin
fit_analysis = robustfit(B(4:i,1),B(4:i,2)); % start from the fourth point;
r2(i)=rsquared(B(4:i,2),fit_analysis(1)+B(4:i,1).*fit_analysis(2)); % calculate each r-square;
end
x=7:nbin; % plot start from the 7 point;
scatter(x,r2(7:nbin)); % scatter plot the r2;
line(x,r2(7:nbin)); % add a line to the scatter plot;
title('r-square to find the roll-off point','FontSize',25); % Add the title;
xlabel('The serial number of the break point in the fit', 'FontSize',25); % Add the x-label;
xt=get(gca,'XTick');
set(gca,'FontSize',20); % enlarge the XTick fontsize;
ylabel('r^2','FontSize',25); % Add the y-label;
xt=get(gca,'YTick');
set(gca,'FontSize',20); % enlarge the YTick fontsize;

% Use two RMSE to find the frequency roll-off;
for i=6:38
fit_analysis1 = robustfit(B(4:i,1),B(4:i,2));
fit_analysis2 = robustfit(B(i:40,1),B(i:40,2));
rmse1(i)=rmse(B(4:i,2),fit_analysis1(1)+B(4:i,1).*fit_analysis1(2));
rmse2(i)=rmse(B(i:40,2),fit_analysis2(1)+B(i:40,1).*fit_analysis2(2));
rmse_t(i)=rmse1(i)+rmse2(i);
end
x=6:38;
scatter(x,rmse_t(6:38));  % scatter plot the rmse;
line(x,rmse_t(6:38)); % add a line to the scatter plot;
title('rmse to find the roll-off point','FontSize',25);
xlabel('The serial number of the break point in the fit', 'FontSize',25);
xt=get(gca,'XTick');
set(gca,'FontSize',20);
ylabel('rmse','FontSize',25);
xt=get(gca,'YTick');
set(gca,'FontSize',20);


% Determine the H parameter
fit_data= robustfit(B(4:34,1),B(4:34,2)); % Fit the background line;
beta=-fit_data(2,1); % Calculate the slope of the fit line-Beta;
D=(8-beta)/2; % Calculate the Fractal Dimension;
H=3-D; % Calculate the Hurst Parameter/Surface Roughness;
%% Background spectral (Artifial randomly rough surfaces generation)
Variance=sum(Pm(:)); % the sum of the total Pdft array;
std=sqrt(Variance); % Calculate the standard deviation;
size_Z=size(fm); % get size of the frequency matrix;
n=size_Z(1,1); % get the number of pixels in y;
m=size_Z(1,2); % get the number of pixels in x;
Lx=size_Z(1,2); % get the size of the topo in x-direction;
size_f=size(fv); % get the size of the frequency vector;
Pv_t=zeros(size_f(1,1),1); % Create vacant matrix to hold the total/sum Pv of generated randomly surfaces;
Pm_t=zeros(size_Z(1,1),size_Z(1,2)); % Create vacant matrix to hold the total/sum Pm of generated randomly surfaces;
loop=5; % generate 5 randomly rough surface;
Pm_art=NaN(size_Z(1,1),size_Z(1,2)); % Create vacant matrix to hold the Pm of each generated randomly surface;
Pv_art=NaN(size_f(1,1),1); % Create vacant matrix to hold the Pv of each geenrated randomly surface;
z_art_mat = zeros(size_Z(1,1), size_Z(1,2), loop); % Create vacant matrix to hold the DEM of each generated randomly surface;
% Generate severral random surfaces;
for i=1:loop
[z_art, PixelWidth, PSD] = artificial_surf(std, H, Lx, m, n); 
z_art_mat(:,:,i)=z_art;
 pad = 0;   
 window = 0; 
[P_BG, Pm_BG, fm_BG, Pv_BG, fv_BG] = fft2D(z_art_mat(:,:,i),dx,dy,pad,window); 
 Pm_art(:,:,i)=Pm_BG;
 Pm_t=Pm_t+Pm_art(:,:,i);
 Pv_art(:,:,i)=Pv_BG;
 Pv_t=Pv_t+Pv_art(:,:,i);
end
Pm_average=Pm_t./loop; % calculate the average 2D power spectrum of the generated randomly surfaces;
Pv_average=Pv_t./loop; % calculate the average 1D power vector of the generated randomlu surfaces;

% Calculate the average artificial topo and shaded plot it
z_art_mat_total=zeros(size_Z(1,1),size_Z(1,2));
for i=1:loop
    z_art_mat_total=z_art_mat_total+z_art_mat(:,:,i);
end
z_art_mat_average=z_art_mat_total/loop;

% Cut the average generated surface to the size the same as the original
% size and plot it;
Cut_averaged_generated_surface=z_art_mat_average((size_Z(1,1)/2)-((Ny/2)-1):((size_Z(1,2)/2)+Ny/2),(size_Z(1,1)/2)-((Nx/2)-1):((size_Z(1,2)/2)+Nx/2));
figure('Name','Fig. 5: Average generated surface','NumberTitle','off')
ShadePlot(1:Nx,1:Ny,Cut_averaged_generated_surface); % Shade plot the cutted averaged generated random surfaces;

% Plot the 1D power spectrum of background;
s2d = figure('Name','Fig. 6: 1D Spectrum of background','NumberTitle','off');
Pv_average(Pv_average<10^-15)=NaN; % Get rid of the very small value of Pv of generated randomly surfaces;
[axf_BG, axw_BG] = SpecPlot1D(fv_BG,Pv_average); % Plot 1D power spectrum figure of the generated randomly surfaces;
hold(axw_BG,'on')
nbin_BG = 40; % Number of bins
B1 = bin(log10(fv),log10(Pv_average),nbin_BG,0); % Calculate the statistics associated with each bin;

% Plot the binned values
hold(axf_BG,'on')
plot(axf_BG,10.^B1(:,1),10.^B1(:,2),'ok','markerfacecolor','w') % Plot the binned values;

% Fit a trend with the form P ~ 1/f^n, and plot it
fit_Background = robustfit(B1(:,1),B1(:,2)); % returns the 2 coefficient. A1+A2logx=logy;
plot(axf_BG,10.^B1(:,1),10^fit_Background(1)*(10.^B1(:,1)).^fit_Background(2),'k'); % 10^(A1+A2logx)=y;
%% Normalized spectra (normalized by generated surfaces)
Pv_average(Pv_average<10^-15)=NaN; % get the rid of the samll Pv_average;
Pv_Nor=Pv./Pv_average; % Calculate the Normalized Pv;

alpha=0.05; % set the 95% threshold;
th_5 = 0.5*icdf('Chisquare',1-alpha,2); % calcualte the threshold;
s4d = figure('Name','Fig. 7: Normalized Power','NumberTitle','off');
[axf_Nor,axw_Nor] = NormaSpecPlot1D(fv,Pv_Nor); % plot the Normalized spectra power;
hold(axf_Nor,'on');
plot(axf_Nor,[10^-04 1],[th_5 th_5],'LineWidth',2); % Plot the threshold;
%% 2D Power Spectrum
% Use this fit to normalize the 2D spectrum
Pm_average(Pm_average<10^-15)=NaN; % get rid of the small value of the Pm_average;
Pmn = Pm./Pm_average; % Normalized 2D power spectrum;
p=chi2cdf(Pmn,2); % use the chi-square distribution with degree of 2 to find the corresponding probability;
figure('Name','Fig. 8: Original 2D Spectrum','NumberTitle','off')
SpecPlot2D(fm,p); % Plot the 2D normalized power spectrum probability;

% Smooth the 2D power spectrum
N = 60; % The size of the kernel;
sigma = 3; % The sigma;
imagesc(gaussian2d_kernel(N,sigma)); % plot the kernel;
filtered_signal=conv2(p,gaussian2d_kernel(N,sigma),'same'); % filter the signal;
filtered_signal = filtered_signal/max(filtered_signal(:));  % Normalize the filtered signal;
filtered_signal(filtered_signal<.50) = nan; % Assume the threshold is 0.90;
figure('Name','Fig. 9: Filtered 2D Spectrum','NumberTitle','off')
SpecPlot2D(fm,filtered_signal()); % Plot the 2D filtered power spectrum;
%% Calculate the corresponding significant wavelength
Index=~isnan(filtered_signal); % The logic index of the 2D filtered power spectrum;
filtered_fm=Index.*fm; % Apply the logic index to the original frequency matrix;
filtered_fm(filtered_fm==0)=nan; % transfer the 0 to nan; 
filtered_wavelength=1./filtered_fm; % Calculate the significant wavelength;
figure('Name','Fig. 10: Significant wavelength','NumberTitle','off')
SpecPlot2D(fm,filtered_wavelength()); % Plot the 2D filtered wavelength;
%% Use the filtered 2D power spectra to reconstruct the surface and compare with the original topo
% reconstruct the original surface using the before filtered 2D power
% spectrum;
Origin_Surface_P=ifftshift(P); % Shift back of the original P to the original FFT output;
reconstruc_original_topo=ifft2(Origin_Surface_P); % inverse fourier transform to reconstruct the surface;
figure('Name','Fig. 11: reconstructed original surface','NumberTitle','off')
ShadePlot(topopiece_vector.x,topopiece_vector.y,reconstruc_original_topo); % plot the reconstrcted surface;

% Plot the detrended and windowed topography;
figure('Name','Fig. 12: The detrended and windowed original surface','NumberTitle','off')
ShadePlot(topopiece_vector.x,topopiece_vector.y,M); % plot the reconstrcuted surface;

% Calculate the relative difference between the detrended and windowed
% topo and the ifft2 topo of detrended and windowed topo;
figure('Name','Fig. 13: The difference between the detrended and windowed original surface and the ifft2 of the detrended and windowed original surface fft2 coefficient','NumberTitle','off')
topo_dif=reconstruc_original_topo-M; % Substract the two DEM to get the difference;
Relative_Dif=topo_dif./M; % Calculate the relative difference between the two topo/surface;
% Contour plot the difference;
contour(Relative_Dif);
colorbar;

% reconstruct the surface based on filtered 2D power spectrum
filtered_P=Index.*P; % Filter all the other unsignificant frequencies;
filtered_P_shiftback=ifftshift(filtered_P); % Shift back to the original FFT output;
reconstruc_topo=ifft2(filtered_P_shiftback); % inverse fourier transform to reconstruct the surface;
figure('Name','Fig. 13: reconstructed surface','NumberTitle','off')
ShadePlot(topopiece_vector.x,topopiece_vector.y,reconstruc_topo); % plot the reconstrcuted filtered surface;

% Calculate the relative difference between the reconstructed 
% topo based on identified significant frequencies and the ifft2 topo of detrended and windowed topo;
figure('Name','Fig. 13: The difference between the reconstructed surface based on identified significant frequencies and the ifft2 of the detrended and windowed original surface fft2 coefficient','NumberTitle','off')
topo_dif_1=reconstruc_topo-reconstruc_original_topo; % Substract the two DEM to get the difference;
Relative_Dif_1=topo_dif_1./reconstruc_original_topo; % Calculate the relative difference between the two topo/surface;
% Contour plot the difference between two topo/surfaces;
contour(Relative_Dif_1);
colorbar;



