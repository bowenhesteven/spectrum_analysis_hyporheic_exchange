%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%       SPECTOP 2.0       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%* Matlab function, Spectop, that determines amplitude coefficients from map 
% data of the Spectral(Fourier) Series given by Worman et al., Gephysical 
% Research Letters, Submitted.
%
%* FUNCTION INPUT - DATA ON HYDRAULIC POTENTIAL SURFACE 
%
% The topological data defines discrete points in orthogonal coordinates with 
% an inner area with higher density of points and outer areas with sparser
% points as shown in Figure 1 (Wšrman, A., A. I. Packman, L. Marklund, J. W. 
% Harvey, and S. H. Stone (2006), Exact three-dimensional spectral solution 
% to surface-groundwater interactions with arbitrary surface topography, 
% Geophys. Res. Lett., 33, L07402, doi:10.1029/2006GL025747).
% 
% The following 3 variables should be available in the functional call: 
% 1. X = Cartesian coordinate with dimensions (1:M_x)
% 2. Y = Cartesian coordinate with dimensions (1:M_y)
% 3. Z = Square matrix with function values data corresponding to Z(X,Y) and 
% the dimensions (1:M_y,1:M_x). Z can be thought of as topographical (map)
% data. Note that the Matlab convention is that x correpsonds to the columns
% of Z and y corresponds to the rows of Z. 
%
% The density of map points does not have to be uniform as can be seen in
% Figure 1 in the paper by Worman et al.
% 
%* FUNCTION INPUT - NUMERICAL PARAMETERS AND SPECTRUM DEFINITION
%
% The code uses a direct method to determine Fourier coefficients from the 
% known topography Z(Y,X). 
% The amplitudes, paramhat, are determined by the initially set of wavelengths, 
% lambda, by a non-linear distribution:
% lambda = f2 + maxlambda*(i/N)^f1, where 
% f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
% f2 = a lag avoiding the smallest wavelengths, i.e. smaller than the map 
% resolution. Normally, f2 should cover at least two map points (Nyquist frequency). 
% maxlambda=f3*((max(X)-min(X))+(max(Y)-min(Y)))/2
% f3 = A factor that gives the maximum wavelength as a factor of the mean
% size of the domain. As compared with the "fundamental Fourier frequency", 
% f3 should be about 1 without "mirroring" and 2 with "mirroring", 
% which results in a maximum wavelength equal to the size of the domain. 
% N = The total number of wavelengths in each of the x- and y-directions.
% periodic = Variable that defines if the Fourier series should be periodic
% with a wavelength equal to the size of the domain. periodic = 1 causes
% such setting of the wavelengths of each term. If the series is not
% periodic this could cause undesireable functional shapes outside the
% domain. These can be counteracted using the mirroring technique. See
% below. If periodic is set to unity, it is not necessary to use the
% mirroring technique.                                                                                                                                       
%
% One option is to assign the same wavelengths in both x- and y-directions,
% whereby each wavelength in the x-direction is combined with all wavelengths
% in the y-directions. This implies a wide range of wavelength ratios in the
% Fourier terms. The wavelength ratio defines the aspect ratio of the hill-type
% harmonic function. 
% Another option is to limit the aspect ratios of the hill functions. Spectral
% analyses of landscapes indicate that the predominant hill-type harmonic functions 
% are limited to a certain interval that carries most of the primary information.
% The following parameters define the two options:
% cross = If cross equals 0, the first option is activated involving hill-type
% functions of all apsect ratios. If cross > 0 (integer), cross defines the minimum 
% number of y-wavelengths associated with each wavelengths in the x-direction.
% ratio_low = For the second option, defines the lower limit of lambda_y/lambda_x that
% are used to define the hill-type functions, where lambda is wavelength
% ratio_high = For the second option, defines the upper limit of lambda_y/lambda_x
%
% This version includes a mirroring procedure. That is described in the
% paper by Worman et al. This technique impose restrictions on the
% amplitude coefficients in a manner that avoids undesirable functional
% shapes outside the area of interest. 
% To "turn on" this option, the variable mirror should be unity in the
% functional call.
%
%* FUNCTION OUTPUT - FOURIER SURFACE AND NUMERICAL PARAMETERS
%
% The 6 return parameters of the routine comprise
% 1. paramhat = Amplitude coefficients with dimension (1:NxN)
% 2. lambda = Wavelength (same in both x- and y-directions) with dimension
% (1xN)
% 3. kx = Wave number, 2*pi/lambda, in x-direction with dimension (1:NxN)
% 4. ky = Wave number in  y-direction with dimension (1:NxN)
% 5. h = Fourier fitted surface that corresponds to the coordinates X and
% Y
% 6. mean_error = Mean error relative to the difference (max(Z) - min(Z)).
%
%  Copyright: Anders Worman
%             Environmental Physics Group, SLU
%             email: anders.worman@bt.slu.se
%  10 Jan 2006, MATLAB 7.1 version 
%  IF YOU PUBLISH WORK BENEFITING FROM THIS M-FILE, PLEASE CITE IT AS:
%    Worman, A. (2006) Spectop: A matlab function for numerical 
%    spectral analysis of landscape topography, 
%    ftp://ftp.agu.org/apend/" (Username = "anonymous", Password = "guest")


%%% 0: Start by checking dimensions of matrices

% Row or column vectors?
function [paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X,Y,Z,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic)
if size(X,2)<size(X,1), X=X';end % X=X' means flip the X vector;
if size(Y,2)<size(Y,1), Y=Y';end % Y=Y' means flip the Y vector;

% Check dimensions of X and Y versus Z
if length(X)~=size(Z,2) | length(Y)~=size(Z,1) % if the length of X is not the length of 1:M_x or the length of Y is not the length of 1:M_y;
error('Dimensions of x(1:M_x) and y(1:M_y) used in the call to Spectop.m do not fit the columns and rows of Z(1:M_y,1:M_x)') % show the error of the size of the matrix Z(Y,X);
end        

%%% A) Optional mirroring of data
% First calculate resolution for later use; find the smallest resolution of
% the data;
minvalue=1000000; % set up a min value of the resolution;
for i=1:length(X)-1 % for each position of the x direction;
    A=abs(X(i+1)-X(i)); % for each difference between two adjacent measurement grid in x-direction;
    if A<minvalue     % find that the smallest value/smallest resolution;
        minvalue=A;
    end
end

resol=2*minvalue; % let the resolution be the twice of the smallest resolution found above; based on Nyquist frequency theorm;
f2=resol; % let f2 be this resolution just calculated above;
if mirror==1
    X_orig=X; % save the original X to X_orig;
    Y_orig=Y; % save the original Y to Y_orig;
    Z_orig=Z; % save the original Z to Z_orig;
    NX=size(X,2); % the number of the grid in x-direction;
    NY=size(Y,2); % the number of the grid in y-direction;
    step=10; % Resolution of mirroring (Number of in mirroring = floor(NX/step))

    % First, twice in positive x-direction;
    %%% Mirroring towards positive X, X=151,152,153,154....165;
    max_X=max(X); % get the maxmum value of X;
    min_X=min(X); % get the minimum value of X;
    count=0;
    for i=1:step:NX-1 % 1:10:149?
    	count=count+1;
        X(NX+count)=2*max_X-X(NX-i);% X(150+1)=2*max_X-X(149), X(150+2)=2*max_X-X(139).....X(150+15))=2*max_X-X(9);
        Z(1:NY,NX+count)=Z(1:NY,NX-i);% Z(1:150,151)=Z(1:150,149), Z(1:150,152)=Z((1:150,139).....Z(1:150,165)=Z((1:150,9);
    end

    %%% Mirroring towards negative X
    % Number of first mirrored nodes (i.e. count after loop or if step=1,      
    % NX_2=NX-1)
    NX_2=size(X,2)-NX; % the difference between the length of new X vector and the original X vector;
    max_X_2=max(X); % get the new maximum value of X;
    % Leave room in the X-vector for the lower X-values
    X(2*NX_2+1:3*NX_2+NX)=X(1:NX+NX_2); % X(31:195)=X(1:150+15);
    Z(1:NY,2*NX_2+1:3*NX_2+NX)=Z(1:NY,1:NX+NX_2); % Z(1:150,31:195)=Z(1:150,1:150+15);
    % Since the first mirroring may have applied a resolution different from 1,
    % the original and the first mirrored field are mirrored separately
    count=0;
    for i=2*NX_2+1:step:2*NX_2+NX-1   % 31:10:179; 
    	count=count+1;
    	X(2*NX_2+1-count)=2*min_X-X(i); % X(30)=2*min_X-X(31), X(29)=2*min_X-X(41), X(28)=2*min_X-X(51)....X(16)=2*min_X-X(171);
    	Z(1:NY,2*NX_2+1-count)=Z(1:NY,i); % Z(1:150,30)=Z(1:150,31), Z(1:150,29)=Z(1:150,41), Z(1:150,28)=Z(1:150,51)....Z(1:150,16)=Z(1:150,171);
    end
    count=0;
    for i=2*NX_2+NX+1:3*NX_2+NX % 181:195;
    	count=count+1;
    	X(NX_2+1-count)=2*min_X-X(i); % X(15)=2*min_X-X(181), X(14)=2*min_X-X(182), X(13)=2*min_X-X(183), X(12)=2*min_X-X(184).....X(1)=2*min_X-X(195);
    	Z(1:NY,NX_2+1-count)=Z(1:NY,i); % Z(1:150,15)=Z(1:150,181),  Z(1:150,14)=Z(1:150,182), Z(1:150,13)=Z(1:150,183), Z(1:150,12)=Z(1:150,184)...Z(1:150,1)=Z(1:150,195);
    end
    
    % Twice in y-direction
    %%% Mirroring towards positive Y
    max_Y=max(Y); % get the maxmum value of Y;
    min_Y=min(Y); % get the minmum value of Y;
    count=0;
    for i=1:step:NY-1 % 1:10:149;
    	count=count+1;
    	Y(NY+count)=2*max_Y-Y(NY-i); % Y(151)=2*max_Y-Y(149), Y(152)=2*max_Y-Y(139), Y(153)=2*max_Y-Y(129).....Y(165)=2*max_Y-Y(9);
    	Z(NY+count,1:NX)=Z(NY-i,1:NX); % Z(151,1:150)=Z(149,1:150), Z(152,1:150)=Z(139,1:150), Z(153,1:150)=Z(129,1:150).... Z(165,1:150)=Z(9,1:150);
    end

    %%% Mirroring towards negative Y
    %Number of first mirrored nodes (i.e. count after loop or if step=1, 
    % NX_2=NX-1)
    NY_2=size(Y,2)-NY; % the difference between the length of new Y vector and the original Y vector;
    max_Y_2=max(Y); % get the new maximum value of Y;
    % Leave room in the X-vector for the lower X-values
    Y(2*NY_2+1:3*NY_2+NY)=Y(1:NY+NY_2); % Y(31:195)=Y(1:150+15);
    Z(2*NY_2+1:3*NY_2+NY,1:NX)=Z(1:NY+NY_2,1:NX); % Z(31:195,1:150)=Z(1:165,1:150);
    % Since the first mirroring may have applied a resolution different from 1,
    % the original and the first mirrored field are mirrored separately
    count=0;
    for i=2*NY_2+1:step:2*NY_2+NY-1  % 31:10:179; 
    	count=count+1;
    	Y(2*NY_2+1-count)=2*min_Y-Y(i); % Y(30)=2*min_Y-Y(31), Y(29)=2*min_Y-Y(41), Y(28)=2*min_Y-Y(51)....Y(16)=2*min_Y-Y(171);
    	Z(2*NY_2+1-count,1:NX)=Z(i,1:NX); % Z(30,1:150)=Z(31,1:150), Z(29,1:150)=Z(41,1:150), Z(28,1:150)=Z(51,1:150)....Z(16,1:150)=Z(171,1:150);
    end
    count=0;
    for i=2*NY_2+NY+1:3*NY_2+NY % 181:195;
    	count=count+1;
    	Y(NY_2+1-count)=2*min_Y-Y(i);  % Y(15)=2*min_Y-Y(181), Y(14)=2*min_Y-Y(182), Y(13)=2*min_Y-Y(183).....Y(1)=2*min_Y-Y(195);
    	Z(NY_2+1-count,1:NX)=Z(i,1:NX); % Z(15,1:150)=Z(181,1:150),  Z(14,1:150)=Z(182,1:150), Z(13,1:150)=Z(183,1:150), Z(12,1:150)=Z(184,1:150)...Z(1,1:150)=Z(195,1:150);
    end

end

%% B) Definition of lambda-spectrum included in Fourier series
av=mean(mean(Z)); % Average of Z, the mean value of the new domain topo;
maxlambda=f3*((max(X)-min(X))+(max(Y)-min(Y)))/2; % 2*((2002500-312740)+(6979500-5260300))/2=3408960;
% Start by selecting wavelengths that are integer parts of maxlambda
maxnum=floor(maxlambda/resol); % maximum times one wavelength can fit into the domain exactly, floor(5681.6)
for i=1:maxnum % 1:5681;
    lambda_pot(maxnum+1-i)=(maxlambda+resol)/i; % lambda_pot(5681)=(2.5641e+06+600)/1,lambda_pot(5680)=(2.5641e+06+600)/2.....lambda_pot(1)=(2.5641e+06+600)/4273
end
%% Distribute wavelengths
if maxnum<N % N = The total number of wavelengths in each of the x- and y-directions.
    output=['maximum number of wavelength that produce a periodic function on the domain is'...
        ,int2str(maxnum),'. The program uses this number and continue'];
    warning(output)
end
%% Assign lambda values and select those that coincide with closest value in lambda_pot. 
count=0;
% Calculate the lambda_test, total of N=28;
for i=1:N % 1:28;
	lambda_test(i)=resol+maxlambda*(i/N)^f1; % lambda_test(1)=600+3.4090e+06*(1/28)^3.5.....lambda_test(28)=600+3.4090e+06*(28/28)^3.5
    % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
    if periodic==1
    k=1;
    while lambda_pot(k)<lambda_test(i)
        k=k+1;
    end
    if k==1
        count=count+1; % count=0+1;
        lambda(count)=lambda_pot(k); % lambda(1)=lambda_pot(1);
    elseif abs(lambda_pot(k-1)-lambda_test(i))<abs(lambda_pot(k)-lambda_test(i))
        count=count+1;
        lambda(count)=lambda_pot(k-1);
        if count>1 
            if lambda(count)==lambda(count-1),count=count-1;
            end
        end
    else
        count=count+1;
        lambda(count)=lambda_pot(k);
        if count>1
            if lambda(count)==lambda(count-1),count=count-1;
            end 
        end
    end
    else
        lambda(i)=lambda_test(i); % when periodic==0, only this sentence works; that means lambda=lambda_test;
    end
end

% The routine allows the same lambda-value in several positions. These are
% now deleated
count=1;
clear lambda_test
lambda_test=lambda;
clear lambda
lambda(count)=lambda_test(1);
for i=2:length(lambda_test) % 2:28;
    if lambda_test(i)==lambda(count)
    else
        count=count+1;
        lambda(count)=lambda_test(i);
    end
end
%if lambda(length(lambda))==lambda(length(lambda)-1),
%    help(1:length(lambda)-1)=lambda(1:length(lambda)-1);clear lambda,lambda=help;clear help
%end
N=length(lambda); % N=28;
for i=1:N % 1:28;
    wavenum(i)=2*pi/lambda(i); % wavenum(1)=2*pi/lambda(1),wavenum(2)=2*pi/lambda(2)....wavenum(28)=2*pi/lambda(28);
end
% sets wave numbers
if cross>0 %cross = If cross equals 0, the first option is activated involving hill-type
% functions of all apsect ratios. If cross > 0 (integer), cross defines the minimum 
% number of y-wavelengths associated with each wavelengths in the x-direction.
	if f2>0
		N_2=round(N*ratio_high/(maxlambda/f2));
		if N_2<cross, N_2=cross;end
		if N_2>N, N_2=N;end
	else
		N_2=cross;
	end
	for i=1:N
		lam_min=lambda(i)*ratio_low;
		lam_max=lambda(i)*ratio_high;
		for j=1:N_2
		kx((i-1)*N_2+j)=wavenum(i);
		lambda_2=lam_min+lam_max*(j/N_2)^f1;
		wavenum_2=2*pi/lambda_2;
		ky((i-1)*N_2+j)=wavenum_2;
        end
	end
else % cross=0; ONLY this following part works;
	for i=1:N % 1:28;
		for j=1:N % 1:28;
		kx((i-1)*N+j)=wavenum(i); % kx(1)=wavenum(1), kx(2)=wavenum(1)...kx(28)=wavenum(1);
        % kx(28+1)=wavenum(2),kx(28+2)=wavenum(2)...kx(28+28)=wavenum(2);
        % kx(28*28);
		ky((i-1)*N+j)=wavenum(j); % ky(1)=wavenum(1), ky(2)=wavenum(2)...ky(28)=wavenum(28);
        % ky(28+1)=wavenum(1),ky(28+2)=wavenum(2)...ky(28+28)=wavenum(28);
        % ky(28*28);
		end
	end
end

%%% C) Calculation of coefficient matrix-H
% Optional transformation that eliminates negative amplitudes if pp=2; 
% If pp = 1, no transformation is applied. pp should be > 1. 
% Transformation:
% h-av=(A^pp)^(1/pp) Sin(kx x) Cos(ky y) => ((h-av)^pp)^(1/pp) = 
% A (Sin(Kx x)^pp)^(1/pp) (Cos(ky y)^pp)^(1/pp) ;  
% In this way, (h-av), cosines and sines are squared before square rooted. 
% A is squared first to avoid square roots of negative values.
pp=1; 


NX=size(X,2);
NY=size(Y,2);
for j=1:NY % 1:195;
	vect((j-1)*NX+1:(j-1)*NX+NX)=((Z(j,1:NX)-av).^pp).^(1/pp); % vect=Z;
    % vect(1:195)=((Z(1,1:195)-av).^1).^(1/1),
    % vect(196:195+195)=((Z(2,1:195)-av).^1).^(1/1),
    % vect(2*195+1:2*195+195)=((Z(3,1:195)-av).^1).^(1/1).....
    % vect(194*195+1:194*195+195)=((Z(195,1:195)-av).^1).^(1/1)
	coef((j-1)*NX+1:(j-1)*NX+NX,1:size(kx,2))=((sin(transpose(X(1:NX))*...
        kx(1:size(kx,2))).^pp).^(1/pp)).*((cos(Y(j)*ones(NX,1)*ky(1:size(kx,2))).^pp).^(1/pp)); % coef=C;
    % coef(1:195,1:784)=((sin(transpose(X(1:195))*kx(1:784))).^1).^(1/1)).*((cos(Y(1)*ones(195,1)*ky(1:784)).^1).^(1/1)),[195*784]
    % coef(1*195+1:1*195+195,1:784)=((sin(transpose(X(1:195))*kx(1:784))).^1).^(1/1)).*((cos(Y(2)*ones(195,1)*ky(1:784)).^1).^(1/1)),[(2*195)*784]
    % coef(2*195+1:2*195+195, 1:784)=((sin(transpose(X(1:195))*kx(1:784))).^1).^(1/1)).*((cos(Y(3)*ones(195,1)*ky(1:784)).^1).^(1/1))
    %....
    % coef(194*195+1:194*195+195, 1:784)=((sin(transpose(X(1:195))*kx(1:784))).^1).^(1/1)).*((cos(Y(195)*ones(195,1)*ky(1:784)).^1).^(1/1))
end


%% D) Solution to the linear problem with direct method
paramhat=((coef')*coef)\((coef')*vect');
paramhat=(paramhat.^pp).^(1/pp); % This is the final calculated H, amplitudes hm;

%%% If mirroring has been applied:
if mirror==1 
    clear X Y Z NX NY
    X=X_orig;Y=Y_orig;NX=size(X,2);NY=size(Y,2);
    Z=Z_orig;
end
%% E) Calculates surface and error estimates
for j=1:NY % 1:150;
	X2D(j,1:NX)=X(1:NX); % X2D(1,1:150)=X(1:150), X2D(2,1:150)=X(1:150).....X2D(150,1:150)=X(1:150);
end
for i=1:NX %1:150;
	Y2D(1:NY,i)=transpose(Y(1:NY)); % Y2D(1:150,1)=transpose(Y(1:150)),Y2D(1:150,2)=transpose(Y(1:150))...Y2D(1:150,150)=transpose(Y(1:150))
end

h=av*ones(NY,NX); % h=av*ones(150,150);
%h=av*ones(NX,NY);

for p=1:size(paramhat,1) % p=1:784;
	h=h+paramhat(p)*sin(kx(p)*X2D).*cos(ky(p)*Y2D); % h=average_h+paramhat(1)*sin(kx(1)*X2D).*cos(ky(1)*Y2D), this is the Fourier Series Fitted Topography;
end
clear X2D Y2D

aaa=h-Z;
mean_error=mean(mean(abs(aaa)))/abs(min(min(Z))-max(max(Z)));