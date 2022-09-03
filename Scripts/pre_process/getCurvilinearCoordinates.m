function[] = getCurvilinearCoordinates(pathData,nameData)

pathData = '/Users/gomezvjd/Downloads/dataBedforms/';
nameData = 'topoVector_subset.mat';

data = load(strcat(pathData,nameData)); 
data = data.topo_vector_subset;

x = data.x;
y = data.y;
z = data.z;
ind = 1:length(x);

n_intervals = 1000;
[xc,yc] = getCentroids(x,y,n_intervals);
xc_all = mean(x);
yc_all = mean(y)

xv_opt = fitCircle(x,y,xc_all,yc_all);
R_opt2 = mean((x-xv_opt(1)).^2 + (y-xv_opt(2)).^2);
[xunit,yunit] = circle(xv_opt(1),xv_opt(2),sqrt(R_opt2));

[x_l,ix] = min(x);
y_l = (y(ix));

[y_r,ix] = min(y);
x_r = (x(ix));

theta_start = atan((y_l - xv_opt(2))/(x_l - xv_opt(1)))
theta_end = atan((y_r - xv_opt(2))/(x_r - xv_opt(1)))

theta = atan((y - xv_opt(2))./(x - xv_opt(1)));

fig1 = figure(1);
plot(x,y,'.k')
hold on
%plot(xc,yc,'or')
plot(xc_all,yc_all,'om')
plot(xv_opt(1),xv_opt(2),'og')
plot(xunit,yunit,'-g')
plot(x_l,y_l,'og')
plot(x_r,y_r,'og')

hold off

fig2 = figure(2);
clf(2)
plot(x,theta,'.k')
hold on



XX = [x y];
[coeff,score,latent,tsquared,explained,mu] = pca(XX);
xtmp = linspace(min(x),max(x),100);
ytmp1 = (coeff(2,1)/coeff(1,1))*xtmp;
ytmp2 = (coeff(2,2)/coeff(1,2))*xtmp;
 
figure(10)
clf(10)
plot(x,y,'ok'), hold on
plot(xtmp+mean(x),ytmp1+mean(y),'-b')
plot(xtmp+mean(x),ytmp2+mean(y),'-g')
plot(mean(x),mean(y),'or')
axis('equal')
hold off

x_t = score(:,1);
y_t = score(:,2);
n_intervals = 100;
[x_tc,y_tc] = getCentroids(x_t,y_t,n_intervals);

figure(11)
clf(11)
plot(x_t,y_t,'.k')
hold on
plot(x_tc,y_tc,'or')
hold off



%     K = cov(XX); % Covariance matrix
%     [V,D] = eig(K);
    
%     
%     
%     L_pca_min(i) = max(score(:,2))-min(score(:,2));
%     L_pca_max(i) = max(score(:,1))-min(score(:,1));
 
%     xtmp = linspace(pond_Xmin(i),pond_Xmax(i),100)-mean(xb);
%     ytmp1 = (coeff(2,1)/coeff(1,1))*xtmp;
%     ytmp2 = (coeff(2,2)/coeff(1,2))*xtmp;



end

function[xc,yc] = getCentroids(x,y,n_intervals)

[xs,ix] = sort(x);
ys = y(ix);


xbreaks = min(x):(round(max(x) - min(x))/n_intervals):max(x);

xc = nan(1,n_intervals);
yc = nan(1,n_intervals);

for(ii = 1:n_intervals)
    
    in = find(x>=xbreaks(ii) & x <= xbreaks(ii+1));
    
    xc(ii) = median(x(in));
    yc(ii) = median(y(in));
    
end



end

function[xunit,yunit] = circle(xcenter,ycenter,r)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + xcenter;
yunit = r * sin(th) + ycenter;
end
