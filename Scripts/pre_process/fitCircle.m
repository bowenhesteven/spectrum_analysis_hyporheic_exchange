function[xv_opt] = fitCircle(x,y,xc_all,yc_all)

xini = [xc_all,yc_all]; % Starting guess
[xv_opt,fval,exitflag,output] = fminunc(@(xv)circleFun(xv,x,y),xini);

end 

function[CV] = circleFun(xv,x,y)

R_est2 = (x-xv(1)).^2 + (y-xv(2)).^2;
CV = nanstd(R_est2)/nanmean(R_est2);

end

