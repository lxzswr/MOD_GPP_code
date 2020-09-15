%MODIS scalar for VPD
function y=fvpd(x,vmin,vmax)

    y=1-(x-vmin)./(vmax-vmin);
    y(y<0)=0;
    y(y>1)=1;

end