%MODIS scalar for minimum temperature
function y=ftair(x,tmin,tmax)

    y=(x-tmin)./(tmax-tmin);
    y(y<0)=0;
    y(y>1)=1;

end