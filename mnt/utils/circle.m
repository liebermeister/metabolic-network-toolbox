function circle(x,y,r,col)

% circle(x,y,r,col)

xlist = x+r*cos(0:0.01:2*pi);
ylist = y+r*sin(0:0.01:2*pi);

fill(xlist,ylist,col,'EdgeColor',col);