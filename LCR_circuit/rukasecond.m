function [ x_it,y_it ] = rukasecond( x, y, t, h, func)
%RK_4_2 Step of 4th order Runge-Kutta algorithm for second order ODE
%   Detailed explanation goes here
    k_y1=func(x, y, t);
    k_x1=y;
    k_y2=func(x+0.5*h*k_x1, y+0.5*h*k_y1, t+0.5*h);
    k_x2=y+0.5*h*k_y1;
    k_y3=func(x+0.5*h*k_x2, y+0.5*h*k_y2, t+0.5*h);
    k_x3=y+0.5*h*k_y2;
    k_y4=func(x+h*k_x3, y+h*k_y3, t+h);
    k_x4=y+h*k_y3;
    
    x_it=x+h*(k_x1+2*(k_x2+k_x3)+k_x4)/6;
    y_it=y+h*(k_y1+2*(k_y2+k_y3)+k_y4)/6;
end

