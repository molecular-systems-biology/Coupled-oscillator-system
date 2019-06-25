function [dth_dt] = kuramoto(t,th,params)
%KURAMOTO Summary of this function goes here
%   t: time 
%   th: theta 
%   params: Parameters 

    dth_dt = zeros(4,1);
    dth_dt(1) = params(1) - params(5)*sin(th(1) - th(2)) - params(6)*sin(th(1) - th(3)) - params(7)*sin(th(1) - th(4));
    dth_dt(2) = params(2) - params(8)*sin(th(2) - th(1)) - params(9)*sin(th(2) - th(3)) - params(10)*sin(th(2) - th(4));
    dth_dt(3) = params(3) - params(11)*sin(th(3) - th(1)) - params(12)*sin(th(3) - th(2)) - params(13)*sin(th(3) - th(4));
    dth_dt(4) = params(4) - params(14)*sin(th(4) - th(1)) - params(15)*sin(th(4) - th(2)) - params(16)*sin(th(4) - th(3));
    
end

