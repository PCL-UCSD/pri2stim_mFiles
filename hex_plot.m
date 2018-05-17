function ax = hex_plot(data,nrows,ax)
% TCS 10/24/14 - for plotting channel responses, weights, etc in hex grid

% TODO: add support for custom color limits, rotating points, etc

if nargin < 3
    figure;
    ax = axes;
end

ax.Color = [0.5 0.5 0.5];
hold on;


[xGrid, yGrid] = make_hex(nrows);

marker_size = 2/(2*nrows);

chan_range = [min(data) max(data)];
n_colors = 1000;
cmap = parula(n_colors+1);


xx = cos(linspace(0,2*pi,501))*marker_size;
yy = sin(linspace(0,2*pi,501))*marker_size;

for ii = 1:length(data)
    
    this_c = cmap(1+floor((data(ii)-chan_range(1))/(chan_range(2)-chan_range(1))*n_colors),:);
    
    patch(xx+xGrid(ii),yy+yGrid(ii),this_c);
    
    
    clear this_c;
    
end

axis equal




return