function cmap = custom_colormap()

% Please do not change this code

% Creates an RGB colormap that starts with red ([1 0 0]), continues to
% white ([1 1 1]) and ends with blue ([0 0 1])

numColors = 100;
r = [ones(numColors/2,1); linspace(0.8,0,numColors/2-1)'];
g = [linspace(0,0.8,numColors/2-1)'; 1; linspace(0.8,0,numColors/2-1)'];
b = [linspace(0,0.8,numColors/2-1)'; ones(numColors/2,1)];
cmap = [r g b];

end