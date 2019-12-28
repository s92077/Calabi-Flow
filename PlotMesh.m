%% PlotMesh
%  Display the triangular mesh.
%
%% Syntax
%   PlotMesh(F, V, Vrgb)
%
%% Description
%  F    : double array, nf x 3, faces of mesh
%  V    : double array, nv x 3, vertices of mesh
%  Vrgb : double array, nv x 3, RGB color of vertices (optional)
%
%% Contribution
%  Author : Mei-Heng Yueh
%  Created: 2016/09/06
% 
%  Copyright 2016 Mei-Heng Yueh
%  http://scholar.harvard.edu/yueh

function P = PlotMesh(F, V, Vrgb)
if nargin == 2
    e = ones(size(V,1),1);
    switch size(V,2)
        case 2
            P = patch('Faces', F, 'Vertices', V, 'FaceVertexCData', zeros(size(V,1),3), 'EdgeColor', 'interp', 'FaceColor', [0.6, 0.8, 1]);
        case 3
            Vrgb = [176/255*e, 224/255*e, 230/255*e];
            P = patch('Faces', F, 'Vertices', V, 'FaceVertexCData', Vrgb, 'EdgeColor', 'interp', 'FaceColor', 'interp');
            shading interp
    end
    light('Position',[0 1 1]);
else
    P = patch('Faces', F, 'Vertices', V, 'FaceVertexCData', Vrgb, 'EdgeColor','interp','FaceColor','interp');
end
set(gcf, 'color', [0 0 0]);
view([0, 0, 1]);
axis equal off

