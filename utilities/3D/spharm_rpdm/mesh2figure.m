function mesh2figure(vs,fs,figtitle,filename,dpi)
% generate (and save) 3D figure from mesh
% vs = vertices
% fs = faces
% figtitle = (optional) title to put on figure
% filename = (optional) filenname to save figure to (default not saved)
% dpi = (optional) dots per inch % figtitle = (optional) title to put on figure (default none)
    figure, patch('vertices', vs, 'faces', fs, 'FaceVertexCData',jet(size(vs,1)),'FaceColor','interp');
    view([45, 45]);
    if exist('figtitle','var')
        title(figtitle);
    end
    if exist('filename','var') && ~isempty(filename)
        if ~exist('dpi','var') || isempty(dpi) dpi = 150; end
        export_fig(filename, '-opengl', '-png', '-a1', ['-r', num2str(dpi)]);
    end
    close
end
