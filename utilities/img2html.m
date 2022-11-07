function img2html ( fileID, imagepath, imagepatht, caption)
    
if isdeployed
    
    fprintf ( fileID, '\n');

else
    
    fprintf( fileID, '\n' );
    fprintf( fileID, '<figure>\n' );
    fprintf( fileID, ['<a href="' imagepath '"><img src="' imagepatht '" /></a>\n'] );
    fprintf( fileID, ['<figcaption>' caption '</figcaption>\n'] );
    fprintf( fileID, '</figure>\n' );
    fprintf( fileID, '\n' );

end