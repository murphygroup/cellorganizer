function header2html ( fileID, heading)

    fprintf( fileID, '\n' );
    fprintf( fileID, '<h2>' );
    fprintf( fileID, heading );
    fprintf( fileID, '</h2>' );
    fprintf( fileID, '\n' );


end