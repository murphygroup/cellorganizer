function text2html ( fileID, string)

    %fprintf( fileID, '\n' );
    fprintf( fileID, '<p>' );
    fprintf( fileID, string );
    fprintf( fileID, '</p>' );
    %fprintf( fileID, '\n' );

end