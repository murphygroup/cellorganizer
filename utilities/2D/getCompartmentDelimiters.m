function delimiters = getCompartmentDelimiters( data )

delimiters = [];
for i=1:1:length( data )
    if( strfind( data{i}, 'compartment' ) )
        delimiters = [ delimiters, i ];
    end
end
end%getCompartmentDelimiters