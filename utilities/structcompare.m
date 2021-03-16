function test = structcompare(struct1,struct2)
%STRUCTCOMPARE Helper function that compares two structures to see if they have the same fields.

%create cell array of both
array1 = struct2cell(struct1);
array2 = struct2cell(struct2);

if length(array1)~=length(array2)
    test = 0;
    fprintf('number of elements in struct did not match');
    array1
    array2
    return
end

test = looparray(array1,array2);

end


function test = looparray(array1,array2)

test = cell(1,length(array1));
for i = 1:length(array1)
    for j = 1:length(array2)
    %if it's still a struct recurse to top
    if isstruct(array1{i})&&isstruct(array2{j})
        test{i} = structcompare(array1{i},array2{j});
    elseif (~isstruct(array1{i})&&isstruct(array2{j}))||(isstruct(array1{i})&&~isstruct(array2{j}))
        test{i} = 0;
%         fprintf('struct elements of array did not match in location %s',num2str(i));
%         array1
%         array2
    elseif iscell(array1{i})&&iscell(array2{j})
        %if the field is an array recurse to looparray
        test{i} = looparray(array1{i},array2{j});
    elseif (~iscell(array1{i})&&iscell(array2{j}))||(iscell(array1{i})&&~iscell(array2{j}))
        test{i} = 0;
        fprintf('cell array elements of array did not match in location %s',num2str(i));
%         array1
%         array2
    elseif ischar(array1{i})&&ischar(array2{j})
        %if not a struct or cell, compare values
    %if strings
        test{i} = strcmp(array1{i},array2{j});
    elseif (~ischar(array1{i})&&ischar(array2{j})) || (ischar(array1{i})&&~ischar(array2{j}))
        test{i} = 0;
        fprintf('string elements of array did not match in location %s',num2str(i));
        array1
        array2
    elseif isnumeric(array1{i})&&isnumeric(array2{j})
        %number
        test{i} = all(array1{i}(:)==array2{j}(:));
    elseif (isnumeric(array1{i})&&isnumeric(array2{j})) || (isnumeric(array1{i})&&isnumeric(array2{j}))
        test{i} = 0;
        fprintf('numeric elements of array did not match in location %s',num2str(i));
        array1
        array2
    else
        test{i} = 0;
        fprintf('unrecognized type in arrays in location %s',num2str(i));
        array1
        array2
    end
    end
%     if j == length(array2)
%         
%     end
end

end
