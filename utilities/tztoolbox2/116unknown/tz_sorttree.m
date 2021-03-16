function newtree=tz_sorttree(tree)

%
if(size(tree,1)<=1)
    newtree=tree;
    return;
end

curnode{1}=tree(1,1);

