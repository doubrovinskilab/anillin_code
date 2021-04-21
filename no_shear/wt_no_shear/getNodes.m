function SetOfNodes = getNodes(list1,list2)

SetOfNodes = [];

for i=1:1:length(list1)-1
    if( list1(i)==list1(i+1) | list1(i)==list2(i+1) )
        SetOfNodes = [SetOfNodes list2(i)];
    else
        SetOfNodes = [SetOfNodes list1(i)];
    end
end

if(  list1(length(list1))==list1(length(list1)-1) | list1(length(list1))==list2(length(list1)-1) )
    SetOfNodes = [SetOfNodes list1(length(list1))];
    SetOfNodes = [SetOfNodes list2(length(list2))];
else
    SetOfNodes = [SetOfNodes list2(length(list2))];
    SetOfNodes = [SetOfNodes list1(length(list1))];
end

end