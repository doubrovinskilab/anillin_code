function res = do_intersect(a,b)

res = 0;

for i=1:1:length(a)
    if( sum(a(i)==b)>0 )
        res = 1;
        
        break;
    end
end

end