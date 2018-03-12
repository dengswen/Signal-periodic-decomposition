function N = calEuler(q)
    %% 1) factor q base on the fundamental theorem of arithmetic
    
    if q ==1
        N = 1;
        return;
    end
    
    f = factor(q);
    L = length(f);
    s = sparse(q,1,0);
    for i =1:length(f)
        p = f(i);
        s(p) = s(p)+1;
    end
    ps = find(s>0);
   
    for i = 1:length(ps)
        p = ps(i);
        a= s(p,1);
        Ns(i)=calEuler_pa(p, a);
    end
    N = prod(Ns);
end
function N = calEuler_pa(p, a)
    if a ==1 
        N= p-1;
    else
        N=p^a*(1-1/p);
    end
end




