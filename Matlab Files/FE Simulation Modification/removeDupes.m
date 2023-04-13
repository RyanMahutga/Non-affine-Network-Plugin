function [ fibers ] = removeDupes( fibers )
%removeDupes removes duplicate fibers from the periodic network

unique_fibs = unique(fibers ,'rows') ;
    
    fibers =[];
    fibers=unique_fibs;
    
    m=0;
    while m < length(fibers)
        m=m+1;
        opposite = [fibers(m,2),fibers(m,1),-fibers(m,3:5)];
        n=m;
        while n < length(fibers)
            n=n+1;
            if fibers(n,:)==opposite
                fibers(n,:)=[];
            end
        end
        
        if fibers(m,1) == fibers(m,2)
            fibers(m,:)=[];
        end
    end

end

