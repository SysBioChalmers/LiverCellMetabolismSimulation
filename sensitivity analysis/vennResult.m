function  vennResult(result, tresh)
    lbMatch = abs((result(:,3)./result(:,2))-1)<tresh;
    ubMatch = abs((result(:,4)./result(:,2))-1)<tresh;
    
    %Adjust for negative flux, i.e. lb <-> ub for negative reactions
    for i = 1:length(lbMatch)
       if sign(result(:,2)) == -1
            tmp = lbMatch(i);
            lbMatch(i) = ubMatch(i);
            ubMatch(i) = tmp;
       end
    end
    
    
    both = and(lbMatch, ubMatch);
    either = or(lbMatch, ubMatch);
    unconstrained = or(result(:,3)>900, result(:,2)<-900);
    l=sum(lbMatch);
    u=sum(ubMatch);
    b=sum(both);
    e=sum(either);
    a=length(result);
    un = sum(unconstrained);
    fprintf('Det\t%i\n', b);
    fprintf('Lb\t%i\n', l-b);
    fprintf('Ub\t%i\n', u-b);
    fprintf('Con\t%i\n', a-un-e);
    fprintf('Un\t%i\n', un);
    fprintf('Tot\t%i\n', a);
    
    
    
end

