function Xregress = getXestimate(t,tX,X)
    Xregress = zeros(length(t),1);
    for i = 1:length(tX)
        Xregress(t==tX(i)) = X(i);
    end
end

