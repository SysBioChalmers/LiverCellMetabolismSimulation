function a = fitMetmodel(t,S,tX,X)
    Xregress = zeros(length(t),1);
    for i = 1:length(tX)
        Xregress(t==tX(i)) = X(i);
    end
    a = fitlm(Xregress, S);
end

