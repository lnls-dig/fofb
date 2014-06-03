function mimo = siso2mimo(siso,n)

if n < 1
    mimo = [];
elseif n > 1
    mimo = tf(num2cell(zeros(n)), num2cell(ones(n)));

    if length(siso) == 1
        siso_ = siso;
    else
        error('n must be equal to 1 or length of siso.');
    end
        
    for i=1:n
        if size(siso,1) == n
            siso_ = siso(i,1);
        elseif size(siso,2) == n
            siso_ = siso(1,i);
        end
        mimo(i,i) = siso_;
    end
else
    mimo = siso;
end