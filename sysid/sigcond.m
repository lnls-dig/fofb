function y = sigcond(x, seglen, iir, fir, ndiscard)

if isempty(seglen)
    seglen = size(x,1);
end

if ~isempty(iir)
    x = filter(iir, x);
end

if ~isempty(fir)
    x = filter(fir, x);
end

if isempty(ndiscard)
    ndiscard = [0 0];
elseif length(ndiscard) < 2
    ndiscard(2) = 0;
end

x = x(ndiscard(1)+1:end-ndiscard(2), :);

% Ensure the number of samples is a multiple of the segment length
x = x(1:floor(size(x,1)/seglen)*seglen, :);

y = zeros(seglen, size(x,2));
for i=1:size(x,2)
    x_stacked = reshape(x(:,i), seglen, size(x,1)/seglen);
    y(:,i) = mean(x_stacked, 2);
end