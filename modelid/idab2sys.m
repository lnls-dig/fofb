function sys = idab2sys(a, b, Ts)

if nargin < 3
    Ts = -1;
end

for i=1:length(a)
    sys{i} = tf(b{i}, a{i}, Ts);
end