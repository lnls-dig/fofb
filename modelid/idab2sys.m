function sys = idab2sys(a,b)

for i=1:size(a,1)
    sys{i} = tf(b(i,:),a(i,:),-1);
end