function Normal = normalization(a);

l= length(a)
sum = 0;
for i=1:l
  sum = sum + a(i);
end

if sum == 0
    error ('The sum of the raw is zero')
end

for i=1:l 
  a(i) = a(i)/sum;
end
  
Normal = double(a);
end
  