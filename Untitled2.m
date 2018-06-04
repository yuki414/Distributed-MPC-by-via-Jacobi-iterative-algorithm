A=[0 1;-1 0];
Ad = expm(A*0.1);
data = [];
x = [0,1]';
for t = 0:0.1:15
    x = Ad*x;
    data=[data;t,x'];
end
plot(data(:,1),data(:,2))