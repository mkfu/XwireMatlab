
N = 5;

dist  =40;
interval = 10;
in  = zeros(dist./interval,N);
out = zeros(dist./interval,N);
step = 256;
speed =2000*step;
for i = 1:N
    m.move(0.246);
    for j = 1:dist/interval
        [in(j,i),trav] = m.move(interval);
        trav
        
    end
    m.move(-0.246);
    for j = 1:dist/interval
        [out(j,i),trav] = m.move(-interval);
        trav
        
    end
    
end

