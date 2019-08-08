function [charact, F] = resetpop(charact, x, F, edg, pop, feat)

% Reset connections 

F(x, :) = 0;
F( :,x) = 0;
                        % total
                        y = length(x);
                       % group, crime, power 12 13
                      charact(x,[1 5 9 12 13]) = 0;
                       charact(x,2) = randn(y,1);
                       
                       % find married
                       a2 = [ones(round(feat(1)*y),1); -1.*ones(y-round(feat(1)*y),1)];
                       charact(x,3) = a2(randperm(y));
                       
                       % wealth
                       charact(x,4) = gamrnd(feat(2),1,y,1)./2 - 0.5;
                       
                       %  % exposure
                       charact(x,6) = feat(9);
                       
                       
                       % Education
                       a2 = [-2.*ones(round(feat(4)*y),1); -1.*ones(round(feat(5)*y),1); zeros(round(feat(6)*y),1); ones(round(feat(7)*y),1); 2.*ones(y-round(feat(4)*y)- round(feat(5)*y) - round(feat(6)*y) -round(feat(7)*y),1)];
                       charact(x,7) = a2(randperm(y));
                       
                       % Religion
                       a2 = [zeros(round(feat(3)*n),1); randi([1 ,2 ], n-round(feat(3)*n) , 1)];
                       charact(x,8) = a2(randperm(y));
                       
                       %  % Bias Terrorist
                       charact(x,10) = 0.75;
                       
                       %  % Bias Military
                       charact(x,11) = 0.5;
                       
                       % Creating new edges for 
                       
                       edg = y*round(edg/pop) + y;
                       c = zeros(4.*edg, 2);
                       c(:,2)   = randi([1, pop],  4.*edg, 1);
                       c(:,1)  = randi([1, y],  4.*edg, 1);
                       c(:,1) = x(c(:,1));
                    c = sort(c,2); 
                    c = unique(c, 'rows');
                    c(c(:,1) - c(:,2) ==0, :) = [];
                    c = c(randperm(length(c), edg),:);

% These are the weights of the connection, which actual barely get used
c2 = [randi([10, 49], round(edg), 1); randi([51, 120], round(5*edg/2 + 1), 1 ); randi([ 121, 199], 5*edg- round(5*edg/2 + 1) -round(edg), 1)];

c2 = c2(randperm(length(c2)));

for i = 1:edg
    F(c(i,1) , c(i, 2)) = c2(i);
    F(c(i, 2) , c(i, 1)) = c2(i);
end

for i = 1:y
    a = find(F(x(i),:) >0, 1 );
    if( ~isempty(a))
        b = find(F(a,:) >0);
        nn = min(length(b), length(c2));
    F(x(i) , b(1:nn)) = c2(randperm(length(c2),nn));
    F(b(1:nn) , x(i)) = c2(randperm(length(c2),nn));
    F(x(i), x(i)) = 0;
    else
        b = randi([1,pop], 2,1);
     [~ , b2] = find( F(b,:) > 0);
       b = unique([b; b2]);
       
        nn = min(length(b), length(c2));
    F(x(i) , b(1:nn)) = c2(randperm(length(c2),nn));
    F(b(1:nn) , x(i)) = c2(randperm(length(c2),nn));
    end
    
end

if(length(F) ~= pop)
    msg = 'Error occurred.';
error(msg)
end

end