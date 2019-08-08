function [A2, cha, F, n] = popDem(A)

% This is Country Set up, Run this with specific country information


% This sets the seed, so that we can reproduce, the same results
s = RandStream('mt19937ar','Seed',210, 'NormalTransform', 'Ziggurat');
RandStream.setGlobalStream(s)



% Creating Terrorist Network
%There is only a Leader and Financer in the begining, but you can add more
lt = 2;
n2=A(10);

% 1: Financer, Leader
T = [  0 201 zeros(1,n2); 201  0 200.*ones(1,n2)];

% Foot soldiers, which are not existant
T = [T; zeros(n2,lt), 50.*ones(n2,n2)];

T = T- diag(diag(T));
% Number of civilian/military + terrorist households
n = 4000 +n2 +lt;



% Famility network connections
F = zeros(n,n);

F = graph(F);

% Number of connections (neighbor, friend)
n3 = 2.*n + 5;
c = randi([1, n],  3.*n3, 2);
c = [c; 1,2];
c = sort(c,2); 
c = unique(c, 'rows');
c(c(:,1) - c(:,2) ==0, :) = [];

c = c(randperm(length(c), n3),:)';

% These are the weights of the connection, which actual barely get used
c2 = [randi([10, 49],1, round(n3/5)), randi([51, 120],1, round(n3/2 + 1) ), randi([ 121, 199],1, n3- round(n3/2 + 1) -round(n3/5))];

c2 = c2(randperm(length(c2)));

F = addedge(F,c(1,:), c(2,:), c2 );
close all

c3 = find(degree(F)>1 & degree(F) < 6);

%
c = randi([1, length(c3)],  1, round(0.1*n));



for i = 1:length(c)
c1 = neighbors(F,c(i));
for j = 1: length(c1)
 F = addedge(F,c1([1:(j-1), j+1:end]), c1(j), c2(randperm(length(c1)-1)) );   
end
F = simplify(F);

end




% % to make sure everyone has at least two connection


ff = find(degree(F)==0 |degree(F)==1);
if( length(ff)/5 > 7)
    for i = 1: ceil(length(ff)/5)
        ff2 = ff(i:5:end);
        for j = 1:length(ff2)
    F = addedge(F, ff2([1:(j-1), j+1:end]), ff2(j), c2(1:length(ff2)-1));
        end
    end
else
    for i = 1: ceil(length(ff)/4)
        ff2 = ff(i:4:end);
    
        for j = 1:length(ff2)
    F = addedge(F, ff2([1:(j-1), j+1:end]) , ff2(j), c2(1:length(ff2)-1));
        end
    end
end



F = simplify(F);




% characteristic matrix
cha = zeros(n, 13);

% Group identity Civilian: 0, Terrorist leader: 4 Terrorist financer: 3
% Terrorist foot soilder:2 military/police: 1
cha(randperm(n -n2, round(.05*(n-n2)+2)) +n2, 1) = 1;

cha(1,1) = 3;
cha(2,1) = 4;
cha((lt+1):(lt+n2),1) = 2;

% personal achievement
cha(:,2) = randn(n,1);
cha(2,2) = sum(maxk(cha(:,2),4))/4;
cha(1,2) = cha(2,2);


% Not married -1, married 1
% UN data
z = [ ones(1, round(A(1)*n)) , -1.*ones(1, n - round(A(1)*n))];
cha(:,3) = z(randperm(n));

cha(1:2,3) = 1;




% wealth
% 1/3 below poverty line CIA factbook
% I attempt to assign it so that neighbors are more likely to be of
% similiar wealth
z = gamrnd(A(2), 1, n+2,1)./2 - 0.5;


z = sort(z, 'descend');
p = 1:n;

z2 = neighbors(F,1);
% 
z2l = round(length(z2)/2 +1);
cha(z2(2:(z2l)), 4) = z(1:(z2l-1));
cha(1,4)=  z(z2l);

p = setdiff(p,[z2(2:z2l); 1]);
%
for i= 2:z2l
    
    z3 = intersect(neighbors(F,z2(i)),p);
    
z3l = round(length(z3)/2);
cha(z3(1:z3l), 4) = z((z2l+1):(z2l+z3l));
z2l = z2l + z3l;
p = setdiff(p,z3(1:z3l));
end
%
%
p2 = find(cha(:,4) == 0);
p2 = p2(randperm(length(p2),1));

while(z2l < n)
    
z3 = neighbors(F,p2);
while( ~isempty(p) &&  ~isempty(z3))
z2 = z3;
for i= 1:length(z2)
    
    z3 = intersect(neighbors(F,z2(i)),p);
    
z3l = length(z3);
cha(z3, 4) = z((z2l+1):(z2l+z3l));
z2l = z2l + z3l;

p = setdiff(p,z3);
end

end

p2 = find(cha(:,4) == 0);
if(~isempty(p2))
p2 = p2(randperm(length(p2),1));
cha(p2,4) = z((z2l+1));
p = setdiff(p,p2);
z2l = z2l +1;
end

end

%
%

% Commited Criminal event

cha(randi([1 n],  5, 1), 5) = randi([1 3], 5,1);

% Exposure to violent event

cha(:,6) = A(9);

% Education  8% illiterate
% 29.3 % tertiary schooling, 29.3 % secondary (15+), 25.2 % primary (15+) 
% 22.1 % no education (15+)
% illiterate: -2, low: -1, medium: 0, high school completed/some college: 1,
% bachelors: 2 more than bachelors: 3

m1 = min(cha(:,4));
ran = max(cha(:,4)) -m1;
z2 = find(cha(:,4) < (0.2*ran + m1));
z2 = z2(randperm(length(z2), round(A(4)*n)));
cha(z2,7) = -2;

z2 = find(cha(:,4) < (0.4*ran + m1) & cha(:,7) == 0  );
z2 = z2(randperm(length(z2), round(A(5)*n)));
cha(z2,7) = -1;



z2 = find(cha(:,4) < (0.75*ran + m1)  & cha(:,7) == 0   );
z2 = z2(randperm(length(z2), round(A(6)*n)));
cha(z2,7) = 0;

z2 = find( cha(:,7) == 0   );
z2 = z2(randperm(length(z2), round(A(7)*n)));
cha(z2,7) = 1;


cha(cha(:,7) == 0 ,7) = 2;




% This is when someone was religious, 1 at later age, 2 at early age, 0
% not muslim 

z = randi([1 2], n- round(A(3)*n),1);
if(A(3) > 0.9)
l = length(z);
z2 =  neighbors(F,2);
z2 = setdiff(z2,[1,2]);
l2 = length(z2);
if((l2 +4) <= l)

    cha(z2,8) = z(1:l2);
    p = randperm(n);
    p = setdiff(p,z2);
    p = setdiff(p, find(cha(:,4) > median(cha(:,4))));
    l3 = randperm(length(p), l - l2  );
    l3 = p(l3);
    cha(l3,8) = z(l2+1:end);
    
    l5 = round(length(l3)/2);
    F = addedge(F, l3(1:l5), 2, c2(1:l5));
    
     for jy = 1:l5
    F = addedge(F, l3(1:l5), l3(jy), c2(1:l5));
     F = rmedge(F, l3(jy), l3(jy));
     end
else
    cha(z2(2:l-3),8) = z(1:(l-4));
    
     p = randperm(n);
    p = setdiff(p,z2);
    p = setdiff(p, cha(:,4) > median(cha(:,4)));
    l3 = randperm(length(p), 4 );
    l3 = p(l3);
    cha(l3,8) = z((l-3):end);
    F = addedge(F, l3, 2, c2(1:length(l3)));
    
     for jy = 1:length(l3)
    F = addedge(F, l3, l3(jy), c2(1:length(l3)));
     F = rmedge(F, l3(jy), l3(jy));
     end
end

else
    z = [ z; zeros(round(A(3)*n) ,1)];
cha(:,8) = z(randperm(n))';
end
% Power is 9


% Terrorist bias
cha(:,10) = 0.75;

% Military bias
cha(:,11) =0.5;




% The goal of the cell leader (2nd person, type 4) is to recruit others (type 2), 
% then when there is enough the attack is planned

% The goal of the military (type 1) is to get terrorist

% financier wants enough money to support terrorist efforts.

% All connections
A2 = adjacency(simplify(F), 'Weighted');
A2(1:(lt+n2),1:(lt+n2)) =  T;
%Ag = digraph(A);

% Then a randomly select its out going connection
%



end