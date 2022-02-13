clear
clc
%Whale Optimization Algorithm created by A. Lewis and S. Mirjalili (2016)
%is used to find the most optimum design for a bridge girder.

%This for loop is created to compute the algorithm ten times, if the user
%do not want to do this. running time can be entered as 1.
NumberOfRunnings=10; %the number of runnings conducted for algorithm


%boundries are defined
LB=[1 0.25 0.01 0.01];
UB=[5 2.5 0.1 0.1];

%parameters are defined
Max_iter=100;
search_agent=100;
%penalty value is defined
r_pen=10^5;
r_pen_increase=1.2;
%paremeter for the shape of the spiral
b=1;
%number of variables
numberofvariables=4;

for z=1:NumberOfRunnings
%initialization
x_star=zeros(1,numberofvariables);
f_star=inf;
collect(1)=f_star;


%whale population is initialized
 for i=1:numberofvariables  
     x(:,i)=rand(search_agent,1).*(UB(i)-LB(i))+LB(i);
 end

%collection array was defined
collect=zeros(1,Max_iter);

%iteration counter was defined
t=1;

while t <= Max_iter
    
    % loop is conducted for every search agent
    for i=1:search_agent
    
    %fitness is computed;
    fitness=FuncValue(x(i,:),r_pen);
    
    if fitness<f_star
        f_star=fitness;
        x_star=x(i,:);
    end    
        
    %a and a2 was defined in terms of max_iter and iteration counter;
    a=2-2*t/Max_iter;
    a2=-1+t*((-1)/Max_iter);
        %variables such as r2 r1 p b A C l are defined;
    r2=rand;
    r1=rand;
    p=rand;
    A=2*a*r1-a;
    C=2*r2;
    l=(a2-1)*rand+1;

    for j=1:numberofvariables
        if p< 0.5
            if abs(A)<1
                D=abs(C*x_star(j)-x(i,j));
                x(i,j)=x_star(j)-A*D;
            elseif abs(A)>=1
                 x_rand(j) = LB(j) + (UB(j)-LB(j)).*rand;    
                 D=abs(C*x_rand(j)-x(i,j));
                 x(i,j)=x_rand(j)-A*D;    
            end
        elseif p>=0.5
    
        D_prime=abs(x_star(j)-x(i,j));
        x(i,j)=D_prime*exp(b.*l)*cos(l.*2*pi)+x_star(j);

        end

     %Upper and lowers bounds of the variables are checked, and they are
     %corrected if they went beyond boundries
    while (x(i,j)>UB(j))||(x(i,j)<LB(j))
       x_rand(1,j) = LB(j) + (UB(j)-LB(j)).*rand;    
       x(i,j)=x_rand(1,j);  
    end  
        
        
        
        
        
end


    end
r_pen=r_pen*r_pen_increase;  
t=t+1;
collect(t)=f_star;

end
%the values in the collect array are rounded by 4 decimals to find when its
%converges
RoundedCollect=round(collect,4);
LastResult=round(collect(length(collect)),4);




%printing the optimum values and points
x_opt=x_star(1,:);
fprintf(['\n','optimum point is ' num2str(x_opt)]);

f_opt=FuncValue(x_opt,r_pen);

%optimum for each running is stored
optimum4eachrun(z,1)=f_opt;
%index of the collect where convergence occurred are found
index=find(RoundedCollect==LastResult);
%using this index convergence time for each running is stored
convergenceTime4eachrun(z,1)=index(1);
%the values computed are printed out
fprintf(['\n','optimum value is ' num2str(f_opt) ]);
    
f_wop = (x_opt(1)*x_opt(4)+2*x_opt(2)*x_opt(3))*35;
fprintf(['\n','optimum value without penalty is ' num2str(f_wop) ]);
if round(f_opt,4)==round(f_wop,4)
    fprintf(['\n','all constraints are satisfied ']);
    %the result is stored in satisfaction array.
    satisfaction(z,1)="Satisfied";
else
    fprintf(['\n','some constraints are not satisfied increase penalty value ']);
     %the result is stored in satisfaction array.
    satisfaction(z,1)="Not Satisfied";
end

end
%objective function
function fi = FuncValue(x,r)
    G = ConstraintVal(x);
    h=x(1);
    b=x(2);
    tf=x(3);
    tw=x(4);
    
    L=35;

    penalty=r*(max(0,G(1,1))^2+max(0,G(2,1))^2+max(0,G(3,1))^2+max(0,G(4,1))^2+max(0,G(5,1))^2+max(0,G(6,1))^2+max(0,G(7,1))^2+max(0,G(8,1))^2+max(0,G(9,1))^2+max(0,G(10,1))^2+max(0,G(11,1))^2+max(0,G(12,1))^2+max(0,G(13,1))^2+max(0,G(14,1))^2);

    fi = (h*tw+2*b*tf)*L+penalty;
 
 %normalized constraints
function G = ConstraintVal(x)
    
    h=x(1);
    b=x(2);
    tf=x(3);
    tw=x(4);
    
    
    L=35;
    Pm=104;
    Ps=155;
    E=200;
    A=h*tw+2*b*tf;
    I=1/12*tw*h^3+2/3*b*tf^3+1/2*b*tf*h*(h+2*tf);
    w=19+77*A;
    M=L/8*(2*Pm+w*L);
    sigma=M/(1000*I)*(0.5*h+tf);
    sigma_f=72845*(tf/b)^2;
    sigma_w=3648276*(tw/h)^2;
    S=0.5*(Ps+w*L);
    D=(L^3/(384*10^6*E*I))*(8*Pm+5*w*L);
    tao=S/(1000*h*tw);
    
    %limitations
    sigma_y=250;
    sigma_a=0.55*sigma_y;
    tao_a=0.33*sigma_y;
    sigma_t=243;
    D_a=L/800;
    
    G(1,1) = sigma/sigma_a-1;
    G(2,1) = sigma/sigma_f-1;
    G(3,1) = sigma/sigma_w-1;
    G(4,1)= tao/tao_a-1;
   G(5,1) = D/D_a-1;
    G(6,1)= sigma/(0.5*sigma_t)-1;
    G(7,1)= h/5-1;
    G(8,1)= 1-h/1;
    G(9,1)= b/2.5-1;
    G(10,1)= 1-b/0.25;
    G(11,1)= tf/0.10-1;
    G(12,1)=1-tf/0.01;
    G(13,1)= tw/0.10-1;
    G(14,1)=1-tw/0.01;

end
end