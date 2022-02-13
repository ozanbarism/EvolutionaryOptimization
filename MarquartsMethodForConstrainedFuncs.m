clear
clc

r=600; %penalty value
r_update=1.0;

%loop is created to compute the algorithm ten times, if the user
%do not want to do this. running time can be entered as 1.
NumberOfRunnings=1; %the number of runnings conducted for algorithm

%If the user wants to use a different initial estimate randomly chosen
%within the boundries specified, It should define yes here. If not, this
%can be left empty.
RandomInitialEstimates="NO";

M =500; %maximum number of iterations
x = [2.5 1 0.5 0.5]; %initial x
e = 10^-3; %allowed error
lamda = 10^3;

eps_xi=10^-7; %eps/xi ratio

for z=1:NumberOfRunnings
    
    if RandomInitialEstimates=="YES"
%a different set of initial estimates were conducted
LB=[1 0.25 0.01 0.01];
UB=[5 2.5 0.1 0.1];
x(1,:)=rand.*(UB(:)-LB(:))+LB(:);
x_initial=x(1,:);
    end
iter=0; %iteration counter
k = 1; % there are cases where k is not increased and loop is conducted one more, so another iteration counter is defined

%a while loop is conducted for until maximum number of iterations are
%reached

while k <= M
    %iteration increased one more
    iter=iter+1;
   
    

    %epsilon is defined as the 0.0001 times of the corresponding x value
    eps = eps_xi*x(k,:);    
    %a for loop is conducted to compute gradient
    for i=1:size(x,2)
        x_up=x(k,:);
        x_down=x(k,:);
        %upper function defined where +0.5eps added
        x_up(i)=x(k,i)+0.5*eps(i);
        %lower function defined where -0.5eps added
        x_down(i)=x(k,i)-0.5*eps(i);
        df_dx(i,1)=(FuncValue(x_up,r)-FuncValue(x_down,r))/eps(i);
        %variables are turned back to what they were at the beginning for
        %the next loop
        x_up=x(k,:);
        x_down=x(k,:);
    end
    
    %another for loop for hesian matrix
    for i=1:size(x,2)
        for j=1:size(x,2)
            %if same variable is derived twice, eps added twice for the
            %first term
            if i==j
            x_i_up=x(k,:);
            x_i_up(i)=x(k,i)+2*eps(i);
            x_i_up_2=x(k,:);
            x_i_up_2(j)=x(k,j)+eps(j);
            
            x_i_down=x(k,:);
            x_i_down(i)=x(k,i)+eps(i);
            x_i_down_2=x(k,:);
            else
                %for other cases, to go up eps added, and in the minus
                %term, the x leaved as it was before.
            x_i_up=x(k,:);
            x_i_up(i)=x(k,i)+eps(i);
            x_i_up(j)=x(k,j)+eps(j);
            x_i_up_2=x(k,:);
            x_i_up_2(j)=x(k,j)+eps(j);
            
            x_i_down=x(k,:);
            x_i_down(i)=x(k,i)+eps(i);
            x_i_down_2=x(k,:);
            end
            %computing the first and second term, the hessian matrix are
            %created
            first_term=((FuncValue(x_i_up,r)-FuncValue(x_i_up_2,r))/eps(i));
            second_term=((FuncValue(x_i_down,r)-FuncValue(x_i_down_2,r))/eps(i));
            df2dx2(i,j)=((first_term-second_term)/eps(j));

            %df2dx2(i,j)=(((FuncValue(x_i_up,r)-FuncValue(x_i_up_2,r))/eps(i))-((FuncValue(x_i_down,r)-FuncValue(x_i_down_2,r))/eps(i)))/(eps(j));
            
        end
    end
    %condition to break the iteration are introduced
     if (norm(df_dx) < e)||(k >= M) 
        break
     end  
   
    %the value to increase the x for the next iteration is computed.
    I =eye(4);
    sk = -inv(df2dx2+lamda*I)*df_dx;
    
   
 
    
    % defining the next point
    for i=1:size(x,2)
    x(k+1,i) = x(k,i)+sk(i);
    end
   
    %current point and next iteration point are computed in objective func.
    fk0 = FuncValue(x(k,:),r);
    fk1 = FuncValue(x(k+1,:),r);
   
    %checking the condition 
    
    %if the new point is a better solution, lamda divided by two and k is
    %increased one more.
    if fk1 < fk0
        lamda = lamda/2;
        k = k + 1;  
    else
        %if not lamda is increased twice, and the same iteration is
        %conducted again.
        lamda = 2*lamda;  

    end
    %r is updated
     r=r*r_update;
end

%printing the optimum values and points
x_opt=x(length(x),:);


%optimum value is computed
f_opt=FuncValue(x_opt,r);
%value of the optimum point without penalty function is computed
f_wop = (x_opt(1)*x_opt(4)+2*x_opt(2)*x_opt(3))*35;
Iter4eachrun(z,1)=k;
opt4eachrun(z,1)=f_opt;
%values are printed out
fprintf(['\n','optimum point is ' num2str(x_opt)]);
fprintf(['\n','optimum value is ' num2str(f_opt) ]);
fprintf(['\n','optimum value without penalty is ' num2str(f_wop) ]);
fprintf(['\n','iteration number is ' num2str(k) ]);
%a if case is conducted to see if all constraints are satisfied by checking
%the f values for objective funciton and objective function without the
%penalty function.
x_opt=round(x_opt,4);
if round(f_opt,3)==round(f_wop,3)
    fprintf(['\n','all constraints are satisfied ']);
     satisfaction(z,1)="Satisfied";
else
    fprintf(['\n','some constraints are not satisfied increase penalty value ']);
    satisfaction(z,1)="Not Satisfied";
end
if length(x)==M
    fprintf(['\n','max number of iterations reached']);
end

end


%objective function
function fi = FuncValue(x,r)
    G = ConstraintVal(x);
    L=35;
    h=x(1);
    b=x(2);
    tf=x(3);
    tw=x(4);

    penalty=r*(max(0,G(1,1))^2+max(0,G(2,1))^2+max(0,G(3,1))^2+max(0,G(4,1))^2+max(0,G(5,1))^2+max(0,G(6,1))^2+max(0,G(7,1))^2+max(0,G(8,1))^2+max(0,G(9,1))^2+max(0,G(10,1))^2+max(0,G(11,1))^2+max(0,G(12,1))^2+max(0,G(13,1))^2+max(0,G(14,1))^2);

    fi = (h*tw+2*b*tf)*L+penalty;
 
 %normalized constraints
function G = ConstraintVal(x)
    
    h=x(1);
    b=x(2);
    tf=x(3);
    tw=x(4);

    
    %parameters are defined, and computed
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
    %constraint values are computed
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
