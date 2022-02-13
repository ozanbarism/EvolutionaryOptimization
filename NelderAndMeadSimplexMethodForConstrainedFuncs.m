clear
clc

%STEP 1
% initial H and D estimates ;
x(1,:)=[2.5,1,0.5, 0.5];
x(2,:)=[1.5,0.5,0.3, 0.7];
x(3,:)=[3.5, 0.5, 0.25, 0.1];
x(4,:)=[1.5,2,0.1,0.6];
x(5,:)=[2,1.5,0.25,0.25];

%If the user wants to use a different initial estimate randomly chosen
%within the boundries specified, It should define yes here. If not, this
%can be left empty.
RandomInitialEstimates="NO";

if RandomInitialEstimates=="YES"
%for different initial estimates
LB=[1 0.25 0.01 0.01];
UB=[5 2.5 0.1 0.1];
%Bounds for very large starting points
LB=[20 10 5 5];
UB=[40 20 10 10];

 for i=1:4
     
     x(:,i)=rand(5,1).*(UB(i)-LB(i))+LB(i);
 end
 x=round(x,2);  
x_initial=x;
end

%parameters such as reflection, contraction, expansion, scaling
r=1;
c=0.5;
e=0.5;
s=0.5;
eps = 0.01;
%penalty parameter is defined
r_pen=10^2;
r_pen_increase=1;


%values of those points are computed
f_i(1)=FuncValue(x(1,:),r_pen);
f_i(2)=FuncValue(x(2,:),r_pen);
f_i(3)=FuncValue(x(3,:),r_pen);   
f_i(4)=FuncValue(x(4,:),r_pen);   
f_i(5)=FuncValue(x(5,:),r_pen); 
f_maxes_i=maxk(f_i,5);
%STEP 2
%the correct xh xl xs are defined here.
x_h=x(find(max(f_i)==f_i),:);
x_2h=x(find(f_maxes_i(2)==f_i),:);
x_3h=x(find(f_maxes_i(3)==f_i),:);
x_4h=x(find(f_maxes_i(4)==f_i),:);
x_l=x(find(min(f_i)==f_i),:);


    

  
    
iter=0;
while true
    %iteration increased one more
iter=iter+1;    



f_h=FuncValue(x_h,r_pen);
f_2h=FuncValue(x_2h,r_pen);
f_3h=FuncValue(x_3h,r_pen);
f_4h=FuncValue(x_4h,r_pen);
f_l=FuncValue(x_l,r_pen);

%STEP 3
%center point of the points excluding x_h
x_o=(x_l+x_2h+x_3h+x_4h)/4;
f_o=FuncValue(x_o,r_pen);
%STEP 4
%reflection is conducted and its value is computed
x_r=x_o+r*(x_o-x_h);
f_r=FuncValue(x_r,r_pen);



%STEP 5 COMPARISIONS ARE DONE
if (f_r < f_l) && (f_r>0) && (f_l>0)
    x_e=x_r+e*(x_r-x_o);
    
    f_e=FuncValue(x_e,r_pen);
    
    if f_e<f_r
     
        x_h=x_e;
    else 
    
        x_h=x_r;
        
    end
    %step 7
    std=sigma(f_h,f_2h,f_3h,f_4h,f_l,f_o);
    if std <= eps
        break
    end
    
elseif (f_r > f_h)&&(f_r>0)&&(f_h >0)
    
    x_c=x_o-c*(x_o-x_h);
    
    f_c=FuncValue(x_c,r_pen);
    %step 6
      if f_c < f_h
      
        x_h=x_c;
 
      else
    
        x_h=x_l+s*(x_h-x_l);
        x_2h=x_l+s*(x_2h-x_l);
        x_3h=x_l+s*(x_3h-x_l);  
        x_4h=x_l+s*(x_4h-x_l);
        
      end

      %step 7
            std=sigma(f_h,f_2h,f_3h,f_4h,f_l,f_o);
             if std <= eps
                 break
             end
      
      
      
end   

if (f_r < f_h)&&(f_r > f_2h)&&(f_r>0)&&(f_h >0)&&(f_2h >0)  
    
    x_c=x_r-c*(x_r-x_o);
    
    f_c=FuncValue(x_c,r_pen);
    %step 6
      if (f_c < f_h)&&(f_c>0)&&(f_h >0)
        x_h=x_c;
      else
        x_h=x_l+s*(x_h-x_l);
        x_2h=x_l+s*(x_2h-x_l);
        x_3h=x_l+s*(x_3h-x_l);  
        x_4h=x_l+s*(x_4h-x_l);
        
      end
      
      %step 7
            std=sigma(f_h,f_2h,f_3h,f_4h,f_l,f_o);
             if std <= eps
                 break
             end
      
      
end

if (f_l < f_r) && (f_r < f_2h) && (f_r>0) && (f_l >0) && (f_2h >0)    
    x_h=x_r;
            %step 7
            std=sigma(f_h,f_2h,f_3h,f_4h,f_l,f_o);
             if std <= eps
                 break
             end
end
  
            %step 7
          std=sigma(f_h,f_2h,f_3h,f_4h,f_l,f_o);
            if std <= eps
                break
            end

x_new(1,:)=x_h(1,:);
x_new(2,:)=x_2h(1,:);
x_new(3,:)=x_3h(1,:);
x_new(4,:)=x_4h(1,:);
x_new(5,:)=x_l(1,:);
%values of the functions are computed again for the accuracy of the next
%iterations
for i=1:5
f(i)=FuncValue(x_new(i,:),r_pen);
end
[val idx]=sort(f);

%since x_h will change in every iteration, the values should be corrected
%again.

x_h=x_new(idx(5),:);
x_2h=x_new(idx(4),:);
x_3h=x_new(idx(3),:);
x_4h=x_new(idx(2),:);
x_l=x_new(idx(1),:);    
    
r_pen=r_pen*r_pen_increase;  
end

x_opt=x_l;
%values are rounded to achieve logical results
x_opt=x_opt.*10^4;
x_opt=ceil(x_opt);
x_opt=x_opt./10^4;


%printing the optimum values and points
fprintf(['\n','optimum point is ' num2str(x_opt)]);

f_opt=FuncValue(x_opt,r_pen);
fprintf(['\n','optimum value is ' num2str(f_opt) ]);
    
f_wop = (x_opt(1)*x_opt(4)+2*x_opt(2)*x_opt(3))*35;
fprintf(['\n','optimum value without penalty is ' num2str(f_wop) ]);
fprintf(['\n','iteration number is ' num2str(iter) ]);


if round(f_opt,3)==round(f_wop,3)
    fprintf(['\n','all constraints are satisfied ']);
else
    fprintf(['\n','some constraints are not satisfied increase penalty value ']);
end

function std=sigma(f_h,f_2h,f_3h,f_4h,f_l,f_o);

sum=(f_h(1,1)-f_o(1,1))^2+(f_2h(1,1)-f_o(1,1))^2+(f_3h(1,1)-f_o(1,1))^2+(f_4h(1,1)-f_o(1,1))^2+(f_l(1,1)-f_o(1,1))^2;
std=sqrt(sum/5);
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
