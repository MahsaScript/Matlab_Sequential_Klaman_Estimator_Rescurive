function [stiff,damp,fkp,fap,fbp]=kcm(n,xk)
x1=xk(1:10);                                  
x2=xk(11:20);                                 
x3=xk(21:31);                                 
x4=xk(32);                                    
x5=xk(33);                                    

load mass.mat;
a=[ 2^0.5 3*pi/4 1 2 0 0;                     
    2^0.5 pi/4   1 2 7 8;
    2   0  0 0 7 8;
    2   0  1 2 3 4;
    2^0.5 3*pi/4 3 4 7 8;
    2   0  7 8 9 10
    2   0  3 4 5 6
    2^0.5 pi/4 3 4 9 10
    2^0.5 3*pi/4 5 6 9 10
    2   0  9 10 0 0
    2^0.5 pi/4 5 6 0 0];

ke(1:11)=x3;
K=zeros(10);
for i=1:11
    p=a(i,:);
    s=a(:,2);
    T=[cos(s(i)) sin(s(i)) 0 0;
        -sin(s(i)) cos(s(i)) 0 0;
        0 0 cos(s(i)) sin(s(i));
        0 0  -sin(s(i)) cos(s(i))];
    k=ke(i)*[1 0 -1 0;
             0 0  0 0;
             -1 0 1 0;
             0  0 0 0];
  smd=T'*k*T;
          for j=1:4          
            for k=1:4
               if p(2+j)~=0 
                 if p(2+k)~=0
                    ii=p(2+j);
                    jj=p(2+k);
                    K(ii,jj)=K(ii,jj)+smd(j,k);
                  
                end
                  end
              end
          end  
end         
stiff=K;

%%%%%%%%%%%%%


damp=mass*x4+stiff*x5;


% % 


for i=1:11
    K1=zeros(10);
    p=a(i,:);
    s=a(:,2);
    T=[  cos(s(i))   sin(s(i))       0         0;
        -sin(s(i))   cos(s(i))       0         0;
              0         0       cos(s(i))    sin(s(i));
              0         0      -sin(s(i))    cos(s(i))];
          k=[1 0 -1 0;
             0 0  0 0;
             -1 0 1 0;
             0  0 0 0];
  smd=T'*k*T;
          for j=1:4          
            for k=1:4
               if p(2+j)~=0 
                 if p(2+k)~=0
                    ii=p(2+j);
                    jj=p(2+k);
                    K1(ii,jj)=K1(ii,jj)+smd(j,k);
                  
                end
                  end
              end
          end  
          fkp(:,i)=K1*x1+K1*x2*x5;
end         

%%%%%%%%%%%%%%%%%%%%%%%%%

fap=mass*x2;
fbp=stiff*x2;



% load s1.txt;load s2.txt;load s3.txt;load s4.txt;load s5.txt;load s6.txt;load s7.txt;load s8.txt;load s9.txt;load s10.txt;load s11.txt;
% x=xk(1:n); fkp1=zeros(n,11);
% fkp1(:,1)=s1*x;fkp1(:,2)=s2*x;fkp1(:,3)=s3*x;
% fkp1(:,4)=s4*x;fkp1(:,5)=s5*x;fkp1(:,6)=s6*x;
% fkp1(:,7)=s7*x;fkp1(:,8)=s8*x;fkp1(:,9)=s9*x;
% fkp1(:,10)=s10*x;fkp1(:,11)=s11*x;

% damp=x4*mass+x5*stiff;
% xp=xk(n+1:2*n); 
% fcp=zeros(n,11);
% fcp(:,1)=s1*xp;fcp(:,2)=s2*xp;fcp(:,3)=s3*xp;
% fcp(:,4)=s4*xp;fcp(:,5)=s5*xp;fcp(:,6)=s6*xp;
% fcp(:,7)=s7*xp;fcp(:,8)=s8*xp;fcp(:,9)=s9*xp;
% fcp(:,10)=s10*xp;fcp(:,11)=s11*xp;
% 
% fkp=fkp1+fcp*x5;