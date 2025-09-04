n=11;
m=10;
A=7.854e-5*ones(1,11);
E=2e7;
RO=7850;
A(8)=A(8);

a=[ 2^0.5 3*pi/4 1 2 0 0;                      
    2^0.5 pi/4   1 2 7 8;
    2   0  0 0 7 8;
    2   0  1 2 3 4;
    2^0.5 3*pi/4 3 4 7 8;
    2   0  7 8 9 10;
    2   0  3 4 5 6;
    2^0.5 pi/4 3 4 9 10;
    2^0.5 3*pi/4 5 6 9 10;
    2   0  9 10 0 0;
    2^0.5 pi/4 5 6 0 0];


K=zeros(m,m);
  for i=1:n
      
    p=a(i,:);
    l=a(:,1);
    s=a(:,2);
    T=[cos(s(i)) sin(s(i)) 0 0;
        -sin(s(i)) cos(s(i)) 0 0;
        0 0 cos(s(i)) sin(s(i));
        0 0  -sin(s(i)) cos(s(i))];
    k=E*A(i)/l(i)*[ 1 0 -1 0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mass=zeros(m,m);

for i=1:n
    p=a(i,:);
    l=a(:,1);
    s=a(:,2);
    T=[cos(s(i)) sin(s(i)) 0 0;
        -sin(s(i)) cos(s(i)) 0 0;
        0 0 cos(s(i)) sin(s(i));
        0 0  -sin(s(i)) cos(s(i))];
    M1=RO*A(1)*l(i)/2*[      1,   0,  0 ,  0;...
                          0,   1,   0,  0;...
                          0,   0,   1,  0;...
                          0,   0,   0,  1];
                  
                  
        M=T'*M1*T;
    
          for j=1:4          
            for k=1:4
               if p(2+j)~=0 
                 if p(2+k)~=0
                    ii=p(2+j);
                    jj=p(2+k);
                    mass(ii,jj)=mass(ii,jj)+M(j,k);
                  
                end
                  end
              end
          end  
end
mass;


%%%%%%%%%%%%%%%%%%%%%%%%
damp=0.3710*mass+0.0232*stiff;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excitation data
%S=input('Force intensity of white noise   ');
S=8;
dt=input('Time step  ');
%dt=0.01;
%Duration=input('Time duration  ');
Duration=40;
SeedNum=input('Input the seed number  ');
%Filter index
Findx=1;
%Filter index =1 -  use a filter to handle the direct pass-through problem. 

%Noise level
nl=input('noise level  ');

% ***** Generate white noise
randn('state',SeedNum);
t=[dt:dt:Duration]';
Nt=length(t);

if Findx==1
   filter_order = 6;
   filter_cutoff = 20; %Hz
   [filt_num,filt_den] = butter(filter_order,filter_cutoff*2*dt);
   Nt2=Nt+2*filter_order;
else
   Nt2=Nt;
end;

ff=S*randn(Nt2,1)./dt^0.5;

if Findx==1
   ff= filter(filt_num,filt_den,ff);
   ff= ff(Nt2-Nt+1:end,:);
end;

force(:,1)=t;
force(:,2)=ff;

ll=length(ff);
force(1,1:2)=0;
force(2:ll+1,1)=t;
force(2:ll+1,2)=ff;


%%%%%%
n=m;
%%%%%%

% State transition matrix
A=zeros(2*n);
A(1:n,n+1:2*n)=eye(n);
A(n+1:2*n,1:n)=-inv(mass)*stiff;
A(n+1:2*n,n+1:2*n)=-inv(mass)*damp;
   
% Excitation localization matrix 
locat=zeros(n,1); locat(8)=1.0;
B=zeros(2*n,1);
B(n+1:2*n)=inv(mass)*locat;

% initial condition
X0=zeros(2*n,1);   % initial conditions
[Y,X]=lsim(A,B,A,B,ff,t,X0);

% Acceleration response
acc=Y(:,n+1:2*n);


% 
%  % Add measurement noise
% % % calculate the rms of all measurements
% ll=length(acc(:,1));
% for i=1:n
%    noise=randn(ll,1);
%    accn(:,i)=acc(:,i)+nl/100*std(acc(:,i))*noise;
% end;
% 
% acc=accn;

save response.mat acc  
save mass.mat mass

fftplot(acc(:,10),0.001);
 ws=eig(stiff,mass); 
   w=sqrt(ws)/2/pi   
      
