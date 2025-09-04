
% Sequential Klaman estimator and rescurive least-square estimation of

clear all

dt=input('Time step  '); % When f(t) is asuumed to be constant in [k*dt (k+1)*dt] dt must be small   
n=10; % The number of DDOF

% Measured acceleration responses
load response.mat
[ll,nn]=size(acc);

% mass matrix: 
 load mass.mat;
l=9; % The number of accelerometers
%dd=eye(n); % Location of Sensors
dd=zeros(l,n); dd(1,2)=1.0; dd(2,3)=1; dd(3,4)=1.0;
dd(4,5)=1.0; dd(5,6)=1.0; dd(6,7)=1.0;
dd(7,8)=1.0; dd(8,9)=1.0; dd(9,10)=1.0;
y=dd*acc'; % Measured acceleration responses

B=eye(10);
% B_un    - excitation influence matrix associated with the r-unknown excitation
Bl=[B(:,8)];      
B_un=inv(mass)*Bl;
G_un=dd*B_un;

% Initial values
X(1:2*n,1)=zeros(2*n,1);                         % Initial values of displacements and velocities
X(2*n+1,1)=1.3*10^3 ;         % Initial values of stiffness 
X(2*n+2,1)=1.3*10^3 ;
X(2*n+5,1)=1.3*10^3;
X(2*n+8,1)=1.3*10^3 ;
X(2*n+9,1)=1.3*10^3;
X(2*n+11,1)=1.3*10^3; 
X(2*n+3,1)=10.5*10^2;
X(2*n+4,1)=10.5*10^2;
X(2*n+6,1)=10.5*10^2;
X(2*n+7,1)=10.5*10^2;
X(2*n+10,1)=10.5*10^2;

X(3*n+2,1)=0.3710;                                  % Initial values of a
X(3*n+3,1)=0.0232;

pk=zeros(3*n+3);
pk(1:2*n,1:2*n)=eye(2*n);                  % Initial values for error covariance of matrix 
pk(2*n+1:3*n+1,2*n+1:3*n+1)=10^6*eye(n+1);      % Initial values of stiffness error covariance matrix
pk(3*n+2,3*n+2)=0.1;
pk(3*n+3,3*n+3)=0.001;
Q=10^-6;
R=1*eye(l);                                  % measurement noise 
  % Recursive Solution

% load force.txt;

% f_un=force(:,2);
f_un(1)=0;

for k=1:ll-1;
% dX/dt=A*X+B_un*force; 
% Y=E*X+G_un*force
% B_un¡¢E¡¢G_un


  A=zeros(3*n+3);                             % State transition matrix
  A(1:n,n+1:2*n)=eye(n);
  [stiff,damp,fkp,fap,fbp]=kcm(n,X(:,k));
  A(n+1:2*n,1:3*n+3)=-inv(mass)*[stiff,damp,fkp,fap,fbp];
  Fi=eye(3*n+3)+A*dt;
   
   
  hk=dd*inv(mass)*(-stiff*X(1:n,k)-damp*X(n+1:2*n,k));
  E=dd*A(n+1:2*n,:);
  
% The predicted extended state vector by numerical integration

  OPTIONS = [];
  pre=ode45(@predict,[dt*k dt*(k+1)],X(:,k),OPTIONS,n,f_un(k),B_un,mass); % Assume f_un is constant in [k*dt (k+1)*dt],dt must be small
  Xbk_1=pre.y(:,end);

  [X(:,k+1),pk_1]=klm(Xbk_1,pk,y(:,k),hk,f_un(k),Fi,E,G_un,R,Q);

   pk=pk_1;
   
 k
 [ X(2*n+1:3*n+3,k)]


 
 
% Estiamte the unkonw input by recusrive least squear estimation

  [stiff,damp,fkp,fap,fbp]=kcm(n,X(:,k+1));
  
  hk_1=dd*inv(mass)*(-stiff*X(1:n,k+1)- damp*X(n+1:2*n,k+1));
  [f_un(k+1)]=rlse(y(:,k+1),hk_1,G_un,R);
% %  f_un=detrend(f_un);

end

