clear all
clc
%addpath("tools")

N=3;
syms q1 q2 q3 q4 dq1 dq2 dq3 dq4 ddq1 ddq2 ddq3 ddq4 ak1 ak2 ak3 au1 au2 au3 real%non saprei interpretare queste a, magari qualcuna non deve dipendere da q
%B=sym('B', [N N])
%a=[a1 a2 a3 a4 a5 a6 a7 a8 a9];
ak=[ak1 ak2 ak3]';
au=[au1 au2 au3]';
a=[ak; au]';
q=[q1 q2 q3]';
dq=[dq1 dq2 dq3]';
ddq=[ddq1 ddq2 ddq3]';
%% Inertia matrix part
%meglio segnare le a sul quaderno
disp("M part")
B(1,1)=au1 + ak1 + 2*au2*q2 + (au3+ak2)*q2^2 + 2*ak3*q2*cos(q3);
B(1,2)=-ak3*sin(q3);
B(2,1)=B(1,2);
B(2,2)=au3 + ak2;
B(1,3)=ak1+ak3*q2*cos(q3);
B(3,1)=B(1,3);
B(2,3)=-ak3*sin(q3);
B(3,2)=B(2,3);
B(3,3)=ak1;
% B(1,4)=0;
% B(4,1)=B(1,4);
% B(2,4)=0;
% B(4,2)=B(2,4);
% B(3,4)=0;
% B(4,3)=B(3,4);
% B(4,4)=0;
B=simplify(B)
pause
%% coriolis part
disp("c part")
c=Coriolis(B, q, dq);
c=simplify(c)
c_control=[2*(au2 + (au3+ak2)*q2 + ak3*cos(q3) )*dq1*dq2 - ak3*q2*sin(q3)*(2*dq1+dq3)*dq3;
    -(au2 + (au3+ak2)*q2 + ak3*cos(q3))*dq1^2 - ak3*cos(q3)*(2*dq1+dq3)*dq3;
    ak3*q2*sin(q3)*dq1^2 + 2*ak3*cos(q3)*dq1*dq2]; %baseline to check the validity of the result
check=simplify(c - c_control)
pause
%% known and unkbown params analysis 
disp("main matrices")
T=B*ddq; 
Ym=simplify(ParamMatrix(T, a));
Yc=simplify(ParamMatrix(c, a));

S=fact_matrix2(B, q, dq);
S=simplify(S)
S_control=[(au2 + (au3+ak2)*q2 + ak3*cos(q3) )*dq2 - ak3*q2*sin(q3)*dq3, (au2 + (au3+ak2)*q2 + ak3*cos(q3) )*dq1,  -ak3*q2*sin(q3)*(dq1+dq3);
    (au2 + (au3+ak2)*q2 + ak3*cos(q3) )*dq1 - ak3*cos(q3)*dq3, 0 , ak3*cos(q3)*(dq1+dq3);
    ak3*q2*sin(q3)*dq1 + ak3*cos(q3)*dq2, ak3*cos(q3)*dq1, 0] %baseline
S_check=simplify(S-S_control) %there are multiple ways to define S such that dM-2S is skew-symmetric
pause

%
disp("known part")
Bk=subs(B, au, zeros(length(au),1));
Bk=simplify(Bk)
Ykm=simplify(ParamMatrix(Bk*ddq, ak))
ck=Coriolis(Bk, q, dq)
ck=simplify(ck);
Sk=simplify(fact_matrix(Bk, q, dq));
Ykc=simplify(ParamMatrix(ck, ak))
Yk=simplify(Ykm+Ykc)
pause


disp("unknown part")
Bu=subs(B, ak, zeros(length(ak),1));
Bu=simplify(Bu)
Yum=simplify(ParamMatrix(Bu*ddq, au))
cu=Coriolis(Bu, q, dq)
cu=simplify(cu);
Su=simplify(fact_matrix(Bu, q, dq));
Yuc=simplify(ParamMatrix(cu, au))
Yu=simplify(Yuc+Yum)
pause

%
%% check skew-symm

disp("check skew-symmetric")
trigger="skew";
dB=dtot(B, q, dq);
% dBk=dtot(Bk, q, dq);
% dBu=dtot(Bu, q, dq);
skew=simplify(dB-2*S)
study_symm(skew, trigger)
skewk=simplify(dBk-2*Sk)
study_symm(skewk, trigger)
skewu=simplify(dBu-2*Su)
study_symm(skewu, trigger)
pause
%
%% g part
%ag
%neglected in this particular case

%% full equation
%atot

%% switch to qr
disp("reference")
syms dqr1 dqr2 dqr3 dqr4 ddqr1 ddqr2 ddqr3 ddqr4 qd1 qd2 qd3 qd4 dqd1 dqd2 dqd3 dqd4 real
qd=[qd1 qd2 qd3];
dqd=[dqd1 dqd2 dqd3 dqd4];
dqr=[dqr1 dqr2 dqr3]';
ddqr=[ddqr1 ddqr2 ddqr3]';
Ykmr=subs_reference(Ykm, ddq, ddqr);
Ykcr=subs_reference(Ykc, dq, dqr);
Yumr=subs_reference(Yum, ddq, ddqr);
Yucr=subs_reference(Yuc, dq, dqr);
Yur=Yumr+Yucr
Ykr=Ykmr+Ykcr

lu=length(Yur);
TAU=eye(lu)
Kp=eye(lu); Kd=eye(lu);
au_est=simplify(TAU*Yur'*(dqr-dq))
tau=simplify( Yur*au_est + Ykr*ak + Kp*(qd-q)*Kd*(dqd-dq) );

%

%% prova


%% functions
% you can comment out this part and use the functions in "tools". To
% comment out a section, use %{ SECTION %}

%return matrix Y
function T=ParamMatrix(V, x)
for j=1:length(V)
    for i=1:length(x)
        T(j,i)=diff(V(j), x(i)); %T=col{ dA/dxi } i=1...n %maybe A(j, :)
    end
end
end

%total differentiation
function res=dtot(A, x, dx)
res=zeros(size(A,1), size(A,2)); %column vector, dimension=number of rows of the matrix/vector A
for i=1:length(x)
    res;
    var=diff(A, x(i))*dx(i);
    res=res+var; %dM/dq=dM/dq1 + dM/dq2 + ... + dM/dqn
end
end 

function matrix=Coriolis(B, q, dq)
q1=q(1); q2=q(2); q3=q(3); %q4=q(4);
b1=B(:,1);
C1=(1/2)*(jacobian(b1,q)+jacobian(b1,q)'-diff(B,q1));
b2=B(:,2);
C2=(1/2)*(jacobian(b2,q)+jacobian(b2,q)'-diff(B,q2));
b3=B(:,3);
C3=(1/2)*(jacobian(b3,q)+jacobian(b3,q)'-diff(B,q3));
% b4=B(:,4);
% C4=(1/2)*(jacobian(b4,q)+jacobian(b4,q)'-diff(B,q4))
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c3=dq'*C3*dq;
% c4=dq'*C4*dq;
c=[c1;c2;c3];
matrix=c;
end

%generate Coriolis terms in Lagrangian formulation, but with less code, and
%more portability
function vector=Coriolis2(B, q, dq)
for i=1:length(q)
   bi=B(:,i);
   Ci=(1/2)*(jacobian(bi,q)+jacobian(bi,q)'-diff(B,q(i)));
   c(i)=dq'*Ci*dq;
end
vector=simplify(c');
end

%factorization matrix S
function A=fact_matrix2(B, q, dq)
for i=1:length(q)
   bi=B(:,i);
   Ci=(1/2)*(jacobian(bi,q)+jacobian(bi,q)'-diff(B,q(i)));
   S(i,:)=dq'*Ci;
end
A=simplify(S);
end

function A=fact_matrix(B, q, dq)
q1=q(1); q2=q(2); q3=q(3);
b1=B(:,1);
C1=(1/2)*(jacobian(b1,q)+jacobian(b1,q)'-diff(B,q1));
b2=B(:,2);
C2=(1/2)*(jacobian(b2,q)+jacobian(b2,q)'-diff(B,q2));
b3=B(:,3);
C3=(1/2)*(jacobian(b3,q)+jacobian(b3,q)'-diff(B,q3));
% b4=B(:,4);
% C4=(1/2)*(jacobian(b4,q)+jacobian(b4,q)'-diff(B,q4))
s1=dq'*C1;
s2=dq'*C2;
s3=dq'*C3;
% c4=dq'*C4*dq;
S=[s1; s2; s3]
A=S;
end

%substitute numerical values to symbolic variables
function R=subs_reference(A, old, new)
R=A;
for i=1:length(new)
   R=subs(R, old(i), new(i));
end
R=simplify(R);
end

%If trigger="symm" study symmetry properties. If trigger="skew", check if
%it's skew symmetric
function value=study_symm(A, trigger)
value="yes";

if trigger=="skew"
for i=1:size(A,1)
    for j=1:size(A,2)
        if i~=j
           if simplify(A(i,j)) ~= - simplify(A(j,i))
               value="no";
               fprintf("not symmetric in: %i, %i\n", i, j)
               difference=A(i,j)+A(j,i);
               difference=simplify(difference)
               return 
           end
        end
    end
end
end %if

if trigger=="symm"
        
for i=1:size(A,1)
    for j=1:size(A,2)
        if i~=j
            if simplify(A(i,j)) ~= simplify(A(j,i))
               value="no";
               fprintf("not symmetric in: %i, %i", i, j);
               difference=A(i,j)-A(j,i);
               simplify(difference)
               return 
               end
        end
    end
 end

end
    
end