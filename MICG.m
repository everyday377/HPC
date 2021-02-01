clear
clc
n=1000;
r=0.1;
%% create A  in crs format
k=1;
A.row_ptr(1)=1;
for i=1:n
    for j=1:n
        if j==i-1
            A.val(k)=r;
            A.col_ind(k)=j;
            k=k+1;
        elseif j==i
            A.val(k)=2;
            A.col_ind(k)=j;
            k=k+1;
        elseif j==i+1
            A.val(k)=1;
            A.col_ind(k)=j;
            k=k+1;
        end 
        A.row_ptr(i+1)=k;
    end
end
%% create A^t  in crs format
k=1;
At.row_ptr(1)=1;
for i=1:n
    for j=1:n
        if j==i-1
            At.val(k)=1;
            At.col_ind(k)=j;
            k=k+1;
        elseif j==i
            At.val(k)=2;
            At.col_ind(k)=j;
            k=k+1;
        elseif j==i+1
            At.val(k)=r;
            At.col_ind(k)=j;
            k=k+1;
        end 
        At.row_ptr(i+1)=k;
    end
end
%% create b 
b=mxv(A,ones(1,n));
%% The BiCG method
x0=zeros(1,n);
r0=b-mxv(A,x0);
r00=ones(1,n);
p0=r0;
p00=r00;

count=1;
residual(1)=1;

while residual(count)>=10^-6
    count=count+1;
    q0=mxv(A,p0);
    q00=mxv(At,p00);
    a0=r00*transpose(r0)/(p00*transpose(q0));
    x1=x0+a0*p0;
    r1=r0-a0*q0;
    r11=r00-a0*q00;
    residual(count)=sqrt(r1*transpose(r1)/(b*transpose(b)));
    beta0=r11*transpose(r1)/(r00*transpose(r0));
    p1=r1+beta0*p0;
    p11=r11+beta0*p00;
    r0=r1;
    r00=r11;
    p0=p1;
    p00=p11;
    x0=x1;
end


