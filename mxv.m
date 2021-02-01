function b=mxv(A,x)
[~,n]=size(A.row_ptr);
n=n-1;
for i=1:n
    b(i)=0;
    for j=A.row_ptr(i):A.row_ptr(i+1)-1
        b(i)=b(i)+A.val(j)*x(A.col_ind(j));
    end
end
        