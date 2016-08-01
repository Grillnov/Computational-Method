%solving CPT6
A=[6,2,1;2,3,1;1,1,1];
A=A-2*[1,0,0;0,1,0;0,0,1];
Ain=inv(A);
z=[1;0;0];
k=0;
while(k<100)
    y=Ain*z;
    val=norm(y,'inf');
    z=y/val;
    k=k+1;
end
disp(val^-1+2);
disp(z);