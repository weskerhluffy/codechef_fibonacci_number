y=1280763703
p=1999996771
y_inv=(1-x)*dos_inv
x=561530634
dos_inv=(p+1)>>1
x_inv=512305481

f=lambda n:((pow(y,n,p)-pow(y_inv,n,p))*x_inv)%p
f1=lambda n:((pow(y,n,p)-pow(y_inv,n,p)))%p
f_par=lambda n:((pow(y,n,p)-pow(y_inv,n,p))*x_inv)%p
