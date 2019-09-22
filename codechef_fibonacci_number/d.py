from gmpy2 import invert
p=1999996771
x=561530634
dos_inv=(p+1)>>1
y=(1+x)*dos_inv
#y=1280763703
y_inv=(1-x)*dos_inv
x_inv=int(invert(x, p))

f=lambda n:((pow(y,n,p)-pow(y_inv,n,p))*x_inv)%p
f1=lambda n:((pow(y,n,p)-pow(y_inv,n,p)))%p
f_par=lambda n:((pow(y,n,p)-pow(y_inv,n,p))*x_inv)%p
