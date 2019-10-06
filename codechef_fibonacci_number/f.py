from sympy import isprime
N=1000
primos_set=set()
for k in range(1,N):
	n1=k*5+1
	n2=k*5-1
	if(isprime(n1)):
		primos_set.add(n1)
	if(isprime(n2)):
		primos_set.add(n2)


primos=sorted(primos_set)
#print(primos)
# https://stackoverflow.com/questions/13326673/is-there-a-python-library-to-list-primes
inicio=6
for p in primos:
	limite=min(100,p)
	tam=(limite>>1)+((p-1)-max(0,p-1-(limite>>1)))
	arch="in{:02d}.txt".format(inicio)
	with open(arch, 'w') as file:  # Use file to refer to the file object
		file.write("{}\n".format(tam))
		for c in range(limite>>1):
	        	file.write("{} {}\n".format(c,p))
		for c in range(p-1,max(0,p-1-(limite>>1)),-1):
		        file.write("{} {}\n".format(c,p))
	inicio+=1
