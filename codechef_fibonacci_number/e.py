p=1999996771
p=1999999609
p=1009
p=11
limite=min(100,p)
print(limite)
for c in range(limite>>1):
	print("{} {}".format(c,p))
for c in range(p-1,max(0,p-1-(limite>>1)),-1):
	print("{} {}".format(c,p))
