p=1999996771
p=1999999609
p=1009
limite=100
print(limite)
for c in range(limite>>1):
	print("{} {}".format(c,p))
for c in range(p-1,p-1-(limite>>1),-1):
	print("{} {}".format(c,p))
