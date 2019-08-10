
data = [[float(j) for j in i.split(" ")[1:]] for i in open('lu980.in','r').read().split("\n")[1:] if i != ""]

def dist(a,b):
	return ((a[0]-b[0])**2 + (a[1]-b[1])**2) ** 0.5

with open('in980','w') as f:
	for i in range(len(data)):
		for j in range(len(data)):
			f.write("%d %d %f\n" % (i, j, dist(data[i], data[j])))
		

