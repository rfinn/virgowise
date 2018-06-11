Xmin = [Xmin,]
Xmax = [Xmax,]
MCsteps = 10
B = 1
for i in MCsteps:
    while(np.random.random()>= np.exp(-B*E)):
