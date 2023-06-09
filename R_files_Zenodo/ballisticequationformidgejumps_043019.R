#using takeoff angles from actual midge jumps, calculate theoretical jump distance and takeoff energy, then compare to real jumps distances and takeoffs

g = 9.81 #m/s^2 acceleration due to gravity

meantakeoffangledegrees = 63 #data from takeoff_kinematicscode; overall mean is 63 degrees = 
meantakeoffangleradians = 63 * pi/180
meanjumpdist = 77.2 * 10^-3  # convert to m; data from table jumps:  average of 77.2 ± 12.2 mm (range: 49-121 mm)
meanmass = 0.0012726 *10^-3 #convert g to kg; data from mean body sizes (sizemeans)

# d = 2 v^2 sin(θ) cos(θ) g^-1  #ballistics equation then adjusted to solve for takeoff v
v = sqrt((0.5*meanjumpdist*g/(sin(meantakeoffangleradians)*cos(meantakeoffangleradians))))
KEballistics = 0.5 * meanmass * v^2
energyperdist = KEballistics/meanjumpdist
energydensdist = energyperdist/meanmass


