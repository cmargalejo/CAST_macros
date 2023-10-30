from unbinned_likelihood_CAST_g4_datasets import *
#In the inner 20 mm side square I have 40*40=1600 "pixels"
c = candidate_weights(dataset_Ar1, 0, 0)*4.0 # 4 cm^2 is the area of my 2x2 cm rectangle
print("Weight of a candidate at (0,0): ", c)
# now I do it for several squares to make it more realistic. Before I was using the interpolated value for the center and applying it to the full area --> >1
# i need to define the coordinates of the center of each square
n = 1000 #number of rectangles in each dimension
area = 4.0/n**2
xMin = -10
xMax = 10
xSize = 2.0 / n
coord = np.linspace(xMin + xSize / 2.0, xMax - xSize / 2.0, n)
sum = 0.0
for x in coord:
    for y in coord:
        sum += candidate_weights(dataset_Xe, x, y)*area
print(sum)
