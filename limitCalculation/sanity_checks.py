from unbinned_likelihood_CAST_g4_datasets import *
#In the inner 20 mm side square I have 40*40=1600 "pixels"
c = candidate_weights(dataset_Ar1, 0, 0)*4.0 # 4 cm^2 is the area of my 2x2 cm rectangle going form 10mm to 10 mm in X and Y
print("Weight of a candidate at (0,0) if applied to the full area: ", c)
c1 = candidate_weights(dataset_Ar1, 0, 0)*0.02*0.02 #pixel of 2x2 mm^2
print("Weight of a candidate at (0,0) if applied to a small pixel of 0.02x0.02 cm: ", c1)
c2 = candidate_weights(dataset_Ar1, 0, 0)*0.01*0.01 # 1e-4 cm^2 is the area of the original 0.1x0.1mm pixels
print("Weight of a candidate at (0,0) if applied to a small pixel of 0.01x0.01 cm, i.e., the original pixel size: ", c2)
# now I do it for several squares to make it more realistic. Before I was using the interpolated value for the center and applying it to the full area --> >1
# i need to define the coordinates of the center of each square
n = 1000 #number of rectangles in each dimension
areaPerPixel = 4.0/n**2
xMin = -10
xMax = 10
xSize = 2.0 / n
c3 = candidate_weights(dataset_Xe, 0, 0)*areaPerPixel # 4 cm^2 is the area of my 2x2 cm rectangle going form 10mm to 10 mm in X and Y
print("Weight of a candidate at (0,0) if applied to a pixel out of ", n, "x", n, " pixels: ", c3)
coord = np.linspace(xMin + xSize / 2.0, xMax - xSize / 2.0, n)
sum = 0.0
for x in coord:
    for y in coord:
        sum += candidate_weights(dataset_Xe, x, y)*areaPerPixel
print(sum)
#check also with convolution

"""
If I check with n = 10000 it takes very long, so paste here the result so that I don't need to repeat it.

Weight of a candidate at (0,0) if applied to the full area:  [112.35234794]
Weight of a candidate at (0,0) if applied to a small pixel of 0.02x0.02 cm:  [0.01123523]
Weight of a candidate at (0,0) if applied to a small pixel of 0.01x0.01 cm, i.e., the original pixel size:  [0.00280881]
Weight of a candidate at (0,0) if applied to a pixel out of  10000 x 10000  pixels:  [1.12352348e-06]
[0.99975832]
"""
