from unbinned_likelihood_CAST_g4_datasets import *
import matplotlib.pyplot as plt

g4 = 1e-40 # assume a value of g4
dataset = dataset_Ar1
c = candidate_weights(dataset_Xe,0,0)
#print(c)
#cs = [()]
candidate_pos_x = np.array(dataset.candidates["xs"])
candidate_pos_y = np.array(dataset.candidates["ys"])
cEnergies = np.array(dataset.candidates["Es"]) 

df = dataset.candidates
#df = df[(df["xs"] > -2 and df["xs"] < 2 and df["ys"] > -5 and df["ys"] < 5)] 
#df = df[(df["xs"] > -2)]
#df = df[(df["xs"] < 2)]
#df = df[(df["ys"] > -5)]
#df = df[(df["ys"] < 5)] 
print(df)
plt.scatter(df["xs"],df["ys"])
plt.show()

# I want to evaluate all s and b terms to have numbers, so that the likelihood function is in the end
# just a function of g^4, because after evaluating all the other terms they just yield a constant.
# s_tot(g4) = alpha*g4
# b_tot =  beta
# s_i(x,y,g4,E) = gamma(E) * cw(x,y) * g4, where gamma is the energy dependcence, cw is the position dependendce
# b_i(E) = delta(E) 
# L = e^(-alpha g4 - beta)*(gamma(E_1) * cw(x_1,y_1) * g4 + detla(E_1)) * ... * (gamma(E_n) * cw(x_n,y_n) * g4 + detla(E_n)) , for n candidates
print("***********************************")
print("Extended likelihood terms")
print("***********************************")
beta = totalBackground(dataset)
print("Total background (beta) = ", beta, " counts.")
alpha = totalSignal(dataset, g4) / g4 # becasue I want it to be independent of g4, but totalSignal is multiplied by g4
print("Total signal (alpha) = ", alpha)
print("***********************************")
print("Sum (if logL) or multiplication (if L) terms (one per candidate)")
print("***********************************")

cEnergies = np.array(df["Es"])
candidate_pos_x = np.array(df["xs"])
candidate_pos_y = np.array(df["ys"])
print(cEnergies)
for candidate in range(len(cEnergies)): # column names in df are xs, ys, Es
    E = cEnergies[candidate]
    x = candidate_pos_x[candidate]
    y = candidate_pos_y[candidate]
    cw = candidate_weights(dataset, x, y)
    gamma = signal(dataset, E, g4, x, y) / cw / g4
    delta = background(dataset, E)
    if cw > 0.0:
        print(f"Energy = {E}, x = {x}, y = {y}, cw = {cw}, gamma = {gamma}, delta = {delta}")

# I have computed the likelihood function by hand. I write it here to plot it. Only for dataset_Ar1 for now.
def explicitLikelihood(g4):
    minusLogL = 5.8e40 * g4 + 25.67 - (np.log(2.18e38 * g4 + 1.55) + np.log(3.59e36 * g4 + 1.9))
    return minusLogL

g4Lin = np.linspace(-3.0e-40, 1.0e-40, 1000)
l = explicitLikelihood(g4Lin)
plt.plot(g4Lin, l)
plt.savefig("fig.pdf")
print(l)

likelihood = np.exp(-l)
plt.plot(g4Lin, likelihood)

def testLimit():
    g4Lin = np.linspace(0, 1.0e-40, 1000)
    l = explicitLikelihood(g4Lin)
    likelihood = np.exp(-l)
    LCumSum = np.cumsum(likelihood)          # cumulative sum
    LMax = LCumSum.max()               # maximum of the cumulative sum
    LCdf = LCumSum / LMax              # normalize to get (empirical) CDF
    limitIdx = bisect_left(LCdf, 0.95) # limit at 95% of the CDF
    print("limit = ", g4Lin[limitIdx] ** 0.25)
testLimit()