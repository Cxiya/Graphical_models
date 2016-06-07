
#get data from txt and change character to number
import numpy as np
from os import walk
import os
os.chdir('/Users/Claire/Desktop/graphical_model/HW1c/proteins')
mypath = '/Users/Claire/Desktop/graphical_model/HW1c/proteins'  # use path to data files
_, _, filenames = next(walk(mypath), (None, None, []))

mSeq = len(filenames)        # read in each sequence
o,x = [],[]
for i in range(1,mSeq):
    f = open( filenames[i] , 'r')
    o.append( f.readline()[:-1] )  # strip trailing '\n'
    x.append( f.readline()[:-1] )
    f.close()


xvals, ovals = set(),set()  # extract the symbols used in x and o
for i in range(0,mSeq-1):
    xvals |= set(x[i])
    ovals |= set(o[i])
xvals = list( np.sort( list(xvals) ) )
ovals = list( np.sort( list(ovals) ) )
dx,do = len(xvals),len(ovals)

for i in range(0,mSeq-1):       # and convert to numeric indices
    x[i] = np.array([xvals.index(s) for s in x[i]])
    o[i] = np.array([ovals.index(s) for s in o[i]])
    
    
#From your data, estimate p(x0), the initial state distribution for a given sequence.
p0 = np.zeros(dx)

for i in range(mSeq-1):
  p0[x[i][0]] +=1
  
p0 = p0/p0.sum()

print p0

#estimate the transition probability matrix
T = np.zeros((dx,dx))

for i in range(mSeq-1):
  length = len(x[i])
  for j in range(1,length):
    T[x[i][j-1]][x[i][j]] +=1
    
for i in range(dx):
  T[i] = T[i]/T[i].sum()
  if i < 5:
    print T[i][0:5]


#Estimate the emission probability matrix, Oik = Pr[ot = k|xt = i]. (Again, print Oik for i, k in the  rst  ve values.)
O = np.zeros((dx,do))

for i in range(mSeq-1):
  length = len(x[i])
  for j in range(length):
    O[x[i][j]][o[i][j]] +=1
    
for i in range(dx):
  O[i] = O[i]/O[i].sum()
  if i < 5:
    print(O[i][0:5])
    

#Complete the code to compute the marginal probability of the state xt given the observation sequence O = [o1,...,oL], for each t
def markovMarginals(x,o,p0,Tr,Ob):
  dx,do = Ob.shape
  L = len(o)
  f = np.zeros((L,dx))
  r = np.zeros((L,dx))
  p = np.zeros((L,dx))
  f[0,:] = p0 * Ob[:,o[0]]
  log_pO = np.log(f[0,:].sum())
  f[0,:] /= f[0,:].sum()
  
  for t in range(1,L):
    f[t,:] = f[t-1,:].dot(Tr) * Ob[:,o[t]]
    log_pO += np.log(f[t,:].sum())
    f[t,:] /=f[t,:].sum()
    
  r[L-1,x[L-1]] = 1.0
  p[L-1,:] = r[L-1,:]

  for t in range(L-2,-1,-1):
    tmp = Ob[:,o[t+1]].transpose() * r[t+1,:]
    r[t,:] = tmp.dot(T.transpose())
    r[t,:] /=r[t,:].sum()
    p[t,:] = f[t,:] * r[t,:]
    p[t,:] /=p[t,:].sum()
  
  return log_pO, p

#verify some outcome
log_pO0,p_0 = markovMarginals(x[0],o[0], p0,T,O)
print p_0[6]

log_pO2,p_2 = markovMarginals(x[2],o[2], p0,T,O)
print p_2[9]

log_pO4,p_4 = markovMarginals(x[4],o[4], p0,T,O)
print log_pO4


