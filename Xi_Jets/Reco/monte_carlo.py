import numpy as np
import math
import random
import matplotlib.pyplot as plt

def f(x):
  return 1


#N = 100
fbar = 0

Ns = [100,500,1000]


monte_carlo(N):
  fbar = 0
  for i in range(0,N):
    x = random.randrange(-1,1)
    y = random.randrange(-1,1)
    if(x**2 + y**2 <= 1):
      fbar+=1

  return fbar/N

for N in Ns:
  print(monte_carlo(N))
