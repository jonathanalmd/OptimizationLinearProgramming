'''
@file 	sat-mod-v2.py
@author	Jonathan Mendes de Almeida
@email	jonathanalmd@gmail.com / jonathan@aluno.unb.br
@page	jonyddev.github.io
@date	04/29/2019 
@info	MSc Research at Computer Networks Lab (COMNET) -- University of Brasília (UnB)
@brief	Full Normal Model using scipy.stats - An initial model for antenna saturation
		100
'''

# Imports
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math

# Urban scenario
# Set normal parameters
mu_urb = 10
variance_urb = 2
sigma_urb = math.sqrt(variance_urb)
# Create X axis values
x_urb = np.linspace(0, 23.99, 100)
# Create normal probability density function
normal_urb = stats.norm.pdf(x_urb, mu_urb, sigma_urb)

# Residential scenario
# Set normal parameters
mu_res = 18
variance_res = 2
sigma_res = math.sqrt(variance_res)
# Create X axis values
x_res = np.linspace(0, 23.99, 100)
# Create normal probability density function
normal_res = stats.norm.pdf(x_res, mu_res, sigma_res)


# Plot
# Set X and Y limitsr
plt.xlim(0,24)
plt.ylim(0,1)
# Plot distributions
plt.plot(x_urb, normal_urb, 'r')
plt.plot(x_res, normal_res, 'b')
# Show the plot
plt.show()