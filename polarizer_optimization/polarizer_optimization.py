# This script performs an analysis on the optimal angle between 
# the two polarizers of the setup, checking which angle produces
# the highest amount of contrast between fringes

# NOTE: THE SYSTEM MATRIX MULTIPLIES TO THE LEFT!!!  x --> Sx NOT x --> xS !! 

import numpy as np
import matplotlib.pyplot as plt 

wavelength       = 700 * 10**(-9)          # The laser wavelength
biref_thickness  = 0.1                     # The birefringent block length
blockArea        = 0.0001                  # The area of the face to which stress is applied
idenitity_matrix = np.array([[1,0],[0,1]]) # 2x2 idenitity matrix
C_upper          = 1.0*10**(-10)           # Upper estimate for C
C_lower          = 5.3*10**(-12)           # Lower estimate for C
mass_max         = 100                     # Maximum mass to put on the block
g                = 9.81                    # Gravitational acceleration

def add_polarizer(smatrix, angle):
    '''Adds a polarizer to the system matrix
    
    Arguments:
    smatrix -- np array showing the current system matrix
    angle -- variable showing the polarizer angle    
    '''
    # Construct the polarizer matrix
    polarizerMatrix = np.array([[np.cos(angle)*np.cos(angle), np.cos(angle)*np.sin(angle)],
                                [np.cos(angle)*np.sin(angle), np.sin(angle)*np.sin(angle)]])
    
    # Return the system matrix with the polarizer multiplied to the left
    return polarizerMatrix @ smatrix 

def add_birefingence(smatrix, thickness, index_x = None, index_y = None, index_d = None):
    '''Adds a birefingent block to the system matrix
    
    Arguments:
    smatrix   -- np array showing the current system matrix
    thickness -- variable showing the thickness of the block 

    Keyword arguments:
    index_x   -- the index of refraction in the x direction
    index_y   -- the index of refraction in the y direction
    index_d   -- the difference between indices of refraction (y - x)
    '''

    global wavelength

    # Compute the difference y - x, check for misuse of the function
    if  ((index_x is None ) or (index_y is None)) and (index_d is not None):
        pass
    elif((index_x is not None ) or (index_y is not None)) and (index_d is None):
        index_d = index_y - index_x
    else:
        raise Exception("Please input either both x,y indices or the difference")

    # Get the offset in the y term caused by the block
    offset =  (2*np.pi / wavelength) * index_d * thickness

    # Construct the birefringent matrix
    birefingentMatrix=np.array([[1, 0],
                               [0, np.e**(complex(imag = 1) * offset)]])
    
    return birefingentMatrix @ smatrix

def get_intensity (smatrix, input_laser):
    '''Gets the (normalized) intensity of the laser
    
    Arguments:
    smatrix     -- np array showing the current system matrix
    input_laser -- np array showing the incoming laser matrix
    '''
    # Compute the outgoing laser matrix
    going_laser = smatrix @ input_laser
    
    # Compute the outgoing/incoming laser intensities
    h_conj_buff = going_laser.transpose()
    h_conj_buff = h_conj_buff.conjugate()
    going_intensity = h_conj_buff @ going_laser

    h_conj_buff = input_laser.transpose()
    h_conj_buff = h_conj_buff.conjugate()
    input_intensity = h_conj_buff @ input_laser

    # Divide the two to normalize the going laser, return the value
    norm_intensity = going_intensity[0][0]/input_intensity[0][0]
    if norm_intensity.imag == 0:
        norm_intensity = norm_intensity.real
    else:
        print("Something went wrong: complex normalized intensity")
    
    return norm_intensity

def auto_compute_system(pangle, rindex):
    global initial_laser
    global idenitity_matrix

    sys = add_polarizer(idenitity_matrix, np.pi/4)                         
    sys = add_birefingence(smatrix=sys, thickness=biref_thickness, index_d=rindex)
    sys = add_polarizer(sys, pangle)

    return get_intensity(sys, initial_laser)

# HERE BEGINS THE MAIN CODE

# Write some arbitrary initial laser polarization
initial_laser = np.array([[1],[1]])

# Make the system
sys = idenitity_matrix                                       # Nothing
sys = add_polarizer(sys, np.pi/4)                            # P1
sys = add_birefingence(sys, biref_thickness, index_d=0)      # Block
sys = add_polarizer(sys, 3 * np.pi/4)                        # P2

# Print the outgoing intensity
print(get_intensity(sys, initial_laser))

# Get figure and axis objects
fig, ax = plt.subplots(ncols = 2, nrows = 1, sharey=True, constrained_layout=True)
fig.set_size_inches((13,6))


# Loop over some systems
polarizer_angles = np.linspace(3*np.pi/12, 9*np.pi/12, 7)
index_refr_lower = np.linspace(0, C_lower*(mass_max*g)/(blockArea),1000)
index_refr_upper = np.linspace(0, C_upper*(mass_max*g)/(blockArea),1000)

# LOWBALL ESTIMATE
# Meshgrid the polarizer angles and the refractive indices
pangles_lower, rindexs_lower = np.meshgrid(polarizer_angles, index_refr_lower, indexing="ij")
transmission_lower           = np.empty_like(pangles_lower)

# Get the transmission for each angle and index
for i in range(len(polarizer_angles)):
    for j in range(len(index_refr_lower)):
            transmission_lower[i,j]=auto_compute_system(pangles_lower[i,j], rindexs_lower[i,j])

# For each angle, plot the transmission against the mass placed on top of the PMMA block
for angle in polarizer_angles:
    ax[0].plot(index_refr_lower*(blockArea)/(C_lower*g), 
               transmission_lower[list(polarizer_angles).index(angle),:], 
               label = f"P2: {round(angle/np.pi - 0.25,2)}, maxDiff: {round(max(transmission_lower[list(polarizer_angles).index(angle),:])-min(transmission_lower[list(polarizer_angles).index(angle),:]),2)}")

# Refine the figure
ax[0].set_title ("Lowball estimate, T vs mass, different P2-P1 angles")
ax[0].set_ylabel("Transmission T")
ax[0].set_xlabel("Mass on the block, in kg")
ax[0].legend()

# HIGHBALL ESTIMATE
# Meshgrid the polarizer angles and the refractive indices
pangles_upper, rindexs_upper = np.meshgrid(polarizer_angles, index_refr_upper, indexing="ij")
transmission_upper           = np.empty_like(pangles_upper)

# Get the transmission for each angle and index
for i in range(len(polarizer_angles)):
    for j in range(len(index_refr_upper)):
            transmission_upper[i,j]=auto_compute_system(pangles_upper[i,j], rindexs_upper[i,j])

# For each angle, plot the transmission against the mass placed on top of the PMMA block
for angle in polarizer_angles:
    ax[1].plot(index_refr_upper*(blockArea)/(C_upper*g), 
               transmission_upper[list(polarizer_angles).index(angle),:], 
               label = f"P2: {round(angle/np.pi - 0.25,2)}, maxDiff: {round(max(transmission_upper[list(polarizer_angles).index(angle),:])-min(transmission_upper[list(polarizer_angles).index(angle),:]),2)}")

# Refine the figure
ax[1].set_title(f"Highball estimate, T vs mass, different P2-P1 angles")
ax[1].set_xlabel("Mass on the block, in kg")
ax[1].legend()

plt.show()

