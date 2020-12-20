#Coding by @Changyu 2021.12.14
#Monte Carlo Simulation for 2D Ising Model

#Package
import math
import random
import matplotlib.pyplot as plt
import numpy as np

#Several constants
N = 20                                       #lattice size
num_spin = N * N                            #Num. of spins   
lat = np.array([[0]*N]*N)                   #set lattice with all 0
num_mc = 10000                              #Monte Carlo steps
normalization = 1 / ( num_mc * num_spin )   #normalization constant
transient = 1000                            #transient steps

#Fuctions
#Initialize our lattice
def initialize(lattice):
    for i in range(N):
        for j in range(N):
            if random.random() >= 0.5 :
                lattice[i][j] =  1
            else:
                lattice[i][j] = -1

#Choose postion to flip
def choose_pos(lattice,size):
    x_pos = random.randrange(size)
    y_pos = random.randrange(size)
    return x_pos, y_pos


#Calculate energy for one dot
def dot_energy(lattice,size,x,y):
    #Boundary condition

    if x == 0:
        up = size - 1
    else:
        up = x - 1
    
    if x == size - 1:
        down = 0
    else:
        down = x + 1
    
    if y == 0:
        left = size - 1
    else:
        left = y - 1
    
    if y == size - 1:
        right = 0
    else:
        right = y + 1

    energy = ( -1 ) * lattice[x][y] * ( lattice[x][left] + lattice[x][right] + lattice[up][y] + lattice[down][y] )
    return energy

#Determine if we select flip or not
def check_flip(lattice,size,x,y,temp):
    dE = - 2 * dot_energy(lattice,size,x,y)
    if dE <= 0:
        return True
    else:
        pro = random.random()
        if  pro < math.e ** ( - dE / temp ):
            return True
        else:
            return False

#Flip the spin
def do_flip(lattice,x,y):
    lattice[x][y] = -lattice[x][y]

#Initialize first stpe (temprature)
def transient_temp(lattice,size,temp):
    for s in range(transient):
        for i in range(num_spin):
            x, y = choose_pos(lattice,size)
            if check_flip(lattice,size,x,y,temp):
                do_flip(lattice,x,y)

#Calculate E_total
def total_energy(lattice,size):
    E = 0
    for i in range(size):
        for j in range(size):
            E = E + dot_energy(lattice,size,i,j)

    return E

#Calculate magnetic
def total_magnetic(lattice,size):
    M = 0
    for i in range(size):
        for j in range(size):
            M = M + lattice[i][j]
    
    return M

#Main programm

initialize(lat)

#Divide Temp
step = 46


#some variables we need
Temp = np.linspace(5,0.5,step)
E_avg = np.linspace(0,0,step)
E_Sq_avg = np.linspace(0,0,step)
M_avg = np.linspace(0,0,step)
M_Sq_avg = np.linspace(0,0,step)
M_Qua_avg = np.linspace(0,0,step)
M_Abs_avg = np.linspace(0,0,step)

#cycle for Monte Carlo
#Temperature loop
for cyc_1 in range(step):

    transient_temp(lat,N,Temp[cyc_1])

    Ener = total_energy(lat,N)
    M = total_magnetic(lat,N)

    #initialize summation variables at each temp
    Ener_total = 0
    Ener_Square_total = 0
    Magn_total = 0
    Magn_Square_total = 0
    Magn_Abs_total = 0
    Magn_Quar_total = 0

    #Monte Carlo 
    for cyc_2 in range(num_mc):

        #Loop for each step
        for cyc_3 in range(num_spin):

            x, y = choose_pos(lat,N)
            if check_flip(lat,N,x,y,Temp[cyc_1]):

                delta_E = - 2 * dot_energy(lat,N,x,y)
                do_flip(lat,x,y)
                Ener = Ener + 2 * delta_E
                M = M + 2 * lat[x][y]
        
        Ener_total = Ener_total + Ener / 2
        Ener_Square_total = Ener_Square_total + ( Ener / 2 ) * ( Ener / 2 )
        Magn_total = Magn_total + M
        Magn_Square_total = Magn_Square_total + M * M
        Magn_Quar_total = Magn_Quar_total + M * M * M * M
        Magn_Abs_total = Magn_Abs_total + math.sqrt( M * M )

    Ener_avg = Ener_total * normalization
    Ener_Square_avg = Ener_Square_total * normalization
    Magn_avg = Magn_total * normalization
    Magn_Square_avg = Magn_Square_total * normalization
    Magn_Quar_avg = Magn_Quar_total * normalization
    Magn_Abs_avg = Magn_Abs_total * normalization

    E_avg[cyc_1] = Ener_avg
    E_Sq_avg[cyc_1] = Ener_Square_avg
    M_avg[cyc_1] = Magn_avg
    M_Sq_avg[cyc_1] = Magn_Square_avg
    M_Qua_avg[cyc_1] = Magn_Quar_avg
    M_Abs_avg[cyc_1] = Magn_Abs_avg


#Heat Capacity
# C = ( <M^2> - (N ** 2 )( <|M|>^2 ) ) / ( k_B * T )

Heat_Capacity = np.linspace(0,0,step)
for i in range(step):
    Heat_Capacity[i] = ( E_Sq_avg[i] - num_spin * E_avg[i] * E_avg[i] ) / ( Temp[i] * Temp[i] )

#Magnetic Susceptability
# K = ( <E^2> - (N ** 2 )( <E>^2 ) ) / ( k_B * ( T ** 2 ) )
Magnetic_Susceptability = np.linspace(0,0,step)
for i in range(step):
    Magnetic_Susceptability[i] = ( M_Sq_avg[i] - num_spin * M_Abs_avg[i] * M_Abs_avg[i] ) / ( Temp[i] )

print(max(Heat_Capacity))
print(max(Magnetic_Susceptability))

plt.subplot(2,2,1)
plt.plot(Temp,Heat_Capacity,'bo')

plt.subplot(2,2,3)
plt.plot(Temp,Magnetic_Susceptability,'ro')

plt.show()


"""

Reference

[1] Jacques Kotze Introduction to Monte Carlo methods for an Ising Model of a Ferromagnet (2017)

[2] Stephen J.Blundell, Katherine M.Blundell, Concepts in Thermal Physics (2010)


"""