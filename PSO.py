#Bruno Iochins Grisci
#MAY/2017
#Implementation of Particle Swarm Optimation based in the book "Essentials of Metaheuristics" by S. Luke

import random
import numpy as np
import sys
import copy

class PSO:
    
    def __init__(self, swarm_size, dimensions, n_in=6, lower_bounds=None, upper_bounds=None, movement_step=1.0, momentum_coef=0.729, personal_coef=2.05, informant_coef=2.05, global_coef=0.0, minimization=False):
    
        self.swarm_size = swarm_size # Number of particles to be used.
        self.dimensions = dimensions # Number of dimensions of each particle.
        self.number_informants = n_in # Number of particles to take into account.
        
        self.lower_bounds = None # List -> dimensions x 1 (can also be passed as a float, in this case the number is given as boundary to all dimensions)
        self.upper_bounds = None # List -> dimensions x 1 (can also be passed as a float, in this case the number is given as boundary to all dimensions)
        self.l_boundaries = False #True if there are lower boundaries
        self.u_boundaries = False #True if there are upper boundaries
    
        #This segment of the code is adapting the bound arguments
        if lower_bounds is None:      
            self.l_boundaries = False
        elif type(lower_bounds) is list:
            if len(lower_bounds) > 0:
                self.lower_bounds = copy.deepcopy(lower_bounds)
                self.l_boundaries = True
            else:
                self.l_boundaries = False 
        else:
            self.lower_bounds = [lower_bounds] * dimensions
            self.l_boundaries = True
        
        if upper_bounds is None: 
            self.u_boundaries = False
        elif type(upper_bounds) is list:
            if len(upper_bounds) > 0:
                self.upper_bounds = copy.deepcopy(upper_bounds)
                self.u_boundaries = True
            else:
                self.u_boundaries = False
        else:
            self.upper_bounds = [upper_bounds] * dimensions
            self.u_boundaries = True 
        
        self.movement_step = movement_step # How fast the particle moves. If movement_step is large, the particles make big jumps towards the better areas and can jump over them by accident. Thus a big movement_step allows the system to move quickly to best-known regions, but makes it hard to do fine-grained optimization. Just like in hill-climbing. Most commonly, movement_step is set to 1.
        self.momentum_coef = momentum_coef # How much of the original velocity is retained.
        self.personal_coef = personal_coef # How much of the personal best is mixed in. If personal_coef is large, particles tend to move more towards their own personal bests rather than towards global bests. This breaks the swarm into a lot of separate hill-climbers rather than a joint searcher.
        self.informant_coef = informant_coef # How much of the informants best is mixed in. The effect here may be a mid-ground between personal_coef and global_coef. The number of informants is also a factor (assuming theyre picked at random): more informants is more like the global best and less like the particles local best.
        self.global_coef = global_coef # How much of the global best is mixed in. If global_coef is large, particles tend to move more towards the best known region. This converts the algorithm into one large hill-climber rather than separate hill-climbers. Perhaps because this threatens to make the system highly exploitative, global_coef is often set to 0 in modern implementations.
        self.minimization = minimization # If True minimizes, if False maximizes
        
        self.swarm_location = [] # List of lists: swarm_size x dimensions -> location of each particle
        self.swarm_velocity = [] # List of lists: swarm_size x dimensions -> velocity of each particle
        self.swarm_score = []    # List: swarm_size x 1 -> score of each particle
    
        self.swarm_best_location = [] # List of lists: swarm_size x dimensions -> best location found for each particle at any moment
        self.swarm_best_score = [] # List: swarm_size x 1 -> best score found for each particle at any moment
    
        self.best_location = [] # List: 1 x dimensions -> best location found globaly at any moment 
        self.best_score = None # Number -> best score found globaly at any moment
             
        self.swarm_best_location = [None] * self.swarm_size
        if self.minimization:
            self.swarm_best_score = [sys.float_info.max] * self.swarm_size
            self.best_score = sys.float_info.max
        else:
            self.swarm_best_score = [-sys.float_info.max] * self.swarm_size
            self.best_score = -sys.float_info.max
        self.best_location = [None] * self.dimensions 
        
        self.create_swarm()
        
    def create_swarm(self):
        #Initializes the swarm with random locations and randon velocities
        for p in xrange(self.swarm_size):
            new_location = []
            new_velocity = []
            for d in xrange(self.dimensions):
                if self.l_boundaries:
                    lb = self.lower_bounds[d]
                else: 
                    lb = -1000000000000000.0
                if self.u_boundaries:
                    ub = self.upper_bounds[d]
                else:
                    ub = 1000000000000000.0
                #print(lb, ub)
                #print(random.uniform(lb, ub))
         
                if self.l_boundaries or self.u_boundaries:
                    l = random.uniform(lb, ub)
                    d = random.uniform(lb, ub)
                else:
                    l = random.gauss(0.0, 1000000000000000.0)
                    d = random.gauss(0.0, 1000000000000000.0)
                v = (d - l) / 2.0
                new_location.append(l)
                new_velocity.append(v)
            self.swarm_location.append(new_location)
            self.swarm_velocity.append(new_velocity)
        
    def update_velocity(self):
        for p in xrange(self.swarm_size):
            current_velocity = np.array(copy.deepcopy(self.swarm_velocity[p]))
            current_location = np.array(copy.deepcopy(self.swarm_location[p]))
            best_personal_location = np.array(copy.deepcopy(self.swarm_best_location[p]))
            best_informant_location = np.array(self.get_informants_best_location(p))
            best_global_location = np.array(copy.deepcopy(self.best_location))
            
            pc = self.personal_coef * np.random.random_sample((self.dimensions,))
            ic = self.informant_coef * np.random.random_sample((self.dimensions,))
            gc = self.global_coef * np.random.random_sample((self.dimensions,))
            
            new_vel = self.momentum_coef * current_velocity 
            new_vel = new_vel + (pc * (best_personal_location - current_location))
            new_vel = new_vel + (ic * (best_informant_location - current_location))
            new_vel = new_vel + (gc * (best_global_location - current_location))
            
            #This segment of the code changes the velocity if the result does not respect the boundaries
            new_velocity = copy.deepcopy(list(new_vel))        
            for d in xrange(self.dimensions):
                if self.l_boundaries:
                    if self.swarm_location[p][d] + (self.movement_step * new_velocity[d]) < self.lower_bounds[d]:
                        if (self.movement_step * new_velocity[d]) != 0.0:
                            reductor = (self.lower_bounds[d] - self.swarm_location[p][d]) / (self.movement_step * new_velocity[d]) * 0.99
                            new_velocity = [di * reductor for di in new_velocity] #reductor multiplied to all dimensions in order to preserve velocity direction                         
                if self.u_boundaries:                
                    if self.swarm_location[p][d] + (self.movement_step * new_velocity[d]) > self.upper_bounds[d]: 
                        if (self.movement_step * new_velocity[d]) != 0.0:
                            reductor = (self.upper_bounds[d] - self.swarm_location[p][d]) / (self.movement_step * new_velocity[d]) * 0.99
                            new_velocity = [di * reductor for di in new_velocity] #reductor multiplied to all dimensions in order to preserve velocity direction
            
            self.swarm_velocity[p] = new_velocity            
                                       
    def get_informants_best_location(self, particle_index):
        #Uses a ring layout to define informants
        best_index = particle_index
        best_location = copy.deepcopy(self.swarm_best_location[best_index])
        best_score = self.swarm_best_score[best_index]
        for i in xrange(1, self.number_informants):
            current_index = particle_index + i
            if current_index >= self.swarm_size:
                current_index = current_index - self.swarm_size 
            if self.is_better(self.swarm_best_score[current_index], best_score):
                best_index = current_index
                best_location = copy.deepcopy(self.swarm_best_location[best_index])
                best_score = self.swarm_best_score[best_index]
        return best_location  
            
    def update_location(self):
        current_locations = np.matrix(copy.deepcopy(self.swarm_location))
        velocity = np.matrix(copy.deepcopy(self.swarm_velocity))
        self.swarm_location = np.matrix.tolist(current_locations + (self.movement_step * velocity))
                 
    def update_best(self):
        for p in xrange(self.swarm_size):
            if self.is_better(self.swarm_score[p], self.swarm_best_score[p]):
                self.swarm_best_score[p] = self.swarm_score[p]
                self.swarm_best_location[p] = copy.deepcopy(self.swarm_location[p])
                if self.is_better(self.swarm_score[p], self.best_score):
                    self.best_score = self.swarm_score[p]
                    self.best_location = copy.deepcopy(self.swarm_location[p])
        
    def is_better(self, x, y):
        if self.minimization:
            return x <= y
        else:
            return x >= y
            
    def get_locations(self):
        return self.swarm_location
        
    def get_velocities(self):
        return self.swarm_velocity
        
    def get_best_locations(self):
        return self.swarm_best_location
        
    def get_best_location(self):
        return self.best_location
        
    def get_best_score(self):
        return self.best_score
        
    def run_step(self, scores, movement_step=None):
        self.swarm_score = copy.deepcopy(scores)
        if movement_step is not None:
            self.movement_step = movement_step
        self.update_best()
        self.update_velocity()
        self.update_location()  
