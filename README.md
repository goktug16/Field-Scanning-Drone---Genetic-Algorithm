# Field-Scanning-Drone---Genetic-Algorithm

X drones will scan a 9*9 area.
A common start and end location within the space for drones
will be given. They will all start from the same given point and return to the same place.

Drones can move in 8 directions.

X paths to be followed by drones will be found by genetic algorithm.
Paths taken by drones to scan maximum area should be as different from each other as possible.
Since sharp turns need to reduce the speed, the turning angles should be minimized.

The drones will start moving at the same time. 
Collision situations where they are assumed to move at different heights will not be considered.

Each solution proposal/individual can consist of a string of numbers denoting directions of length X*M.

In this assignment I completed using Genetic Algorithm, 
I used 3 fitness functions. These functions are in order:

F1 : I tried to minimize the function by calculating the maximum field trip, that is, 
the number of cells it did not visit, and then taking the inverse of its normalized form.

Since the values that we get from this function, where the no-go area is high, are reversed, the probability of being selected will decrease. 
A large number of squared solutions will increase the probability of being selected.

F2 : The proximity of each drone to the starting point is calculated with the formula of the euclidean distance, 
increasing the probability of selecting the drones that are closest to the start and likely to return to the start.

F3 : Since we want the movements of the drones to be with as few angles as possible, 
the function that deduces the value of drones with low angles is high and increases the probability of choosing these solutions.

Using these 3 functions, first around a main loop as many as generation. For each cycle, all solutions in the population are divided proportionally to the number of drones entered, and fitness values are created for each solution. By summing these values, a general fitness value is obtained. From these fitness values formed for each generation, the most suitable value is selected and put into the array where the best values of the generation are kept. The path used by this index is also cast to the array where the paths of the best values are kept for later use in chart drawing.

By using the roulette wheel method, it is determined which of the individuals to be created for the next generation should be crossover using the individuals that are now found, and the indices with high probability are selected. 

For each round, the best population/2 number of individuals is passed to the next generation. After the conclusion of the generation cycle, the best of the best among the series we keep the best is selected and graphs are drawn for this solution.
