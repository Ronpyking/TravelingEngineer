from math import sqrt, sin, cos, atan2, radians
import random
import matplotlib.pyplot as plt

class Gene:
    """Gene created from two parent genes (or randomly when no arguments are given)"""
    def __init__(self, dominant_parent=None, recessive_parent=None):  # recessive_parent
        self.path = []    # Representation of path using indices
        self.distance_value = float()
        if (dominant_parent == None and recessive_parent == None):  # Create a gene of the first generation
            self.path = list(range(22))
            random.shuffle(self.path)
        else:
            random_initial_idx1, random_initial_idx2 = random.randint(0, 10), random.randint(0, 12)
            self.path[:12] = dominant_parent.path[random_initial_idx1:random_initial_idx1 + 12]
            self.path[12:] = recessive_parent.path[random_initial_idx2:random_initial_idx2 + 10]

            multiples, missing = [], []
            for i in range(22):
                if self.path.count(i) == 2: multiples.append(i)     # Check which ones occur twice
                elif self.path.count(i) == 0: missing.append(i)     # Check which ones are still missing

            random.shuffle(multiples), random.shuffle(missing) # Shuffle the lists to get random mutations

            for element in multiples:
                for i in reversed(range(22)):
                    if self.path[i] == element: # Find multiples replace the last multiple with a missing one
                        self.path[i] = missing.pop(0)
                        break

    def Distance(self):
        """Calculate the distance value of this gene"""
        self.distance_value = 0.0
        self.distance_value += CalculateDistance(locations["Roosendaal"], locations[location_idx_list[self.path[0]]])     # First distance
        for i in range(21): # Intermediate path distance
            self.distance_value += CalculateDistance(locations[location_idx_list[self.path[i]]], locations[location_idx_list[self.path[i+1]]])
        self.distance_value += CalculateDistance(locations[location_idx_list[self.path[-1]]], locations["Roosendaal"])   # Final distance
        return self.distance_value        

def CalculateDistance(location1, location2):
    """Calculate the spherical distance between the two locations in km (using the Haversine)"""
    r = 6371    # Radius of the sphere (Earth)
    a = (sin((radians(location1[0]) - radians(location2[0])) / 2))**2 + cos(radians(location1[0])) * cos(radians(location2[0])) \
        * (sin((radians(location1[1]) - radians(location2[1])) / 2))**2
    distance = r * 2 * atan2(sqrt(a), sqrt(1-a))
    return distance

locations = {   # Locations in degrees longitude (N) and lattitude (E)
    'Roosendaal': (773/15, 89/20), 'AMS': (188309/3600, 5717/1200), 'ANR': (92141/1800, 16057/3600), 'BRU': (36649/720, 1009/225),
    'CLR': (2523/50, 1603/360), 'DHR': (7621/144, 1721/360), 'EIN': (1029/20, 4837/900), 'ENS': (62731/1200, 8267/1200),
    'GLZ': (92821/1800, 2959/600), 'GRQ': (425/8, 79/12), 'KJK': (182947/3600, 11543/3600), 'LEY': (188857/3600, 9949/1800),
    'LGG': (182291/3600, 9797/1800), 'LID': (313/6, 5303/1200), 'LUX': (14887/300, 1396/225), 'LWR': (191623/3600, 10369/1800),
    'MST': (61099/1200, 20797/3600), 'OBL': (184553/3600, 713/150), 'OST': (46079/900, 644/225), 'RTM': (37409/720, 3197/720),
    'UDE': (92983/1800, 6829/1200), 'UTC': (62561/1200, 19019/3600), 'WOE': (11576/225, 521/120)}
location_idx_list = [key for key in locations if key != "Roosendaal"]  # Ordered list of the locations

nr_of_parents = 5000        # The number of genes that are kept and from which children will be made
nr_of_children = 10         # The number of children per pair of parents

gene_pool = []      # List containing all genes
for i in range(nr_of_parents): gene_pool.append(Gene())    # Create the first generation of genes

gene_pool.sort(key=Gene.Distance)     # Sort from best to worst
current_best = gene_pool[0].Distance()  # Find the current best
alltime_best = current_best             # Set all-time best
random.shuffle(gene_pool)               # Shuffle the gene_pool to create more variation

is_converged = False        # Used to determine when to stop
convergence_counter = 0     # Used to determine when to stop
it_nr = 0                   # Iteration (generation) number
while not is_converged:
    it_nr += 1
    previous_best = current_best
    for couple in range(int(len(gene_pool)/2)):     # Create a new generation
        for child in range(nr_of_children): gene_pool.append(Gene(gene_pool[2*couple], gene_pool[2*couple+1]))

    gene_pool.sort(key=Gene.Distance)       # Sort from best to worst
    gene_pool = gene_pool[:nr_of_parents]   # Only keep the 40 best genes
    current_best = gene_pool[0].Distance()  # Find the current best
    if current_best >= previous_best: convergence_counter += 1  # Improvement check  
    else: convergence_counter = 0    
    random.shuffle(gene_pool)               # Shuffle the gene_pool to create more variation

    if current_best < alltime_best:         # Set all-time best to current best when appropriate
        alltime_best = current_best
        print(f"Generation: {it_nr} \t Distance: {alltime_best:0.01f} km")     # Print out improvement

    if convergence_counter > 75: is_converged = True    # Quit when no improvement for 75 consecutive gens

gene_pool.sort(key=Gene.Distance)       # Sort from best to worst
print(location_idx_list[idx] for idx in gene_pool[0].path)      # Take out the best gene

# Plotting the results
x_range = [locations["Roosendaal"][1], locations[location_idx_list[gene_pool[0].path[0]]][1]]
y_range = [locations["Roosendaal"][0], locations[location_idx_list[gene_pool[0].path[0]]][0]]
plt.plot(x_range, y_range)  # Plot the initial leg
for i in range(21):     
    x_range = [locations[location_idx_list[gene_pool[0].path[i]]][1], locations[location_idx_list[gene_pool[0].path[i+1]]][1]]
    y_range = [locations[location_idx_list[gene_pool[0].path[i]]][0], locations[location_idx_list[gene_pool[0].path[i+1]]][0]]
    plt.plot(x_range, y_range)  # Plot the intermediate path
x_range = [locations[location_idx_list[gene_pool[0].path[-1]]][1], locations["Roosendaal"][1]]
y_range = [locations[location_idx_list[gene_pool[0].path[-1]]][0], locations["Roosendaal"][0]]
plt.plot(x_range, y_range)  # Plot the final leg
for location in locations:  # Plot the locations with their tags
    plt.plot(locations[location][1], locations[location][0], 'co')
    plt.annotate(location, (locations[location][1] + 0.02, locations[location][0] - 0.1))
plt.show()
