"""
This program reads in the data from the CSV datafile containing information about the ALX1 gene for 4 finch species, computes the gene and beak distances from the finch species to the outgroup individual and plots the computed results.
"""

# species named in the data file finches.csv 
finch1 = 'G.conirostris_Espanola'
finch2 = 'G.conirostris_Genovesa'
finch3 = 'G.difficilis'
finch4 = 'G.magnirostris'
outgroup = 'L.noctis'

def read_input(fileName):
    ''' This function reads in the data from the CSV datafile and builds the data list. '''
    data_list = []
    finch_list = []
    
    import csv

    # Open the file to read
    with open(fileName, 'r') as csvfile:
        finch_reader = csv.reader(csvfile)

        # Merge two lists containing information for the alleles A and B and then return the list of merged lists
        for row in finch_reader:
            if row[2] == 'A':
                data_list += row
                data_list.remove('A')
                data_list[3] = float(data_list[3])

            elif row[2] == 'B':
                data_list.insert(3, row[3])
                data_list.insert(5, row[4])
                
                data_list[5] = float(data_list[5])
                finch_list.append(data_list)
                data_list = []
    
    return finch_list

def allele_dist(gene1, gene2):
    ''' This function computes and returns the Hamming distance between the two allele gene sequences given as the arguments. '''
    dist_count = 0

    # Add to count if two characters do not match
    for i in range(len(gene1)):
        if gene1[i] != gene2[i]:
            dist_count += 1
    
    return dist_count

def gene_dist(finch1, finch2):
    ''' This function computes and returns the average Hamming distance between the genes of the two individual finches given as arguments. '''
    # Call allele_dist() 4 times and average the results from the 4 calls
    sum = 0

    # Calculate distance between alleles A and B from finch1 and alleles A and B from finch2
    for i in range(2, 4):
        allele_dist_A = allele_dist(finch1[2], finch2[i])
        allele_dist_B = allele_dist(finch1[3], finch2[i])
        sum += (allele_dist_A + allele_dist_B)

    # Return average distance between the genes from finch1 and finch2
    return (sum / 4)

def beak_dist(finch1, finch2):
    ''' This function returs the average difference between the beak scores of the two individual finches given as arguments. '''
    # Compute the difference 4 times and average the results
    total_diff = 0

    for i in range(4, 6):
        allele_A_beak_diff = abs(finch1[4] - finch2[i])
        allele_B_beak_diff = abs(finch1[5] - finch2[i])
        total_diff += (allele_A_beak_diff + allele_B_beak_diff)
    
    # Return the average beak distance between the two finches
    return (total_diff / 4)


def outgroup_distance(finches, speciesName, outgroupName):
    ''' This function takes 3 arguments: finches list, name of the finch and outgroup species. It calls to gene_dist() and beak_dist() to compute the gene and beak distances between the finch and outgroup species and returns them as lists in a tuple. '''
    geneDist = []
    beakDist = []

    outgroup_list = []
    species_list = []

    # Search through the finches list to find information about the outgroup species and the finch species given as an argument
    for finch_list in finches:
        if outgroupName in finch_list:
            outgroup_list += finch_list
        
        if speciesName in finch_list:
            species_list.append(finch_list)

    # Compute and return the gene and beak distances between the outgroup species and the finch species
    for species in species_list:
        geneDist.append(gene_dist(species, outgroup_list))
        beakDist.append(beak_dist(species, outgroup_list))

    return geneDist, beakDist


def plot_data(fileName):
    ''' This function reads in the data from the CSV datafile and uses the data to plot the graph. '''
    import matplotlib.pyplot as plt

    # Read in the information from the CSV datafile
    finch_info = read_input(fileName)

    # Call to the outgroup_distance() to get the gene and beak distances between the outgroup and finch species
    finches_list = [finch1, finch2, finch3, finch4]
    marker_colors = ['b', 'y', 'r', 'm']
    i = 0

    for finch in finches_list:
        geneDist, beakDist = outgroup_distance(finch_info, finch, outgroup)
        plt.scatter(geneDist, beakDist, color=marker_colors[i], marker='D')
        i += 1

    # Graph info
    plt.title('The Gene Distance vs. The Beak Distance\nto the outgroup from the finch species')
    plt.xlabel('The Gene Distance to the outgroup')
    plt.ylabel('The Beak Distance to the outgroup')
    plt.legend(finches_list)

    # Save the figure in the .png file
    plt.savefig('finches.png')

''' Add code here to call/use/verify above functions and plot results '''

# Calls to test the functions
finch_info = read_input('finches.csv')
print(finch_info)
print()
print('Average distance between genes:', gene_dist(finch_info[1], finch_info[3]))
print('Average beak distance between the two finches:', beak_dist(finch_info[2], finch_info[5]))
print()
print('The gene and beak distances between the outgroup and finch species:', outgroup_distance(finch_info, finch2, outgroup))

# Plotting the graph
plot_data('finches.csv')

''' 
The sets of points of G.conirostris_Espanola and G.magnirostris are overlapping, or are mixed together, as well as the sets of points of G.conirostris_Genovesa and G.difficilis. However, these two couples are visually separate from each other, meaning they both have different gene and beak distance to the outgroup species. The couple of G.conirostris_Espanola and G.magnirostris have higher beak distance to the outgroup species, while G.conirostris_Genovesa and G.difficilis have higher gene distance to the outgroup species.
'''