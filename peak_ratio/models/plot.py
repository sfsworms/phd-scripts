import matplotlib.pyplot as plt
import os

# Constants for plotting
MAX_HEIGHT = 2500
MIN_HEIGHT = -100
RESOLUTION_DPI = 96
FORMAT_FILE = "png"
BASE_RANGE = 6
CHROMATOGRAM_COLORS = {"A": "tab:green",
                       "C": "tab:blue",
                       "G": "tab:gray",
                       "T": "tab:red"}


def chromatogram(sample, read, ylim):
    # Get the location (indices) of mutations in the read
    mutation_locations = read.get_locations(sample.sample_locations)
    # Define the range around the mutation to be plotted
    lower_location = mutation_locations[0]
    upper_location = mutation_locations[-1]

    # Determine the range for plotting the chromatogram around the mutation
    pic_lower_location = read.base_location[lower_location - BASE_RANGE]
    pic_upper_location = read.base_location[upper_location + BASE_RANGE]

    # Extract a subsequence from the read around the mutation for plotting going for BASE_RANGE in both direction around the 
    sequence = list(read.sequence[lower_location - BASE_RANGE:upper_location + BASE_RANGE + 1])

   # Define x-ticks for the plot based on nucleotide positions around the mutation
    xticks = []
    for base_location in range(lower_location - BASE_RANGE, upper_location + BASE_RANGE + 1):
        xticks.append(read.base_location[base_location])

    # Create a new figure and axis for plotting
    fig, ax = plt.subplots(figsize=(5, 4))
    # Plot each chromatogram (one for each nucleotide: A, C, G, T)
    for key in read.chromatograms.keys(): # For each nucleotide
        ax.plot(read.chromatograms[key], color=CHROMATOGRAM_COLORS[key])

    # Set plot title, limits, and x-ticks
    plt.title(read.name)
    plt.ylim(ylim[0], ylim[1])
    plt.xlim(pic_lower_location - 1, pic_upper_location + 1)
    plt.xticks(xticks, sequence)
    # Highlight the mutation position in red (this is based on the BASE_RANGE) 
    # TODO changer BASE_RANGE par une list
    ax.get_xticklabels()[BASE_RANGE].set_color("red")
    # Adjust the layout for optimal display
    plt.tight_layout()
    # Check if the "plot" subfolder exists, if not create it
    if not os.path.exists('plot'):
        os.makedirs('plot')
    # Save the plot to a file with the read's name
    plt.savefig(fname=os.path.join('plot', read.name + '.png'), dpi=RESOLUTION_DPI, format=FORMAT_FILE)

    # Close the plot
    plt.close()


def chromatogram_read(sample, ylim=(MIN_HEIGHT, MAX_HEIGHT)):
    chromatogram(sample=sample, read=sample.read, ylim=ylim)
