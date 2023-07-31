import matplotlib.pyplot as plt


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
    mutation_locations = read.get_locations(sample.sample_locations)
    lower_location = mutation_locations[0]
    upper_location = mutation_locations[-1]

    pic_lower_location = read.base_location[lower_location - BASE_RANGE]
    pic_upper_location = read.base_location[upper_location + BASE_RANGE]

    sequence = list(read.sequence[lower_location - BASE_RANGE:upper_location + BASE_RANGE + 1])

    xticks = []
    for base_location in range(lower_location - BASE_RANGE, upper_location + BASE_RANGE + 1):
        xticks.append(read.base_location[base_location])

    fig, ax = plt.subplots(figsize=(5, 4))
    for key in read.chromatograms.keys():
        ax.plot(read.chromatograms[key], color=CHROMATOGRAM_COLORS[key])

    plt.title(read.name)
    plt.ylim(ylim[0], ylim[1])
    plt.xlim(pic_lower_location - 1, pic_upper_location + 1)
    plt.xticks(xticks, sequence)
#    for mutation_location in mutation_locations:  TODO changer BASE_RANGE par une list
    ax.get_xticklabels()[BASE_RANGE].set_color("red")
    plt.tight_layout()
#    plt.show()
    plt.savefig(fname=read.name, dpi=RESOLUTION_DPI, format=FORMAT_FILE)
    plt.close()


def chromatogram_read(sample, ylim=(MIN_HEIGHT, MAX_HEIGHT)):
    chromatogram(sample=sample, read=sample.read, ylim=ylim)
