import matplotlib.pyplot as plt
import csv


def process_time_points(column_list, n_time_points):
    """
    Takes the columns of a line as a list, takes the columns that contain frequency data, 
    and converts the frequency data in a a list containing the frequencies for each time point
    """
    freq_list = []
    for time_point in range(n_time_points):
        freqs = list(map(int, column_list[4+time_point].strip().split(',')))
        freq_list.append(freqs)
    print(freq_list)
    return freq_list
    
def is_below_thresh(freq_list, min_reads):
    """Check whether haplotypes have at least MIN_READS number of reads"""
    return any(any(freq < min_reads for freq in time_point) for time_point in freq_list)

def update_bin_freqs(freq_list, bin_freqs):
    num_haplotypes  = len(freq_list[0])
    bin_freqs[num_haplotypes-3] += 1

FILE_PATH = "filtered_sites.tsv"

MIN_READS = 3
MAX_NUM_HAPLOTYPES = 7

al_bins = [0] * (MAX_NUM_HAPLOTYPES-2)


with open(FILE_PATH, 'r', encoding="utf-8") as file:
    for idx, line in enumerate(file):
        columns = line.strip().split('\t')
        if idx == 0:
            num_time_points = len(columns)-4
            freqs_list = process_time_points(columns, num_time_points)
        
        freqs_list = process_time_points(columns, num_time_points)
        
        # Only consider positions where haplotypes have at least MIN_READS number of reads
        if is_below_thresh(freqs_list, MIN_READS):
            print("Was below thresold!")
            continue

        update_bin_freqs(freqs_list, al_bins)



print(al_bins)
x_values = list(range(3,len(al_bins)+3))
print(x_values)

# Write the list to the CSV file
with open("num_alleles_at_site.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write the integers as a single row
    writer.writerow(x_values)

with open("num_sites.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write the integers as a single row
    writer.writerow(al_bins)


plt.bar(x_values, al_bins, width=0.8, color='blue', edgecolor='black')

# Add text annotations (absolute numbers) on top of each bar
for x, y in zip(x_values, al_bins):
    plt.text(x, y + 0.5, str(y), ha='center', va='bottom')  # Adjust 'y + 0.5' for spacing above the bar


# Set y-axis to logarithmic scale
plt.yscale('log')
plt.ylim([0, 50000])

title_str = f"Hard filter: At least {MIN_READS} reads"
# Add labels and title
plt.xlabel('Number of alleles at locus', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
#plt.title(title_str, fontsize=16)

#plt.savefig(f"allele_count_hist_min_{MIN_READS}_reads.png", dpi=300, bbox_inches='tight', transparent=False)


# Show the plot
plt.show()
