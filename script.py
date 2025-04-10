import matplotlib.pyplot as plt
import argparse

HYDROPATHY_INDEX = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5,
    'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
    'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
    'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
DEFAULT_WINDOW_SIZE = 6

def parse_arguments():
    parser = argparse.ArgumentParser(description='Hydropathy plotter')

    parser.add_argument('-s', metavar='widow_size', type=int, default=DEFAULT_WINDOW_SIZE, help='Define a custom window size')
    parser.add_argument('input_file', help='FASTA file with protein sequence')

    args = parser.parse_args()
    return args.s, args.input_file


def readInput(input):
    sequenceName = ''
    with open(input, 'r') as file:
        sequence = []
        firstLine = True
        for line in file:
            if firstLine and line.startswith('>'):
                firstLine = False
                sequenceName = line[1:]
                continue
            sequence.append(line.strip())
    if len(sequenceName) == 0:
        sequenceName = 'unknown sequence\n'
    return ''.join(sequence), sequenceName


def calculate_hydropathy(sequence, window_size):
    hydropathy = []
    for i in range(len(sequence)):
        window = sequence[max(0, i - window_size // 2):min(len(sequence), i + window_size // 2 + 1)]
        values = [HYDROPATHY_INDEX.get(aa, 0) for aa in window]
        hydropathy.append(sum(values) / len(values))
    return hydropathy


def plot_hydropathy(sequence, hydropathy, threshold):
    positions = range(len(sequence))

    plt.figure(figsize=(12, 6))
    plt.plot(positions, hydropathy, label="Hydropathy", color='blue')
    plt.axhline(threshold, color='red', linestyle='--', label=f'Threshold = {threshold}')

    for i, value in enumerate(hydropathy):
        if value >= threshold:
            plt.axvspan(i - 0.5, i + 0.5, color='orange', alpha=0.3)

    plt.title("Kyte-Doolittle Hydropathy Plot")
    plt.xlabel("Amino Acid Position")
    plt.ylabel("Hydropathy")
    plt.legend()
    plt.grid(True)
    plt.show()


def run(window_size, input_file):
    seq, name = readInput(input_file)
    hydropathy_values = calculate_hydropathy(seq, window_size)
    plot_hydropathy(seq, hydropathy_values, threshold=1.6)


window_size, input_file = parse_arguments()
run(window_size, input_file)
