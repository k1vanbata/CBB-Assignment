import pandas as pd
import matplotlib.pyplot as plt

def create_counts_dataframe(file_path):
    counts = []
    repeats = []

    # Open and read the file
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Process each line and extract counts
    for line in lines:
        parts = line.split(', ')
        if len(parts) == 2 and 'Count:' in parts[1]:
            try:
                # Extract repeat and count
                repeat = parts[0].split(': ')[1]
                count = int(parts[1].split(':')[1].strip())
                counts.append(count)
                repeats.append(repeat)  # Store the repeat sequence for example output
            except ValueError:
                # Skip lines with unexpected count format
                print(f"Warning: Could not parse count in line: {line.strip()}")
    
    # Create DataFrame from counts list
    df = pd.DataFrame({'Repeat': repeats, 'Count': counts})
    return df

# File paths for the randomized data
randomized_file_paths = [
    '1_5genes_output.txt',
    '5genes_output.txt',
    '2_5genes_output.txt',
    '3_5genes_output.txt',
    '4_5genes_output.txt'
]

# Load and concatenate the randomized sequence counts
df_randomized = pd.concat([create_counts_dataframe(file) for file in randomized_file_paths], ignore_index=True)

# Calculate statistics for the randomized sequence counts
mean_randomized = df_randomized['Count'].mean()
st_dev_randomized = df_randomized['Count'].std()
upper_limit_randomized = mean_randomized + (3 * st_dev_randomized)
lower_limit_randomized = max(0, mean_randomized - (3 * st_dev_randomized))

# Reference files for analysis
reference_files = {
    "NonCoding_analysis_output.txt": "NonCoding",
    "Coding_analysis_output.txt": "Coding"
}

# Analyze each reference file separately
for reference_file, label in reference_files.items():
    # Load the reference sequence counts
    df_reference = create_counts_dataframe(reference_file)

    # Calculate statistics for the reference sequence counts
    mean_reference = df_reference['Count'].mean()
    st_dev_reference = df_reference['Count'].std()
    upper_limit_reference = mean_reference + (3 * st_dev_reference)
    lower_limit_reference = max(0, mean_reference - (3 * st_dev_reference))

    # Print statistics for the current reference dataset
    print(f"{label} Sequence Statistics:")
    print(f"Mean Count: {mean_reference}")
    print(f"Standard Deviation of Counts: {st_dev_reference}")
    print(f"Upper Limit (Mean + 3*SD): {upper_limit_reference}")
    print(f"Lower Limit (Mean - 3*SD): {lower_limit_reference}\n")

    # Identify repeats that are significantly more or less abundant
    abundant_repeats = df_reference[df_reference['Count'] > upper_limit_randomized]
    rare_repeats = df_reference[df_reference['Count'] < lower_limit_randomized]

    # Print out a few examples of abundant and rare repeats for each reference file
    print(f"Examples of Repeats More Abundant than Expected in {label} (3 SD above mean of randomized):")
    print(abundant_repeats.head(5))  # Show top 5 examples

    print(f"\nExamples of Repeats Less Abundant than Expected in {label} (3 SD below mean of randomized):")
    print(rare_repeats.head(5))  # Show top 5 examples

    # Save all significantly more or less abundant repeats to a .txt file for the current reference file
    output_file = f"significant_repeats_analysis_{label}.txt"
    with open(output_file, 'w') as f:
        f.write(f"Repeats More Abundant than Expected in {label} (3 SD above mean of randomized):\n")
        for index, row in abundant_repeats.iterrows():
            f.write(f"Repeat: {row['Repeat']}, Count: {row['Count']}\n")
        
        f.write(f"\nRepeats Less Abundant than Expected in {label} (3 SD below mean of randomized):\n")
        for index, row in rare_repeats.iterrows():
            f.write(f"Repeat: {row['Repeat']}, Count: {row['Count']}\n")

    print(f"\nAll significantly abundant and rare repeats for {label} have been saved to {output_file}")

    # Create side-by-side subplots for histogram comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Histogram for randomized sequence counts
    axes[0].hist(df_randomized['Count'], bins=10, color='skyblue', edgecolor='black')
    axes[0].axvline(mean_randomized, color='blue', linestyle='dashed', linewidth=1.5, label=f'Randomized Mean ({mean_randomized:.2f})')
    axes[0].axvline(upper_limit_randomized, color='blue', linestyle='dotted', linewidth=1.5, label=f'Randomized +3 SD ({upper_limit_randomized:.2f})')
    axes[0].axvline(lower_limit_randomized, color='blue', linestyle='dotted', linewidth=1.5, label=f'Randomized -3 SD ({lower_limit_randomized:.2f})')
    axes[0].set_xlabel("Repeat Count")
    axes[0].set_ylabel("Frequency")
    axes[0].set_title("Distribution of Repeat Counts in Randomized Sequences")
    axes[0].legend()

    # Histogram for reference sequence counts
    axes[1].hist(df_reference['Count'], bins=10, color='lightcoral', edgecolor='black')
    axes[1].axvline(mean_reference, color='darkred', linestyle='dashed', linewidth=1.5, label=f'{label} Mean ({mean_reference:.2f})')
    axes[1].axvline(upper_limit_reference, color='darkred', linestyle='dotted', linewidth=1.5, label=f'{label} +3 SD ({upper_limit_reference:.2f})')
    axes[1].axvline(lower_limit_reference, color='darkred', linestyle='dotted', linewidth=1.5, label=f'{label} -3 SD ({lower_limit_reference:.2f})')
    axes[1].set_xlabel("Repeat Count")
    axes[1].set_ylabel("Frequency")
    axes[1].set_title(f"Distribution of Repeat Counts in {label} Sequence")
    axes[1].legend()

    plt.suptitle(f"Comparison of {label} vs. Randomized Repeat Count Distributions")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
