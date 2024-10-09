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

# File paths for the data files
file_paths = [
    '1_5genes_output.txt',
    '5genes_output.txt',
    '2_5genes_output.txt',
    '3_5genes_output.txt',
    '4_5genes_output.txt'
]

# Load and concatenate the randomized sequence counts
df_randomized = pd.concat([create_counts_dataframe(file) for file in file_paths], ignore_index=True)

# Hardcoded reference example data
reference_counts = [2, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Example repeat counts for reference
df_reference = pd.DataFrame(reference_counts, columns=['Count'])

# Calculate statistics for the randomized sequence counts
mean_randomized = df_randomized['Count'].mean()
st_dev_randomized = df_randomized['Count'].std()
upper_limit_randomized = mean_randomized + (3 * st_dev_randomized)
lower_limit_randomized = max(0, mean_randomized - (3 * st_dev_randomized))

# Calculate statistics for the reference example counts
mean_reference = df_reference['Count'].mean()
st_dev_reference = df_reference['Count'].std()
upper_limit_reference = mean_reference + (3 * st_dev_reference)
lower_limit_reference = max(0, mean_reference - (3 * st_dev_reference))

# Print statistics for both datasets
print(f"Randomized Sequences Statistics:")
print(f"Mean Count: {mean_randomized}")
print(f"Standard Deviation of Counts: {st_dev_randomized}")
print(f"Upper Limit (Mean + 3*SD): {upper_limit_randomized}")
print(f"Lower Limit (Mean - 3*SD): {lower_limit_randomized}\n")

print(f"Reference Sequence Statistics:")
print(f"Mean Count: {mean_reference}")
print(f"Standard Deviation of Counts: {st_dev_reference}")
print(f"Upper Limit (Mean + 3*SD): {upper_limit_reference}")
print(f"Lower Limit (Mean - 3*SD): {lower_limit_reference}\n")

# Identify and print examples of repeats that are significantly more or less abundant
abundant_repeats = df_randomized[df_randomized['Count'] > upper_limit_randomized]
rare_repeats = df_randomized[df_randomized['Count'] < lower_limit_randomized]

# Print out a few examples of abundant and rare repeats
print("Examples of Repeats More Abundant than Expected (3 SD above mean):")
print(abundant_repeats.head(5))  # Show top 5 examples

print("\nExamples of Repeats Less Abundant than Expected (3 SD below mean):")
print(rare_repeats.head(5))  # Show top 5 examples

# Plot histogram for randomized sequence counts
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.hist(df_randomized['Count'], bins=10, color='skyblue', edgecolor='black')
plt.axvline(mean_randomized, color='green', linestyle='dashed', linewidth=1.5, label=f'Mean ({mean_randomized:.2f})')
plt.axvline(upper_limit_randomized, color='red', linestyle='dashed', linewidth=1.5, label=f'+3 SD ({upper_limit_randomized:.2f})')
plt.axvline(lower_limit_randomized, color='blue', linestyle='dashed', linewidth=1.5, label=f'-3 SD ({lower_limit_randomized:.2f})')
plt.xlabel("Repeat Count")
plt.ylabel("Frequency")
plt.title("Distribution of Repeat Counts in Randomized Sequences")
plt.legend()

# Plot histogram for reference example counts
plt.subplot(1, 2, 2)
plt.hist(df_reference['Count'], bins=10, color='lightcoral', edgecolor='black')
plt.axvline(mean_reference, color='green', linestyle='dashed', linewidth=1.5, label=f'Mean ({mean_reference:.2f})')
plt.axvline(upper_limit_reference, color='red', linestyle='dashed', linewidth=1.5, label=f'+3 SD ({upper_limit_reference:.2f})')
plt.axvline(lower_limit_reference, color='blue', linestyle='dashed', linewidth=1.5, label=f'-3 SD ({lower_limit_reference:.2f})')
plt.xlabel("Repeat Count")
plt.ylabel("Frequency")
plt.title("Distribution of Repeat Counts in Reference Sequence")
plt.legend()

plt.tight_layout()
plt.show()
