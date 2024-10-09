import pandas as pd

def create_counts_dataframe(data):
    counts = []

    # Split the data into lines
    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        parts = line.split(', ')
        if len(lines) == 2:
            count = int(parts[1].split(':')[1])
            counts.append(count)
    
    
    for line in lines:
        # Extract the count
        parts = line.split(', ')
        if len(parts) == 2:
            count = int(parts[1].split(': ')[1])
            counts.append(count)

    # Create a DataFrame
    df = pd.DataFrame(counts, columns=['Count'])
    return df

# Example input data
input_data1 = '1_5genes_output.txt'
input_data2 = '5genes_output.txt'
input_data3 = '2_5genes_output.txt'
input_data4 = '3_5genes_output.txt'
input_data5 = '4_5genes_output.txt'

# Create the DataFrame
df1 = create_counts_dataframe(input_data1)
df2 = create_counts_dataframe(input_data2)
df3 = create_counts_dataframe(input_data3)
df4 = create_counts_dataframe(input_data4)
df5 = create_counts_dataframe(input_data5)

df_combined = pd.concat([df1, df2, df3, df4, df5], ignore_index=True)

# Output the DataFrame
print(df_combined)

# Output mean and st dev as variables

mean = df_combined['Count'].mean()
st_dev = df_combined['Count'].std()
Max = mean + (3 * st_dev)
Min = mean - (3 * st_dev)

print(f'Mean: {mean}')
print(f'Standard Deviation: {st_dev}')
print(f'Max: {Max}')
print(f'Min: {Min}')