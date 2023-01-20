# Create an empty list to store the results
results = []

# Iterate through the columns of Table A
for column in table_a.columns:
    # Merge Table A with Table B on the current column
    result = pd.merge(table_a[[column]], table_b, left_on=column, right_on='Gene', how='left')
    # Append the result to the list
    results.append(result)

# Concatenate the results into a single dataframe
final_result = pd.concat(results)
