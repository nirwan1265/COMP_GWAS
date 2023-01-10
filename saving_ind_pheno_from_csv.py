import csv

# Define the file to process
file_path = "Sorghum_allphospho_africa.csv"

# Open the file
with open(file_path) as file:
    # Use the CSV reader to read the file
    reader = csv.reader(file)
    # Get the header row
    header = next(reader)
    # Loop through the columns
    for column_name in header:
        # Create an empty list to store the values for this column
        values = []
        # Iterate over the rows
        for row in reader:
            # Append the value of this column to the list
            values.append(row[header.index(column_name)])
        # reset the reader to the top of the file
        file.seek(0)
        # discard the header line
        next(reader)
        # Save the values to a file with the column name as the file name
        with open(column_name+".txt", "w") as output:
            output.write("\n".join(values))

