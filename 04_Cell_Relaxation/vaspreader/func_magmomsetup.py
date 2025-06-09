#!/usr/bin/env python3

def split_string(magmom_cell_string):
    import re

    #magmom_cell = '4*0, 6*0 1*3 1*-3'

    # Initialize arrays for counts and values
    counts = []
    values = []

    # Process the input string
    for match in re.finditer(r'(\d+)\*(\-?\d+)', magmom_cell_string):
        count = int(match.group(1))
        value = int(match.group(2))
        counts.append(count)
        values.append(value)

    # Resulting array
    result = [counts, values]
    return result

def map_elements(small_header,magmom_cell):
    # Input arrays
    #quantities = [[4, 6, 1, 1], [0, 0, 3, -3]]
    #elements = [['Ta', 'Al', 'Fe'], [4, 6, 2]]

    quantities = magmom_cell
    elements = small_header
    
    # Extract data
    element_labels = elements[0]
    element_counts = elements[1]

    # Expand the element labels according to their counts
    expanded_elements = []
    for label, count in zip(element_labels, element_counts):
        expanded_elements.extend([label] * count)
        
    count = 0
    new_atomheader = []
    for f in quantities[0]:
        new_atomheader.append(expanded_elements[count])
        count += f

    # Map quantities to expanded elements
    mapped_quantities = []
    for row in quantities:
        current_mapped_row = []
        for idx, quantity in enumerate(row):
            # Match quantity with the corresponding expanded element
            current_mapped_row.append((new_atomheader[idx], quantity))
        mapped_quantities.append(current_mapped_row)


    return mapped_quantities, new_atomheader



def find_duplicates(elements):
    """
    Find indexes of duplicate values in a list.
    
    Args:
        elements (list): The list to check for duplicates.
    
    Returns:
        dict: A dictionary where keys are the duplicate values, and values are lists of indexes.
    """
    from collections import defaultdict

    index_map = defaultdict(list)

    # Populate the dictionary with indexes
    for index, value in enumerate(elements):
        index_map[value].append(index)

    # Filter out non-duplicates
    duplicates = {value: indexes for value, indexes in index_map.items() if len(indexes) > 1}
    return duplicates

