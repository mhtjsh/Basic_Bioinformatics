def find_pattern_with_mismatches(text, pattern, max_mismatches):
    """
    Find all occurrences of a pattern in a text with at most a given number of mismatches.

    Args:
        text (str): The main sequence to search.
        pattern (str): The pattern to find.
        max_mismatches (int): Maximum allowed mismatches.

    Returns:
        list: Starting positions of the pattern in the text with at most max_mismatches.
    """
    locations = []
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            locations.append(i)
    return locations


# Example Input: Multiple DNA strings
sample_strings = [
    "XABACDXABACD",
    "CDEFCDEHABAC",
    "GHIJKLABACDa",
    "PQRSTUABACDX",
    "VWXYZABACDcA",
    "HIJKLMNOPQRST",
    "CDEFGHIJKLMN",
    "QRSTUVWXYZAB",
    "EFGHIJKLMNOP",
    "UVWXYZABCDEF"
]

pattern = "ABACD"  # Pattern to search for
max_mismatches = 2  # Allow up to 2 mismatches

# Process each string to find the pattern with mismatches
for text in sample_strings:
    locations = find_pattern_with_mismatches(text, pattern, max_mismatches)
    if locations:
        print(f"String: {text}, Pattern locations (with at most {max_mismatches} mismatches): {locations}")
    else:
        print(f"String: {text}, Pattern not found with at most {max_mismatches} mismatches.")
