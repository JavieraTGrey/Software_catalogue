# Esteb codigo es para generalizar lo que esta ocurriendo
# en los otros codigos, para así despues utilizarlo sin problema!

import numpy as np


def parse_ranges(range_str):
    ranges = []
    try:
        # Use literal_eval to safely evaluate the string as a Python literal
        from ast import literal_eval
        ranges = literal_eval(range_str)
        if not isinstance(ranges, tuple) and all(isinstance(item, tuple) and len(item) == 2 for item in ranges):
            raise ValueError("Invalid range format")
    except (ValueError, SyntaxError):
        raise ValueError("Invalid range format")

    return ranges


def get_user_input():
    filters_to_combine = input("Enter filters to combine (comma-separated): ").split(',')

    # Parse cutting_edges ranges
    cutting_edges_input = input("Enter cutting edges ranges (comma-separated): ")
    cutting_edges = parse_ranges(cutting_edges_input)

    # Create the dictionary
    user_params = {
        'filters_to_combine': filters_to_combine,
        'cutting_edges': cutting_edges,
        # Add more parameters as needed
    }

    return user_params


if __name__ == "__main__":
    try:
        user_params = get_user_input()
        print("User Parameters:")
        for key, value in user_params.items():
            print(f"{key}: {value}")
    except ValueError as e:
        print(f"Error: {e}")

# Now lets start the cutting part
print('Cutting files...')
cut = user_params['cutting_edges']
column = np.array((cut[0]))
line = np.array((cut[1]))


for i in user_params['filters_to_combine']:
    print(i)
