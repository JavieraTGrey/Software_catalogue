# Esteb codigo es para generalizar lo que esta ocurriendo
# en los otros codigos, para así despues utilizarlo sin problema!

import numpy as np
from weight_images import cutting_images


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
    filters_to_combine = input("filters to combine (name of fits file): ")
    filters_to_combine = filters_to_combine.split(',')

    # Parse cutting_edges ranges
    cutting_edges_input = input("cutting edges ranges (comma-separated):")
    cutting_edges = parse_ranges(cutting_edges_input)

    # Create the dictionary
    user_params = {
        'filters_to_combine': filters_to_combine,
        'cutting_edges': cutting_edges,
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

filters = user_params['filters_to_combine']
want_im = user_params['cutting_edges']

cutting_images(line, column, filters)


# Images/A2744_F277W,Images/A2744_F356W,Images/A2744_F444W
# (4303, 10045),(5280, 9455)

# names = ('weighted', 'UNCOVER_DR2_LW_SUPER_catalog', 'Images/A2744_F356W')
# ceros = np.array([[(4303, 10045), (5280, 9455)], [(4477, 5103), (76, 668)]])
# fwhm = 3.5
# max_sep = 0.15*u.arcsec
# explore_and_graph(names, ceros, fwhm, max_sep)
