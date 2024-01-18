from astropy.coordinates.matching import match_coordinates_sky
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from detection import detection


def match(detcoord, catalogcoord):
    """Matches every detection to the closest detection in the catalog.
    """
    index, sep2d, _ = match_coordinates_sky(detcoord, catalogcoord)
    return index, sep2d


def sep_constraint(sep, constraint):
    ''' Separation constraint on the separation list between sources.'''
    over = sep > constraint
    under = sep < constraint
    return under, over


def cut_catalog(matchcoord, catalogcoord, index, sep):
    '''Cut the catalog according to the match index and the
    separation constraint given'''
    catalog = catalogcoord[index[sep]]
    match = matchcoord[sep]
    return catalog, match


def catalog_comparison(max_sep, matchcoord, catalogcoord):
    """ Compare the catalog matches and the detections by a separation
        constraint. The function looks for the matches between the catalog
        and the detections, and by the given constraint looks for the perfect
        matches, the detections without matches, and the catalog sources
        that are unmacthed.
    """
    index, sep1 = match(matchcoord, catalogcoord)
    index2, sep2 = match(catalogcoord, matchcoord)

    under_sep, over_sep = sep_constraint(sep1, max_sep)
    _, over_switch = sep_constraint(sep2, max_sep)

    # All matches
    catalog_match = catalogcoord[index]

    # Perfect match
    under = cut_catalog(matchcoord, catalogcoord, index, under_sep)

    # Detections without matches
    over = cut_catalog(matchcoord, catalogcoord, index, over_sep)

    # Catalog sources unmatched with detections
    unmatched = cut_catalog(catalogcoord, matchcoord, index2, over_switch)

    return catalog_match, under, over, unmatched


def offset(base_coords, coords_to_offset, offset_sign='-'):
    """Moves a set of coordinates according to a given offset"""
    offset_column, offset_line = base_coords
    column, line = coords_to_offset

    if offset_sign == '+':
        new_column, new_line = column + offset_column[0], line + offset_line[0]
    elif offset_sign == '-':
        new_column, new_line = column - offset_column[0], line - offset_line[0]
    else:
        raise ValueError("Invalid offset sign. Use '+' or '-'.")

    return new_column, new_line


def square_selection(column, line, square_lim):
    """Creates a boolean array for coordinates within the specified square.
    """
    column_select = (column > square_lim[0][0]) & (column < square_lim[0][1])
    line_select = (line > square_lim[1][0]) & (line < square_lim[1][1])

    select_sq = column_select & line_select

    return select_sq


def square_cut(column, line, select_sq):
    """ Select the coordinates within the specified square"""
    column_sq, line_sq = column[select_sq], line[select_sq]
    return column_sq, line_sq


def get_square_on_image(detcoord, realcoord, first_cero, second_cero, wcs_wrd):
    """Gets the selected square on the detection and catalog coordinates """

    # Transform the ICRS coordinates to pixel
    r_column, r_line = wcs_wrd.world_to_pixel(realcoord)
    d_column, d_line = wcs_wrd.world_to_pixel(detcoord)

    # MOves the (0,0) pixel to weighted's zero
    r_object = r_column, r_line
    d_object = d_column, d_line
    r_column_w, r_line_w = offset(first_cero, r_object, '-')
    d_column_w, d_line_w = offset(first_cero, d_object, '-')

    select_sq = square_selection(r_column_w, r_line_w, second_cero)

    # Select the square without changing their pixel values
    d_column_sq, d_line_sq = square_cut(d_column_w, d_line_w, select_sq)
    r_column_sq, r_line_sq = square_cut(r_column_w, r_line_w, select_sq)

    #Changes zero from the square to (0,0)
    r_object_square = r_column_sq, r_line_sq
    d_object_square = d_column_sq, d_line_sq
    r_column_sq, r_line_sq = offset(second_cero, r_object_square, '-')
    d_column_sq, d_line_sq = offset(second_cero, d_object_square, '-')

    return r_column_sq, r_line_sq, d_column_sq, d_line_sq


def get_data(weighted_image, catalog_fits, example_filter):
    weighted = fits.getdata(weighted_image + '.fits')

    cat = fits.open(catalog_fits + '.fits')
    data = cat[1].data
    real = (data['x'], data['y'])

    header = fits.getheader(example_filter + '.fits', 0)
    wcs_world = WCS(header)
    return weighted, real, wcs_world, weighted_image


def explore_square(names, ceros, init, stop, fwhm, thresh, max_sep):
    weighted, catalog, filter_example = names
    weighted, real, wcs_world, weighted_image = get_data(weighted, catalog, filter_example)

    weighted_cero, square_cero = ceros

    objects, std_dv = detection(weighted_image + '.fits', init, stop, fwhm, thresh)
    objects = objects['x'], objects['y']
    detec = offset(weighted_cero, objects, '+')

    # Convertir a ra/dec usando UNCOVER, cualquier filtro

    detcoords = wcs_world.pixel_to_world(detec[0], detec[1])
    realcoords = wcs_world.pixel_to_world(real[0], real[1])

    line, column = square_cero
    weighted = weighted[column[0]:column[1], line[0]:line[1]]

    catalog_match, under, over, unmatched = catalog_comparison(max_sep, detcoords, realcoords)
    match_under, detcoords_under = under
    match_over,  detcoords_over = over
    unmatch_real, unmatch_det = unmatched

    # catalog_0, catalog_1, matched_0, matched_1
    a, b, c, d = get_square_on_image(detcoords_under, match_under, weighted_cero, square_cero, wcs_world)
    e, f, g, h = get_square_on_image(detcoords_over, match_over, weighted_cero, square_cero, wcs_world)
    j, k, n, o = get_square_on_image(detcoords, catalog_match, weighted_cero, square_cero, wcs_world)
    p, q, r, t = get_square_on_image(unmatch_det, unmatch_real, weighted_cero, square_cero, wcs_world)
    return a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv


def graph(a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv):
    fig, ((ax, ax1), (ax3, ax4)) = plt.subplots(2, 2)
    ax.set_title('Match under 0.15 arcsec')
    m, s = np.mean(weighted), std_dv

    ax.imshow(weighted, interpolation='nearest', cmap='gray',
              vmin=m-s, vmax=m+s, origin='lower')

    ax.scatter(a, b, color='red', marker='*', label='UNCOVER')
    ax.scatter(c, d, color='blue', marker='x', label='Toro')
    for i in range(a.size):
        ax.plot([a[i], c[i]], [b[i], d[i]], color='green')
    ax.legend()
    ax1.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')
    ax1.set_title('Toro detection unmatched')
    ax1.scatter(g, h, color='blue', marker='x')

    ax4.set_title('UNCOVER + Toro')
    ax4.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')
    ax4.scatter(j, k, color='red', marker='*', label='UNCOVER')
    ax4.scatter(n, o,  color='blue', marker='x', label='Toro')
    ax4.scatter(r, t, color='red', marker='*')

    ax3.set_title('UNCOVER detection unmatched')
    ax3.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')
    ax3.scatter(r, t, color='red', marker='*')
    return plt.show()


def explore_and_graph(names, ceros, init, stop, fwhm, thresh, max_sep):
    a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv = explore_square(names, ceros, init, stop, fwhm, thresh, max_sep)
    return graph(a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv)


names = ('weighted', 'UNCOVER_DR2_LW_SUPER_catalog', 'Images/A2744_F356W')
ceros = np.array([[(4303, 10045), (5280, 9455)], [(4477, 5103), (76, 668)]])
fwhm = 3.5
init = -5
stop = 1
thresh = 1.2
max_sep = 0.15*u.arcsec
# explore_and_graph(names, ceros, init, stop, fwhm, thresh, max_sep)
