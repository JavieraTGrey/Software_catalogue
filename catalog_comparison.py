from astropy.coordinates.matching import match_coordinates_sky
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from detection import detection
from astropy.io import fits


def match(list1, list2):
    # Compara deteccion entre ambos
    index, sep2d, dist3d = match_coordinates_sky(list1, list2)
    # reordeno realcoords
    return index, sep2d


def sep_constraint(sep, constraint):
    over = sep > constraint
    under = sep < constraint
    return under, over


def cut_catalog(matchcoord, catalogcoord, index, sep):
    catalog = catalogcoord[index[sep]]
    match = matchcoord[sep]
    return catalog, match


def catalog_comparison(max_sep, matchcoord, catalogcoord):
    index, sep1 = match(matchcoord, catalogcoord)
    index2, sep2 = match(catalogcoord, matchcoord)

    under_sep, over_sep = sep_constraint(sep1, max_sep)
    under_switch, over_switch = sep_constraint(sep2, max_sep)

    # whole catalog
    catalog_match = catalogcoord[index]

    # Busca match perfecto!
    under = cut_catalog(matchcoord, catalogcoord, index, under_sep)

    # Busca las que estan en mi catalogo pero no en UNCOVER
    over = cut_catalog(matchcoord, catalogcoord, index, over_sep)

    # Busca las que estan en UNCOVER pero no en Catalogo
    unmatched = cut_catalog(catalogcoord, matchcoord, index2, over_switch)

    return catalog_match, under, over, unmatched


def offset(ceros, column, line, sign):
    offset_column, offset_line = ceros
    if sign == 'm':
        new_column, new_line = column - offset_column[0], line - offset_line[0]
    else:
        new_column, new_line = column + offset_column[0], line + offset_line[0]

    return new_column, new_line


def square_selection(column, line, cero):
    # Ahora buscamos el cuadrado seleccionado

    column_select = (column > cero[0][0]) & (column < cero[0][1])
    line_select = (line > cero[1][0]) & (line < cero[1][1])

    select_sq = column_select & line_select

    return select_sq


def square_cut(column, line, select_sq):
    column_sq, line_sq = column[select_sq], line[select_sq]
    return column_sq, line_sq


def get_square_on_image(detcoord, realcoord, first_cero, second_cero, wcs_wrd):
    # Buscamos ahora el recorte
    # Pasamos de ra/dec a pixeles con wcs de UNCOVER
    r_column, r_line = wcs_wrd.world_to_pixel(realcoord)
    d_column, d_line = wcs_wrd.world_to_pixel(detcoord)

    # Corregimos al cero de la imagen weighted
    r_column_w, r_line_w = offset(first_cero, r_column, r_line, 'm')
    d_column_w, d_line_w = offset(first_cero, d_column, d_line, 'm')

    select_sq = square_selection(r_column_w, r_line_w, second_cero)

    d_column_sq, d_line_sq = square_cut(d_column_w, d_line_w, select_sq)
    r_column_sq, r_line_sq = square_cut(r_column_w, r_line_w, select_sq)

    # Cambiamos a cero de cuadradoselect_UN
    r_column_sq, r_line_sq = offset(second_cero, r_column_sq, r_line_sq, 'm')
    d_column_sq, d_line_sq = offset(second_cero, d_column_sq, d_line_sq, 'm')

    return r_column_sq, r_line_sq, d_column_sq, d_line_sq


def get_data(weighted_image, catalog_fits, example_filter):
    weighted = fits.getdata(weighted_image + '.fits')

    cat = fits.open(catalog_fits + '.fits')
    data = cat[1].data
    real = (data['x'], data['y'])

    header = fits.getheader(example_filter + '.fits', 0)
    wcs_world = WCS(header)
    return weighted, real, wcs_world


def explore_square(names, ceros, fwhm, max_sep):
    weighted, catalog, filter_example = names
    weighted, real, wcs_world = get_data(weighted, catalog, filter_example)

    weighted_cero, square_cero = ceros

    objects, std_dv = detection(weighted, fwhm)
    detec = offset(weighted_cero, objects['x'], objects['y'], sign='p')

    # Convertir a ra/dec usando UNCOVER, cualquier filtro

    detcoords = wcs_world.pixel_to_world(detec[0], detec[1])
    realcoords = wcs_world.pixel_to_world(real[0], real[1])

    line, column = square_cero
    weighted = weighted[column[0]:column[1], line[0]:line[1]]

    catalog_match, under, over, unmatched = catalog_comparison(max_sep, detcoords, realcoords)
    detcoords_under, match_under = under
    detcoords_over, match_over = over
    unmatch_det, unmatch_real = unmatched

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
    ax4.scatter(p, q, color='red', marker='*')

    ax3.set_title('UNCOVER detection unmatched')
    ax3.imshow(weighted, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')
    ax3.scatter(p, q, color='red', marker='*')
    return plt.show()


def explore_and_graph(names, ceros, fwhm, max_sep):
    a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv = explore_square(names, ceros, fwhm, max_sep)
    return graph(a, b, c, d, e, f, g, h, j, k, n, o, p, q, r, t, weighted, std_dv)
