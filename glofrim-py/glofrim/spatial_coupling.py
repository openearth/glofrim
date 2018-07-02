import numpy as np
import rasterio

## coupling functions
def get_rgrid_index(transfrom, x, y):
    # rio transform
    return rows, cols, valid

def get_ugrid_index():
    # KDTree
    return

def get_1d_index():
    # RTREE
    return


def couple_rgrid_to_1d(rgrid_transform, coords_1d, rgrid_mask):
    """
    """
    # get cell indices at 1D coordinates
    rows, cols, valid = get_rgrid_index(rgrid_transform, coords_1d, rgrid_mask)
    i = np.arange(rows.size)

    # create coupled 1-to-1 downstream-to-upstream indices dictionary
    coupled_dict = {i: (r, c) for i, r, c in zip(np.arange(rows.size)[valid], rows[valid], cols[valid])}
    return coupled_dict


def get_coupled_area_frac(coupled_dict, area_other):
    # check valid area other cells
    # get area fractions of coupled cells
    area_frac = {} # [m2/m2]
    for _, idx in list(coupled_dict.items()):
        area_other_sum_idx = np.sum(area_other[idx])
        afs = {i: area_other[i]/area_other_sum_idx for i in idx}
        area_frac.update(afs)
    coupled_area_frac = np.array([area_frac[i] for i in np.arange(area_other.size)])
    return coupled_area_frac

def get_rgrid_mask(coupled_dict, dd, rgrid_shape):
    # get mask with 1) coupled runoff and 2) coupled discharge
    rows, cols = zip(*coupled_dict.keys())
    rows, cols = np.atleast_1d(rows), np.atleast_1d(cols)
    rows_us, cols_us = dd.next_upstream(rows, cols) # find upstream cells where to couple discharge
    # mask of the coupled domain
    coupled_mask = np.zeros(shape)
    coupled_mask[rows_us, cols_us] = 2 # couple discharge to next downstream cell
    coupled_mask[rows, cols] = 1 # couple runoff
    return coupled_mask


def couple_grid_to_grid(self, other):
    """Couple external grid to internal model 2d grid.

    The exact grid coupling method is model specific:
    - CMF: via the inpmat
    Note: only implemented for the CMF model

    Arguments
    ---------
    other : BMI_model_wrapper object
        downstream BMI_model_wrapper object
    """
    if (self.name != 'PCR') or (other.name != 'CMF'):
        msg = 'Grid to grid coupling has only been implemented PCR to CMF (other) coupling'
        raise NotImplementedError(msg)

    # check model grid extends
    if not hasattr(self, 'model_grid_bounds'):
        self.get_model_grid()
    if not hasattr(other, 'model_grid_bounds'):
        other.get_model_grid()
    bounds, res = self.model_grid_bounds, self.model_grid_res

    logger.info('Coupling {} grid to {} grid.'.format(self.name, other.name))
    # set impmat of CMF model based on the model grid bounds and resolution
    if other.name == 'CMF':
        other.set_inpmat_file(bounds, res)

    # deactivate routing in full domain
    self.deactivate_routing('all')
