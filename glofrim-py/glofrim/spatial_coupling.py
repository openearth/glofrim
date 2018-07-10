
def couple_routing_2d(rout_bmi, hyd_bmi):
    """"""
    if (rout_bmi._name == 'CaMa-Flood') and (hyd_bmi.get_grid_type() == 2):
        hyd_bmi.get_grid_transform()
        rout_bmi.set_inpmat(hyd_bmi.get_grid_bounds(), hyd_bmi.get_grid_res())
    else:
        raise NotImplementedError('Unknown routing model and/or hydrology model grid type')

def couple_flood_1d(flood_bmi, bmi):
    """"""
    raise NotImplementedError

def couple_flood_2d(flood_bmi, bmi):
    """"""
    raise NotImplementedError