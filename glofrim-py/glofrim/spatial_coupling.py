from gbmi import GBmiGridType


class SpatialCoupling(object):
    def __init__(self):
        self.method = ''
        self.from_name = ''
        self.from_grid_type = GBmiGridType.UNKNOWN
        self.to_name = ''
        self.to_grid_type = GBmiGridType.UNKNOWN
        self.from_ind = []
        self.to_ind = []
        self.area_frac = []

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "SpatialCoupling object from {} to {} based on {} method".format(self.from_name, self.to_name, self.method)

    def to_file(self, filename):
        raise NotImplementedError()

    def from_file(self, filename):
        raise NotImplementedError()

    def _set_to_bmi(self, to_bmi):
        self.to_name = to_bmi._name
        self.to_grid_type = to_bmi.get_grid_type()

    def _set_from_bmi(self, from_bmi):
        self.from_name = from_bmi._name
        self.from_grid_type = from_bmi.get_grid_type()

    def inpmat_2d(self, to_bmi, from_bmi):
        """"""
        self.method = "inpmat_2d"
        self._set_from_bmi(from_bmi)
        self._set_to_bmi(to_bmi)
        if not self.to_grid_type == 5:
            raise NotImplementedError('{} model has incorrect grid type. must by UNITCATCHMENT'.format(to_bmi._name))
        if not self.from_grid_type == 2: # RECTILINEAR
            raise ValueError('{} model has incorrect grid type, must be RECTILEAR'.format(from_bmi._name))
        to_bmi.set_inpmat(from_bmi.get_grid_bounds(), from_bmi.get_grid_res())

    def upstream_1d(self, flood_bmi, bmi):
        """"""
        raise NotImplementedError

    def all_1d(self, flood_bmi, bmi):
        """"""
        raise NotImplementedError

    def couple_flood_2d(self, flood_bmi, bmi):
        """"""
        raise NotImplementedError