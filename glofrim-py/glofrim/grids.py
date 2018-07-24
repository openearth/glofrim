
import rasterio.transform
import numpy as np 
import rasterio 
from os.path import isfile, join
from abc import ABCMeta, abstractmethod

import glofrim_lib as glib

class GridType(object):
    """
    Enumeration with grid types.
    """
    UNKNOWN = 0
    RGRID = 1
    UGRID = 2
    UCGRID = 3 # Unit-catchment grid

class Grid(object):
    __metaclass__ = ABCMeta
    type = GridType.UNKNOWN
    _1d = None  # placehorder for 1d network
    _dd = None  # placeholder for drainage direction object

    @abstractmethod
    def index(self, x, y):
        """
        Returns the index of the grid cell containing (x, y) given a
        coordinate reference system. In If (x, y) fall outside the domain an
        index value -1 is returned

        Parameters
        ----------
        x : float
            x value in coordinate reference system
        y : float
            y value in coordinate reference system
        Returns
        -------
        numpy array
            index
        """
        pass

    def get_type(self):
        """
        Returns the GridType
        """
        return self.type

    def set_1d(self, nodes, links=None, inds=None, fill_value=-1, ds_idx=None, us_idx=None):
        """
        Sets the 1d network attribute
        """
        self._1d = Network1D(nodes, links=links, inds=inds, fill_value=fill_value, ds_idx=ds_idx, us_idx=us_idx)

    @abstractmethod
    def set_dd(self, fn, **kwargs):
        """
        Sets the drainage direction map attribute. Only implemented for RGrid and UCGrid
        """
        pass

    def has_dd(self):
        return self._dd is not None

    def has_1d(self):
        return self._1d is not None

class RGrid(Grid):
    type = GridType.RGRID
    def __init__(self, transform, height, width, crs=None, mask=None):
        self.transform = rasterio.transform.guard_transform(transform)
        self.height = int(height)
        self.width = int(width)
        self.shape = (self.height, self.width)
        self.bounds = rasterio.transform.array_bounds(self.height, self.width, self.transform)
        self.res = self.transform.a
        self.NtoS = self.transform.e < 0
        self.crs = crs
        self.mask = mask

    def index(self, x, y, flat_index=True, **kwargs):
        x, y = np.atleast_1d(x), np.atleast_1d(y)
        r, c = rasterio.transform.rowcol(self.transform, x, y, **kwargs)
        r, c = np.asarray(r).astype(int), np.asarray(c).astype(int)
        inside = self._inside(r, c)
        if flat_index:
            # calculate 1d index; NOTE invalid indices have value -1
            inds = np.ones_like(inside, dtype=int)
            inds = inds * -1 # fill with -1
            if np.any(inside):
                inds[inside] = self.ravel_multi_index(r[inside], c[inside])
        else:
            r[~inside], c[~inside] = -1, -1
            inds = (r, c) 
        return inds

    def xy(self, ind=None, row=None, col=None, **kwargs):
        if ind is not None:
            ind = np.atleast_1d(ind)
            row, col = self.unravel_index(ind)
        if (row is None) or (col is None):
            raise ValueError('Provide either row and col or flat index')
        row, col = np.atleast_1d(row), np.atleast_1d(col)
        inside = self._inside(row, col)
        x, y = np.ones(row.size)*-1, np.ones(row.size)*-1
        x[inside], y[inside] = rasterio.transform.xy(self.transform, row[inside], col[inside])
        return x, y

    def get_poly_coords(self):
        rows, cols = np.where(self.mask)
        xs, ys = self.xy(row=rows, col=cols)
        r = self.res/2
        bounds = [[x-r, y-r, x+r, y+r] for x,y in zip(xs, ys)]
        cells = [[[s,e], [n,e], [n,w], [s,w], [s,e]] for s, e, n, w in bounds]
        return cells

    def ravel_multi_index(self, r, c):
        return np.ravel_multi_index((r, c), (self.height, self.width))

    def unravel_index(self, flat_ind):
        return np.unravel_index(flat_ind, (self.height, self.width))

    def _inside(self, r, c):
        inside = np.logical_and.reduce((r>=0, r<self.height, c>=0, c<self.width))
        if self.mask is not None:
            inside[inside] = np.logical_and(inside[inside], self.mask[r[inside], c[inside]])
        return inside

    def set_dd(self, fn, ddtype, **kwargs):
        from nb.nb_io import read_dd_pcraster, read_dd_rasterio
        if (ddtype == 'ldd') and fn.endswith('.map'):
            self._dd = read_dd_pcraster(fn, transform=self.transform, **kwargs)
        else:
            self._dd = read_dd_rasterio(fn, ddtype=ddtype, **kwargs)



class UCGrid(RGrid):
    type = GridType.UCGRID
    def __init__(self, transform, height, width, fn_catmxy, crs=None):
        super(UCGrid, self).__init__(transform, height, width, crs=crs)
        self.fn_catmxy = fn_catmxy
        
    def index(self, x, y, flat_index=True, **kwargs):
        """Get CMF indices (1d) of xy coordinates using the catmxy.tif file"""
        # read catmxy temporary into memory
        with rasterio.open(self.fn_catmxy, 'r') as ds:
            if ds.count != 2:
                raise ValueError("{} file should have two layers".format(self.fn_catmxy))
            catmxy_rg = RGrid(ds.transform, ds.height, ds.width, ds.crs)
            # python zero based index for high res CMF grid
            inds = catmxy_rg.index(x, y)
            valid = inds >= 0
            # read low-res CMF fortran one-based index
            # invalid indices have value - 1
            # go from fortran one-based to python zero-based indices
            c, r = np.ones_like(inds)*-1, np.ones_like(inds)*-1
            c[valid] = ds.read(1).flat[inds[valid]] - 1
            r[valid] = ds.read(2).flat[inds[valid]] - 1
        # check if inside low res domain
        inside = self._inside(r, c)
        if flat_index:
            # calculate 1d index; NOTE invalid indices have value -1
            inds = np.ones_like(inside, dtype=int)
            inds = inds * -1 # fill with -1
            if np.any(inside):
                inds[inside] = self.ravel_multi_index(r[inside], c[inside])
        else:
            r[~inside], c[~inside] = -1, -1
            inds = (r, c) 
        return inds

    def set_inpmat(self, bounds, res, mapdir, olat='NtoS'):
        """Set the CMF inpmat file model based on the grid definition of upstream model"""
        if not isfile(join(mapdir, 'generate_inpmat')):
            raise ValueError('{} not found'.format(join(mapdir, 'generate_inpmat')))
        w, s, e, n = bounds
        # generate inpmat
        cmd = './generate_inpmat {} {} {} {} {} {:s}'
        cmd = cmd.format(abs(res), w, e, n, s, olat)
        # print(cmd)
        glib.subcall(cmd, cwd=mapdir)

class UGrid(Grid):
    type = GridType.UGRID

    # TODO: connect to pyugrid package at https://github.com/pyugrid/pyugrid/blob/master/pyugrid/ugrid.py

    def __init__(self, nodes, faces, 
                edges=None, boundaries=None, face_face_connectivity=None, 
                fill_value=-1, inds=None, crs=None,
                face_coordinates=None, edge_coordinates=None, boundary_coordinates=None):
        self.nodes = np.atleast_2d(nodes)
        self.faces = np.ma.masked_equal(np.atleast_2d(faces), fill_value)
        self.nnodes = self.nodes.shape[0]
        self.nfaces = self.faces.shape[0]
        assert np.logical_and(self.faces.max() <= self.nnodes, self.faces.min() >= 0)
        # optional index array
        self.inds = np.arange(self.nfaces) if inds is None else np.array(inds)
        assert self.nfaces == self.inds.size
        # optional links: links[0,:]=from-idx, links[1,:]=to-idx
        self.edges = edges
        self.boundaries = boundaries 
        self.face_face_connectivity = face_face_connectivity
        # optional coordinate array
        self.edge_coordinates = edge_coordinates
        self.face_coordinates = face_coordinates
        self.boundary_coordinates = boundary_coordinates
        # coordinate reference system
        self.crs = crs
        self._tree = None
    
    def get_poly_coords(self):
        n_nodes_per_face = (~self.faces.mask).sum(axis=1)
        ragged = [
            face[:n_nodes].filled()
            for n_nodes, face 
            in zip(n_nodes_per_face, self.faces)
        ]
        fcoords = [self.nodes[np.append(face, face[0])] for face in ragged]
        return fcoords

    def index(self, x, y, **kwargs):
        raise NotImplementedError()
        # if self._tree is None:
        #     self._build_celltree()
        # indices = self._tree.locate(points)
        # return indices

    def build_celltree(self):
        """
        Tries to build the celltree for the current UGrid. Will fail if nodes
        or faces is not defined.
        """
        try:
            from cell_tree2d import CellTree
        except ImportError:
            raise ImportError("the cell_tree2d package must be installed to "
                                "use the celltree search:\n"
                                "https://github.com/NOAA-ORR-ERD/cell_tree2d/")
        self._tree = CellTree(self.nodes, self.faces)

    def set_dd(self):
        raise NotImplementedError("No drainage direction method is implemented for UGrid")

class Network1D(object):
    links = None
    nlinks = 0
    def __init__(self, nodes, links=None, inds=None, fill_value=-1, ds_idx=None, us_idx=None):
        self.nodes = np.atleast_2d(nodes)
        self.nnodes = self.nodes.shape[0]
        # optional links: links[0,:]=from-idx, links[1,:]=to-idx
        if self.links is not None:
            self.links = np.ma.masked_equal(np.atleast_2d(self.nodes), fill_value)
            self.nlinks = self.links.shape[0]
        # optional index array 
        # indices only used for index function only!
        self.inds = np.arange(self.nnodes) if inds is None else np.array(inds)
        self.ds_idx = np.array([]) if ds_idx is None else np.array(ds_idx) # indices of downstream nodes
        self.us_idx = np.array([]) if us_idx is None else np.array(us_idx) # indices of upstream nodes
        assert self.nnodes == self.inds.size
        # will be created later
        self._tree = None

    def get_line_coords(self):
        if self.links is None:
            raise ValueError('links between nodes not provided')
        lcoords = self.nodes[self.links]
        return lcoords

    def get_end_points(self):
        if self.links is None:
            raise ValueError('links between nodes not provided')
        idx, counts = np.unique(self.links, return_counts=True)
        return idx[counts==1]
        
    def get_downstream_points(self):
        end_idx = self.get_end_points()
        self.ds_idx = np.array([p for p in end_idx if p in self.links[0, :]])
        return self.ds_idx

    def get_upstream_points(self):
        end_idx = self.get_end_points()
        self.us_idx = np.array([p for p in end_idx if p in self.links[1, :]])
        return self.us_idx

    def index(self, x, y):
        if self._tree is None:
            self._build_rtree()
        x, y = np.atleast_1d(x), np.atleast_1d(y)
        return np.array([list(self._tree.nearest(xy, 1))[0] for xy in zip(x, y)])

    def xy(self, ind):
        orig_indices = self.inds.copy().argsort()
        ndx = orig_indices[np.searchsorted(self.inds[orig_indices], ind)]
        return self.nodes[ndx].squeeze()

    def _build_rtree(self):
        """Creat a spatial index for the 1d coordinates. A model_1d_index
        attribute funtion is created to find the nearest 1d coordinate tuple"""
        # build spatial rtree index of points
        import rtree
        self._tree = rtree.index.Index()
        for i, xy in zip(self.inds, self.nodes):
            self._tree.insert(i, xy.tolist())
