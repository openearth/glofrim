
import rasterio.transform
import numpy as np 
import rasterio 
from os.path import isfile, join, dirname
from abc import ABCMeta, abstractmethod

import glofrim.glofrim_lib as glib

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

    def plot_2d(self, ax):
        import matplotlib
        cells = self.get_poly_coords()
        gpatches = (matplotlib.patches.Polygon(cell) for cell in cells)
        gpc = matplotlib.collections.PatchCollection(gpatches, edgecolor='grey')
        gpc.set_facecolor('none')
        ax.add_collection(gpc)
        return gpc

    def plot_1d(self, ax, c=None, **kwargs):
        if self.has_1d():
            idx_us = self._1d.links[:,0]
            idx_ds = self._1d.links[:,1]
            x = self._1d.nodes[idx_us,0]
            y = self._1d.nodes[idx_us,1]
            dx = self._1d.nodes[idx_ds,0]-x
            dy = self._1d.nodes[idx_ds,1]-y
            args = (x, y, dx, dy)
            if c is not None:
                args += (c, )
            kwargs.update(scale_units='xy', angles='xy', scale=1)
            q = ax.quiver(*args, **kwargs)
            return q
        else:
            raise ValueError('grid object doens not have a 1D network')

class RGrid(Grid):
    type = GridType.RGRID
    def __init__(self, transform, height, width, crs=None, mask=None, flip_transform=False, **kwargs):
        """
        If flip_transform is set to True, then the transform will be flipped so that grids read from BMI are compatible with
        the North-South orientation of the grid as read from a clone map with rasterio
        """
        self.transform = rasterio.transform.guard_transform(transform)
        self.height = int(height)
        self.width = int(width)
        self.shape = (self.height, self.width)
        self.bounds = rasterio.transform.array_bounds(self.height, self.width, self.transform)
        self.res = self.transform.a
        self.NtoS = self.transform.e < 0
        self.crs = crs
        self.set_mask(mask)
        if flip_transform:
            self.flip_transform()

    def flip_transform(self):
        """
        Flips the transform upside down so that origin is lower left corner. Use if a typical rasterio.read
        on a map gives the up-side-down of that same map when called from BMI (WFlow has this behaviour)
        """
        bounds = rasterio.transform.array_bounds(self.height, self.width, self.transform)
        west, south, east, north = bounds
        # call the from_bounds method, but with north and south flipped!!
        self.transform = rasterio.transform.from_bounds(west, north, east, south, self.width, self.height)

    def index(self, x, y, flat_index=True, src_crs=None, **kwargs):
        """
        Returns the index values within destination grid of coordinates 'x' and 'y'.
        If a coordinate ref. system is provides, 'x' and 'y' will first be transformed to the projection of self object
        before computing index values

        Parameters
        ----------
        x : float
            x value in coordinate reference system
        y : float
            y value in coordinate reference system
        src_crs : string
            rasterio crs object or string defining the crs of the x, y coordinates
        """
        if (src_crs is not None) and (src_crs is not self.crs):
            # coordinate ref systems of src and dst are different, so a transform is needed
            x, y = rasterio.warp.transform(src_crs, self.crs, x, y)
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

    def _inside_bounds(self, x, y):
        xmin, ymin, xmax, ymax = self.bounds
        inside = np.logical_and.reduce((x >= xmin, x <= xmax, y >= ymin, y <= ymax))
        return inside

    def set_dd(self, fn, ddtype, **kwargs):
        from .nb.nb_io import read_dd_pcraster, read_dd_rasterio, read_dd_cmfbin
        if (ddtype == 'ldd') and fn.endswith('.map'):
            self._dd = read_dd_pcraster(fn, transform=self.transform, **kwargs)
        if (ddtype == 'nextxy') and fn.endswith('.bin'):
            self._dd = read_dd_cmfbin(fn, self.transform, self.height, self.width)
        else:
            self._dd = read_dd_rasterio(fn, ddtype=ddtype, **kwargs)

    def get_pits(self):
        if self._dd is not None:
            # set index of pits (locatoins of outlfow)
            row, col = self._dd.get_pits()
            self.pits = self.ravel_multi_index(row, col)
        else:
            raise ValueError("The grid has no drainage direction defined")

    def set_mask(self, mask):
        if mask is not None:
            assert np.all(mask.shape == self.shape)
        self.mask = mask



class UCGrid(RGrid):
    type = GridType.UCGRID
    def __init__(self, transform, height, width, fn_location_txt, crs=None, mask=None):
        super(UCGrid, self).__init__(transform, height, width, crs=crs, mask=mask)
        self.hires_dir = dirname(fn_location_txt)
        self.hires_locs = self._parse_location(fn_location_txt)
        self.hires_rgrids = [RGrid(p['transform'], p['height'], p['width']) for p in self.hires_locs]
        
    def index(self, x, y, flat_index=True, **kwargs):
        """Get CMF indices (1d) of xy coordinates using the catmxy files"""
        x, y = np.asarray(x), np.asarray(y)
        # find hires catmxy file for x, y coords
        c, r = np.ones(x.size, dtype=int)*-1, np.ones(x.size, dtype=int)*-1
        for i, rg in enumerate(self.hires_rgrids):
            inds_hires = rg.index(x, y)
            valid = inds_hires >= 0
            # read low-res CMF fortran one-based index
            # invalid indices have value - 1
            # go from fortran one-based to python zero-based indices
            fn = join(self.hires_dir, '{}.catmxy'.format(self.hires_locs[i]['area']))
            if not isfile(fn):
                raise IOError('hires catmxy file not found {}'.format(fn))
            a = np.memmap(filename=str(fn), dtype='int16', mode="r+", shape=(2, rg.height, rg.width))
            c[valid] = np.array(a[0, :, :].flat[inds_hires[valid]]-1).astype(int)
            r[valid] = np.array(a[1, :, :].flat[inds_hires[valid]]-1).astype(int)
        # check if inside low res domain
        inside = self._inside(r, c)
        if flat_index:
            # calculate 1d index; NOTE invalid indices have value -1
            inds = np.ones(x.size, dtype=int) * -1 # fill with -1
            if np.any(inside):
                inds[inside] = self.ravel_multi_index(r[inside], c[inside])
        else:
            r[~inside], c[~inside] = -1, -1
            inds = (r, c) 
        return inds

    def _parse_location(self, fn):
        """parse location.txt file in CMF map/hires folder to obtain high-res grid definition"""
        rename = {'ny': 'height',
                'nx': 'width',
                'csize': 'res'}
        with open(fn, 'r') as txt:
            locs = []
            lines = txt.readlines()
            for line in lines:
                vals = line.strip().split()
                assert len(vals) == 2, "invalid location.txt file"
                if vals[0] == 'code':
                    code = int(vals[1])
                    locs.append({})
                elif vals[0] == 'area':
                    locs[code-1]['area'] = str(vals[1])
                else:
                    locs[code-1][rename.get(str(vals[0]), str(vals[0]))] = float(vals[1])
        for p in locs:
            p['transform'] = rasterio.transform.from_origin(p['west'], p['north'], p['res'], p['res'])
        return locs

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
        if links is not None:
            self.links = np.atleast_2d(links)
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
