from typing import Optional, Union

from dataclasses import dataclass
from abc import ABC, abstractmethod

import pickle, redis

import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import scipy


def get_len(idx: Union[int, slice, str, list[str], list[int]], max_len: int) -> int:
    if isinstance(idx, int) or isinstance(idx, str):
        return 1

    if isinstance(idx, list):
        return len(idx)

    _obs_start = 0 if idx.start is None else idx.start
    _obs_stop = max_len if idx.stop is None else idx.stop

    return _obs_stop - _obs_start

@dataclass
class RedisData:
    _key: str
    _redis: redis.client.Redis

    def __init__(self, key: str, _redis: redis.client.Redis, data: any = None):
        self._key = key
        self._redis = _redis
        self.data = data

    def __del__(self):
        self._redis.delete(self._key)

    def __repr__(self) -> str:
        return f"RedisData(key={self._key})"

    @property
    def key(self):
        return self._key

    @property
    def data(self):
        return self.get()

    @data.setter
    def data(self, data: any):
        assert self._redis.set(self._key, self.dumps(data)), "Failed to set data"

    def get(self, default=None):
        bts = self._redis.get(self._key)
        if bts:
            return self.loads(bts)
        else:
            return default

    def copy(self):
        return self.get()

    def __getitem__(self, key):
        return self.get()[key]

    def dumps(self, var: any) -> bytes:
        _type = type(var)

        if _type is scipy.sparse.csr_matrix:
            assert var.data.nbytes < 512 * 1024 * 1024, "Data is too large"

        if _type is np.ndarray:
            if var.nbytes >= 512 * 1024 * 1024:
                return pickle.dumps(scipy.sparce.csr_matrix(var))

        if _type in [scipy.sparse.csr_matrix, np.ndarray, pd.DataFrame, pd.Series, dict, list, tuple, set, frozenset, int, float, str, bool]:
            return pickle.dumps(var)

        if _type in [anndata.compat.OverloadedDict]:
            return pickle.dumps(dict(var))

        assert False, _type

    def loads(self, value):
        return pickle.loads(value)


@dataclass
class RedisMatrix(RedisData):
    _shape: tuple[int, int] = None

    def __init__(self, key: str, _redis: redis.client.Redis, data: Union[np.ndarray, scipy.sparse.csr_matrix]):
        super().__init__(key, _redis, data)
        self._shape = data.shape
    
    def __repr__(self) -> str:
        return f"RedisMatrix(key={self._key}, shape={self._shape})"

    @property
    def shape(self):
        return self._shape

    def __getattribute__(self, __name: str):
        try:
            atr = super(RedisMatrix, self).__getattribute__(__name)
        except:
            atr = super(np.ndarray, super(RedisMatrix, self).__getattribute__(
                "data")).__getattribute__(__name)

        return atr

    
@dataclass
class RedisFrame(RedisData):
    def __init__(self, key: str, _redis: redis.client.Redis, data: pd.DataFrame):
        super().__init__(key, _redis, data)

    def __repr__(self) -> str:
        return f"RedisFrame(key={self._key})"

    def __getattribute__(self, __name: str):
        try:
            atr = super(RedisFrame, self).__getattribute__(__name)
        except:
            atr = super(pd.DataFrame, super(RedisFrame, self).__getattribute__(
                "data")).__getattribute__(__name)

        return atr

    

@dataclass    
class RedisTable:
    _data: dict[str, RedisData] = None
    _redis: redis.client.Redis = None
    _keys: list[str] = None
    _redis_key_prefix: str = ""

    def __init__(self, data: dict, _redis: redis.client.Redis, redis_key_prefix: str = ""):
        self._redis = _redis
        self._keys = list(data.keys())
        self._redis_key_prefix = redis_key_prefix

        self._data = {}

        for key in self._keys:
            if type(data[key]) in [np.ndarray, scipy.sparse.csr_matrix]:
                self._data[key] = RedisMatrix(f"{self._redis_key_prefix}_{key}", self._redis, data[key])
            else:
                self._data[key] = RedisData(f"{self._redis_key_prefix}_{key}", self._redis, data[key])

    def __del__(self):
        for key in self._keys:
            del self._data[key]
        
        # del _data
        # del _keys

    def __repr__(self) -> str:
        return f"RedisTable(keys={self._keys}, redis_key_prefix={self._redis_key_prefix})"

    def __str__(self) -> str:
        return self.__repr__()

    def __getitem__(self, key):
        return self._data[key].get()

    def __setitem__(self, key, value):
        if key not in self._keys:
            self._keys.append(key)
            if type(value) in [np.ndarray, scipy.sparse.csr_matrix]:
                self._data[key] = RedisMatrix(f"{self._redis_key_prefix}_{key}", self._redis, value)
            else:
                self._data[key] = RedisData(f"{self._redis_key_prefix}_{key}", self._redis, value)
        else:
            self._data[key].data = value

    def keys(self):
        return self._keys

    def values(self):
        return [self._data[key].get() for key in self._keys]

    def items(self):
        return [(key, self._data[key].get()) for key in self._keys]

    def get(self, key=None, default=None):
        if key is not None:
            return {key: self._data[key].get(default=default) for key in self._keys}

        return self._data[key].get()


@dataclass
class AnnDatabase(anndata.AnnData):
    host: str = "localhost"
    port: int = 6379
    db: int = 0
    _redis: redis.client.Redis = None
    _is_view: bool = False

    _shape: tuple[int, int] = None
    _n_obs: int = None
    _n_vars: int = None

    _X: RedisMatrix = None
    _obs: RedisFrame = None
    _var: RedisFrame = None
    
    _var_keys: RedisData = None
    _obs_keys: RedisData = None
    _var_names: RedisData = None
    _obs_names: RedisData = None

    _obsm_keys: RedisData = None
    _layers_keys: RedisData = None
    _varm_keys: RedisData = None
    _uns_keys: RedisData = None
    _obsp_keys: RedisData = None
    _varp_keys: RedisData = None

    _uns: RedisTable = None
    _obsm: RedisTable = None
    _layers: RedisTable = None
    _varm: RedisTable = None

    _obsp: RedisData = None
    _varp: RedisData = None

    def __post_init__(self):
        # Start Server
        self._redis = redis.Redis(host=self.host, port=self.port, db=self.db)
        
        # self._raw = self
        # self.file = None

    def init(self, adata: anndata.AnnData):
        #
        self._shape = adata.shape
        self._n_obs = adata.n_obs
        self._n_vars = adata.n_vars

        # self.clean()

        self._X = RedisMatrix("X", self._redis, adata.X)

        self._obs = RedisFrame("obs", self._redis, adata.obs)
        self._obs_keys = RedisData("obs_keys", self._redis, list(adata.obs.keys()))
        self._obs_names = RedisData("obs_names", self._redis, list(adata.obs_names))

        self._var = RedisFrame("var", self._redis, adata.var)
        self._var_keys = RedisData("var_keys", self._redis, list(adata.var.keys()))
        self._var_names = RedisData("var_names", self._redis, list(adata.var_names))

        self._varm_keys = RedisData("varm_keys", self._redis, list(adata.varm.keys()))
        self._varm = RedisTable(adata.varm, self._redis, "varm")

        self._obsm_keys = RedisData("obsm_keys", self._redis, list(adata.obsm.keys()))
        self._obsm = RedisTable(adata.obsm, self._redis, "obsm")

        self._layers_keys = RedisData("layers_keys", self._redis, list(adata.layers.keys()))
        self._layers = RedisTable(adata.layers, self._redis, "layers")

        self._uns_keys = RedisData("uns_keys", self._redis, list(adata.uns.keys()))
        self._uns = RedisTable(adata.uns, self._redis, "uns")
    
        self._obsp_keys = RedisData("obsp_keys", self._redis, list(adata.obsp.keys()))
        self._obsp = RedisTable(adata.obsp, self._redis, "obsp")

        self._varp_keys = RedisData("varp_keys", self._redis, list(adata.varp.keys()))
        self._varp = RedisTable(adata.varp, self._redis, "varp")

    def clean(self):
        self._X = None
        self._obs = None
        self._var = None
        self._uns = None
        self._obsm = None
        self._layers = None
        self._varm = None
        self._obsp = None
        self._varp = None

        self._uns = None
        self._var_keys = None
        self._obs_keys = None
        self._var_names = None
        self._obs_names = None
        self._obsm_keys = None
        self._layers_keys = None
        self._varm_keys = None
        self._uns_keys = None
        self._obsp_keys = None
        self._varp_keys = None

    def __del__(self):
        self.clean()

    def to_anndata(self):
        adata = anndata.AnnData(X=self.X.get(), obs=self.obs.get(), var=self.var.get())
        adata.uns = self.uns.get()
        adata.obsm = self.obsm.get()
        adata.layers = self.layers.get()
        adata.varm = self.varm.get()
        adata.obsp = self.obsp.get()
        adata.varp = self.varp.get()
        return adata

    def __getitem__(self, key):
        assert type(key) in [int, slice, str, list, tuple], "key must be a str, slice, int, or a tuple[int/str]"

        if not isinstance(key, tuple):
            obs_idx = key
            var_idx = slice(None, None, None)
        else:
            assert len(key) == 2, "key must be a slice or a tuple[tuple] of length 2"
            obs_idx, var_idx = key

        if type(obs_idx) == str:
            obs_idx = [obs_idx]

        if isinstance(obs_idx, list) and issubclass(obs_idx, list[str]):
            _temp = []
            for idx in obs_idx:
                if isinstance(idx, str):
                    _temp.append(self.obs_names.index(idx))
                else:
                    _temp.append(idx)
        
            obs_idx = _temp
        
        if type(var_idx) == str:
            var_idx = [var_idx]
        
        if isinstance(var_idx, list):
            _temp = []
            for idx in var_idx:
                if isinstance(idx, str):
                    _temp.append(self.var_names.index(idx))
                else:
                    _temp.append(idx)
        
            var_idx = _temp

        return AnnDatabaseView(self, obs_idx, var_idx)
        
    def read(self, path: str):
        # load anndata
        adata = sc.read_h5ad(path)
        self.init(adata)
        del adata

    def obs_keys(self) -> list[str]:
        return self._obs_keys.get()

    def obsm_keys(self) -> list[str]:
        return self.obsm_keys.get()

    def var_keys(self) -> list[str]:
        return self._var_keys.get()

    def varm_keys(self) -> list[str]:
        return self._varm_keys.get()

    @property
    def n_obs(self) -> int:
        return self._n_obs
    
    @property
    def n_vars(self) -> int:
        return self._n_vars
    
    @property
    def shape(self) -> tuple[int, int]:
        return self._shape

    @property
    def var_names(self) -> list[str]:
        return self._var_names.get()
    
    @property
    def obs_names(self) -> list[str]:
        return self._obs_names.get()

    @property
    def X(self) -> RedisData:
        return self._X

    @property
    def obs(self) -> RedisData:
        return self._obs
    
    @property
    def var(self) -> RedisData:
        return self._var

    @property
    def uns(self) -> RedisTable:
        return self._uns

    @property
    def obsm(self) -> RedisTable:
        return self._obsm

    @property
    def layers(self) -> RedisTable:
        return self._layers

    @property
    def varm(self) -> RedisTable:
        return self._varm

    @property
    def obsp(self) -> RedisData:
        return self._obsp

    @property
    def varp(self) -> RedisData:
        return self._varp


@dataclass
class AnnDatabaseView:
    _abase: AnnDatabase
    _obs_idx: Union[int, slice, str, list[str], list[int]]
    _var_idx: Union[int, slice, str, list[str], list[int]]

    def obs_keys(self) -> list[str]:
        return self._abase.obs_keys[self._obs_idx]

    def var_keys(self) -> list[str]:
        return self._abase.var_keys[self._var_idx]

    @property
    def var_names(self) -> list[str]:
        return self._abase.var_names[self._var_idx]

    @property
    def obs_names(self) -> list[str]:
        return self._abase.obs_names[self._obs_idx]

    @property
    def shape(self) -> tuple[int, int]:
        if isinstance(self._obs_idx, slice):
            assert self._obs_idx.step is None, "step is not supported"
        
        if isinstance(self._var_idx, slice):
            assert self._var_idx.step is None, "step is not supported"

        _obs_dim = get_len(self._obs_idx)
        _var_dim = get_len(self._var_idx)
        return (_obs_dim, _var_dim)

    def obs_keys(self) -> list[str]:
        return self._abase._obs_keys.get()[self._obs_idx]

    def obsm_keys(self) -> list[str]:
        return self._abase.obsm_keys.get()[self._obs_idx]

    def var_keys(self) -> list[str]:
        return self._abase._var_keys.get()[self._var_idx]

    def varm_keys(self) -> list[str]:
        return self._abase._varm_keys.get()[self._var_idx]

    @property
    def n_obs(self) -> int:
        return self._abase.shape[0]

    @property
    def n_vars(self) -> int:
        return self._abase.shape[1]

    @property
    def var_names(self) -> list[str]:
        return self._abase._var_names.get()[self._var_idx]

    @property
    def obs_names(self) -> list[str]:
        return self._abase._obs_names.get()[self._obs_idx]

    @property
    def X(self) -> np.ndarray:
        return self._abase._X.get()[self._obs_idx, self._var_idx]

    @property
    def obs(self) -> pd.DataFrame:
        return self._abase._obs.get()[self._obs_idx]

    @property
    def var(self) -> pd.DataFrame:
        return self._abase._var.get()[self._var_idx]

    @property
    def uns(self) -> dict[str, any]:
        assert False, "Not implemented views for uns"
        return self._abase._uns.get()

    @property
    def obsm(self) -> dict[str, any]:
        return self._abase._obsm.get()[self._obs_idx]

    @property
    def layers(self) -> dict[str, any]:
        return self._abase._layers.get()[self._obs_idx, self._var_idx]

    @property
    def varm(self) -> dict[str, any]:
        return self._abase._varm.get()[self._var_idx]

    @property
    def obsp(self) -> dict[str, any]:
        return self._abase._obsp.get()[self._obs_idx, self._obs_idx]

    @property
    def varp(self) -> dict[str, np.ndarray]:
        return self._abase._varp.get()[self._var_idx, self._var_idx]
    

