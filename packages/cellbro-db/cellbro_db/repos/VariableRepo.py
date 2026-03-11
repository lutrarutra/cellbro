import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..types import VariableType, AnnDataLayerType

from ..models import Variable

class VariableRepo:
    @classmethod
    def Get(cls, id: int) -> sa.Select[tuple[Variable]]:
        q =  sa.select(Variable)
        q = q.where(Variable.id == id)
        return q
    
    @classmethod
    def GetByNameAndLayer(cls, name: str, layer: AnnDataLayerType) -> sa.Select[tuple[Variable]]:
        q =  sa.select(Variable)
        q = q.where(Variable.name == name)
        q = q.where(Variable.layer == layer)
        return q
    
    @classmethod
    def Create(cls, name: str, layer: AnnDataLayerType, type: VariableType) -> Variable:
        return Variable(name=name, layer=layer, type=type)
    
    @classmethod
    def Select(
        cls,
        query: sql.Select[tuple[Variable]] = sa.select(Variable),
        name: str | None = None,
        layer: AnnDataLayerType | None = None,
        type: VariableType | None = None,
        layer_in: list[AnnDataLayerType] | None = None,
        type_in: list[VariableType] | None = None,
    ) -> sa.Select[tuple[Variable]]:
        if layer is not None:
            query = query.where(Variable.layer == layer)

        if type is not None:
            query = query.where(Variable.type == type)

        if layer_in is not None:
            query = query.where(Variable.layer.in_(layer_in))

        if type_in is not None:
            query = query.where(Variable.type.in_(type_in))

        if name is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(Variable.name, name.lower()).desc()))
        
        return query